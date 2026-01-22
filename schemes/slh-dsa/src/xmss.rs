//! XMSS (eXtended Merkle Signature Scheme) implementation for SLH-DSA.
//!
//! XMSS provides a Merkle tree structure over WOTS+ public keys, allowing
//! multiple messages to be signed with a single XMSS key pair. Each leaf
//! of the tree is a WOTS+ public key.
//!
//! # Structure
//!
//! An XMSS tree of height h' contains:
//! - 2^h' WOTS+ key pairs at the leaves
//! - Internal nodes computed by hashing pairs of children
//! - A single root node that serves as the XMSS public key
//!
//! # Signatures
//!
//! An XMSS signature consists of:
//! - The WOTS+ signature on the message
//! - An authentication path (sibling nodes from leaf to root)

use crate::address::{Address, AddressType};
use crate::hash::h;
use crate::params::Params;
use crate::wots::{wots_pk_from_sig, wots_pk_gen, wots_sign};

/// XMSS signature containing WOTS+ signature and authentication path.
#[derive(Clone, Debug)]
pub struct XmssSig {
    /// WOTS+ signature (wots_len n-byte strings)
    pub wots_sig: Vec<Vec<u8>>,
    /// Authentication path (h' n-byte strings)
    pub auth: Vec<Vec<u8>>,
}

impl XmssSig {
    /// Creates a new empty XMSS signature structure.
    pub fn new(params: &Params) -> Self {
        XmssSig {
            wots_sig: Vec::with_capacity(params.wots_len),
            auth: Vec::with_capacity(params.hp),
        }
    }
}

/// Computes an XMSS tree node.
///
/// xmss_node(SK.seed, PK.seed, i, z, ADRS) - FIPS 205 Algorithm 7
///
/// Computes the node at height z and position i in the XMSS tree.
///
/// # Arguments
/// * `sk_seed` - Secret seed
/// * `pk_seed` - Public seed
/// * `i` - Position of the node at height z
/// * `z` - Height of the node (0 = leaf level)
/// * `adrs` - Address structure
/// * `params` - Parameter set
pub fn xmss_node(
    sk_seed: &[u8],
    pk_seed: &[u8],
    i: u32,
    z: u32,
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    if z == 0 {
        // Leaf node: compute WOTS+ public key
        adrs.set_type(AddressType::WotsHash);
        adrs.set_keypair(i);
        return wots_pk_gen(sk_seed, pk_seed, adrs, params);
    }

    // Internal node: hash the two children
    let left = xmss_node(sk_seed, pk_seed, 2 * i, z - 1, adrs, params);
    let right = xmss_node(sk_seed, pk_seed, 2 * i + 1, z - 1, adrs, params);

    adrs.set_type(AddressType::Tree);
    adrs.set_tree_height(z);
    adrs.set_tree_index(i);

    h(pk_seed, adrs, &left, &right, params)
}

/// Generates an XMSS signature.
///
/// xmss_sign(M, SK.seed, idx, PK.seed, ADRS) - FIPS 205 Algorithm 8
///
/// Signs message M using the WOTS+ key pair at leaf index idx.
///
/// # Arguments
/// * `msg` - Message to sign (n bytes)
/// * `sk_seed` - Secret seed
/// * `idx` - Index of the leaf to use
/// * `pk_seed` - Public seed
/// * `adrs` - Address structure
/// * `params` - Parameter set
pub fn xmss_sign(
    msg: &[u8],
    sk_seed: &[u8],
    idx: u32,
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> XmssSig {
    // Generate authentication path
    let auth = compute_auth_path(sk_seed, pk_seed, idx, adrs, params);

    // Generate WOTS+ signature
    adrs.set_type(AddressType::WotsHash);
    adrs.set_keypair(idx);
    let wots_sig = wots_sign(msg, sk_seed, pk_seed, adrs, params);

    XmssSig { wots_sig, auth }
}

/// Computes the authentication path for a leaf.
///
/// The authentication path contains the sibling nodes needed to
/// recompute the root from a leaf.
fn compute_auth_path(
    sk_seed: &[u8],
    pk_seed: &[u8],
    idx: u32,
    adrs: &mut Address,
    params: &Params,
) -> Vec<Vec<u8>> {
    let mut auth = Vec::with_capacity(params.hp);

    for j in 0..params.hp {
        // At height j, we need the sibling of the path node
        // If bit j of idx is 0, sibling is at position idx + 2^j
        // If bit j of idx is 1, sibling is at position idx - 2^j
        let k = (idx >> j) ^ 1; // Sibling index at height j
        let sibling = xmss_node(sk_seed, pk_seed, k, j as u32, adrs, params);
        auth.push(sibling);
    }

    auth
}

/// Computes XMSS public key from signature.
///
/// xmss_PKFromSig(idx, SIG_XMSS, M, PK.seed, ADRS) - FIPS 205 Algorithm 9
///
/// Reconstructs the XMSS root from a signature and message.
///
/// # Arguments
/// * `idx` - Index of the signing leaf
/// * `sig` - XMSS signature
/// * `msg` - Signed message (n bytes)
/// * `pk_seed` - Public seed
/// * `adrs` - Address structure
/// * `params` - Parameter set
pub fn xmss_pk_from_sig(
    idx: u32,
    sig: &XmssSig,
    msg: &[u8],
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    // Compute WOTS+ public key from signature
    adrs.set_type(AddressType::WotsHash);
    adrs.set_keypair(idx);
    let mut node = wots_pk_from_sig(&sig.wots_sig, msg, pk_seed, adrs, params);

    // Climb the tree using authentication path
    adrs.set_type(AddressType::Tree);

    for k in 0..params.hp {
        adrs.set_tree_height((k + 1) as u32);

        // Determine if current node is left or right child
        if (idx >> k) & 1 == 0 {
            // Node is left child, auth[k] is right sibling
            adrs.set_tree_index((idx >> (k + 1)) as u32);
            node = h(pk_seed, adrs, &node, &sig.auth[k], params);
        } else {
            // Node is right child, auth[k] is left sibling
            adrs.set_tree_index((idx >> (k + 1)) as u32);
            node = h(pk_seed, adrs, &sig.auth[k], &node, params);
        }
    }

    node
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::SLH_DSA_SHAKE_128F; // Use 128f for faster tests (hp=3)

    #[test]
    fn test_xmss_node_leaf() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let mut adrs = Address::new();
        adrs.set_layer(0);
        adrs.set_tree(0);

        // Compute leaf 0
        let leaf = xmss_node(&sk_seed, &pk_seed, 0, 0, &mut adrs, &params);
        assert_eq!(leaf.len(), params.n);
    }

    #[test]
    fn test_xmss_node_internal() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let mut adrs = Address::new();
        adrs.set_layer(0);
        adrs.set_tree(0);

        // Compute internal node at height 1
        let node = xmss_node(&sk_seed, &pk_seed, 0, 1, &mut adrs, &params);
        assert_eq!(node.len(), params.n);
    }

    #[test]
    fn test_xmss_node_deterministic() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];

        let mut adrs1 = Address::new();
        adrs1.set_layer(0);
        adrs1.set_tree(0);

        let mut adrs2 = Address::new();
        adrs2.set_layer(0);
        adrs2.set_tree(0);

        let node1 = xmss_node(&sk_seed, &pk_seed, 0, 1, &mut adrs1, &params);
        let node2 = xmss_node(&sk_seed, &pk_seed, 0, 1, &mut adrs2, &params);

        assert_eq!(node1, node2);
    }

    #[test]
    fn test_xmss_root() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let mut adrs = Address::new();
        adrs.set_layer(0);
        adrs.set_tree(0);

        // Compute the root of the XMSS tree
        let root = xmss_node(&sk_seed, &pk_seed, 0, params.hp as u32, &mut adrs, &params);
        assert_eq!(root.len(), params.n);
    }

    #[test]
    fn test_xmss_sign_structure() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];
        let mut adrs = Address::new();
        adrs.set_layer(0);
        adrs.set_tree(0);

        let sig = xmss_sign(&msg, &sk_seed, 0, &pk_seed, &mut adrs, &params);

        // Check WOTS+ signature size
        assert_eq!(sig.wots_sig.len(), params.wots_len);
        for s in &sig.wots_sig {
            assert_eq!(s.len(), params.n);
        }

        // Check authentication path size
        assert_eq!(sig.auth.len(), params.hp);
        for a in &sig.auth {
            assert_eq!(a.len(), params.n);
        }
    }

    #[test]
    fn test_xmss_sign_verify_roundtrip() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        // Compute root
        let mut root_adrs = Address::new();
        root_adrs.set_layer(0);
        root_adrs.set_tree(0);
        let root = xmss_node(&sk_seed, &pk_seed, 0, params.hp as u32, &mut root_adrs, &params);

        // Sign at leaf 0
        let mut sign_adrs = Address::new();
        sign_adrs.set_layer(0);
        sign_adrs.set_tree(0);
        let sig = xmss_sign(&msg, &sk_seed, 0, &pk_seed, &mut sign_adrs, &params);

        // Verify
        let mut verify_adrs = Address::new();
        verify_adrs.set_layer(0);
        verify_adrs.set_tree(0);
        let recovered_root = xmss_pk_from_sig(0, &sig, &msg, &pk_seed, &mut verify_adrs, &params);

        assert_eq!(root, recovered_root);
    }

    #[test]
    fn test_xmss_sign_different_leaves() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        // Compute root
        let mut root_adrs = Address::new();
        root_adrs.set_layer(0);
        root_adrs.set_tree(0);
        let root = xmss_node(&sk_seed, &pk_seed, 0, params.hp as u32, &mut root_adrs, &params);

        // Sign and verify with leaf 0
        let mut adrs0 = Address::new();
        adrs0.set_layer(0);
        adrs0.set_tree(0);
        let sig0 = xmss_sign(&msg, &sk_seed, 0, &pk_seed, &mut adrs0, &params);

        let mut verify0 = Address::new();
        verify0.set_layer(0);
        verify0.set_tree(0);
        let root0 = xmss_pk_from_sig(0, &sig0, &msg, &pk_seed, &mut verify0, &params);
        assert_eq!(root, root0);

        // Sign and verify with leaf 1
        let mut adrs1 = Address::new();
        adrs1.set_layer(0);
        adrs1.set_tree(0);
        let sig1 = xmss_sign(&msg, &sk_seed, 1, &pk_seed, &mut adrs1, &params);

        let mut verify1 = Address::new();
        verify1.set_layer(0);
        verify1.set_tree(0);
        let root1 = xmss_pk_from_sig(1, &sig1, &msg, &pk_seed, &mut verify1, &params);
        assert_eq!(root, root1);

        // Signatures should be different
        assert_ne!(sig0.wots_sig, sig1.wots_sig);
        assert_ne!(sig0.auth, sig1.auth);
    }

    #[test]
    fn test_xmss_wrong_message() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];
        let wrong_msg = vec![4u8; params.n];

        // Compute root
        let mut root_adrs = Address::new();
        root_adrs.set_layer(0);
        root_adrs.set_tree(0);
        let root = xmss_node(&sk_seed, &pk_seed, 0, params.hp as u32, &mut root_adrs, &params);

        // Sign correct message
        let mut sign_adrs = Address::new();
        sign_adrs.set_layer(0);
        sign_adrs.set_tree(0);
        let sig = xmss_sign(&msg, &sk_seed, 0, &pk_seed, &mut sign_adrs, &params);

        // Verify with wrong message
        let mut verify_adrs = Address::new();
        verify_adrs.set_layer(0);
        verify_adrs.set_tree(0);
        let bad_root = xmss_pk_from_sig(0, &sig, &wrong_msg, &pk_seed, &mut verify_adrs, &params);

        assert_ne!(root, bad_root);
    }

    #[test]
    fn test_xmss_wrong_index() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        // Compute root
        let mut root_adrs = Address::new();
        root_adrs.set_layer(0);
        root_adrs.set_tree(0);
        let root = xmss_node(&sk_seed, &pk_seed, 0, params.hp as u32, &mut root_adrs, &params);

        // Sign with leaf 0
        let mut sign_adrs = Address::new();
        sign_adrs.set_layer(0);
        sign_adrs.set_tree(0);
        let sig = xmss_sign(&msg, &sk_seed, 0, &pk_seed, &mut sign_adrs, &params);

        // Verify with wrong leaf index
        let mut verify_adrs = Address::new();
        verify_adrs.set_layer(0);
        verify_adrs.set_tree(0);
        let bad_root = xmss_pk_from_sig(1, &sig, &msg, &pk_seed, &mut verify_adrs, &params);

        assert_ne!(root, bad_root);
    }
}
