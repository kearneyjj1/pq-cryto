//! FORS (Forest of Random Subsets) implementation for SLH-DSA.
//!
//! FORS is a few-time signature scheme used to sign the message digest in SLH-DSA.
//! It provides better efficiency than using WOTS+ directly for the message.
//!
//! # Structure
//!
//! FORS consists of k binary trees, each of height a:
//! - Each tree has 2^a leaves
//! - The message digest determines one leaf index per tree
//! - The signature reveals the selected leaves and their authentication paths
//!
//! # Signature
//!
//! A FORS signature contains:
//! - k secret key values (the selected leaves)
//! - k authentication paths (a nodes each)

use crate::address::{Address, AddressType};
use crate::hash::{f, h, prf, t_l};
use crate::params::Params;

/// FORS signature structure.
#[derive(Clone, Debug)]
pub struct ForsSig {
    /// Secret key values for each tree (k values, each n bytes)
    pub sk: Vec<Vec<u8>>,
    /// Authentication paths for each tree (k paths, each with a nodes)
    pub auth: Vec<Vec<Vec<u8>>>,
}

impl ForsSig {
    /// Creates a new empty FORS signature structure.
    pub fn new(params: &Params) -> Self {
        ForsSig {
            sk: Vec::with_capacity(params.k),
            auth: Vec::with_capacity(params.k),
        }
    }
}

/// Generates a FORS secret key element.
///
/// fors_SKgen(SK.seed, PK.seed, ADRS, idx) - FIPS 205 Algorithm 12
fn fors_sk_gen(
    sk_seed: &[u8],
    pk_seed: &[u8],
    adrs: &mut Address,
    idx: u32,
    params: &Params,
) -> Vec<u8> {
    // Save keypair before type change (set_type clears type-specific fields)
    let keypair = adrs.keypair();

    let mut sk_adrs = adrs.copy();
    sk_adrs.set_type(AddressType::ForsPrf);
    sk_adrs.set_keypair(keypair); // Restore keypair after type change
    sk_adrs.set_tree_height(0);
    sk_adrs.set_tree_index(idx);

    prf(pk_seed, sk_seed, &sk_adrs, params)
}

/// Computes a FORS tree node.
///
/// fors_node(SK.seed, PK.seed, i, z, ADRS) - FIPS 205 Algorithm 13
///
/// Computes the node at height z and position i in the specified FORS tree.
#[allow(dead_code)]
fn fors_node(
    sk_seed: &[u8],
    pk_seed: &[u8],
    i: u32,
    z: u32,
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    if z == 0 {
        // Leaf node: hash the secret key value
        let sk = fors_sk_gen(sk_seed, pk_seed, adrs, i, params);

        adrs.set_type(AddressType::ForsTree);
        adrs.set_tree_height(0);
        adrs.set_tree_index(i);

        return f(pk_seed, adrs, &sk, params);
    }

    // Internal node: hash the two children
    let left = fors_node(sk_seed, pk_seed, 2 * i, z - 1, adrs, params);
    let right = fors_node(sk_seed, pk_seed, 2 * i + 1, z - 1, adrs, params);

    adrs.set_type(AddressType::ForsTree);
    adrs.set_tree_height(z);
    adrs.set_tree_index(i);

    h(pk_seed, adrs, &left, &right, params)
}

/// Generates a FORS signature.
///
/// fors_sign(md, SK.seed, PK.seed, ADRS) - FIPS 205 Algorithm 14
///
/// Signs the message digest md using FORS.
pub fn fors_sign(
    md: &[u8],
    sk_seed: &[u8],
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> ForsSig {
    let mut sig = ForsSig::new(params);
    let indices = message_to_indices(md, params);

    for i in 0..params.k {
        let idx = indices[i];

        // Compute the global leaf index for this tree
        let tree_offset = (i as u32) * (1 << params.a);
        let global_idx = tree_offset + idx;

        // Get the secret key value
        let sk_val = fors_sk_gen(sk_seed, pk_seed, adrs, global_idx, params);
        sig.sk.push(sk_val);

        // Compute authentication path
        let mut auth_path = Vec::with_capacity(params.a);
        for j in 0..params.a {
            // At height j, find the sibling node
            let sibling_idx = (idx >> j) ^ 1;

            // Compute the sibling node at height j
            // The node at position sibling_idx at height j covers leaves
            // from sibling_idx * 2^j to (sibling_idx + 1) * 2^j - 1
            let node = fors_tree_node(sk_seed, pk_seed, i, sibling_idx as u32, j as u32, adrs, params);
            auth_path.push(node);
        }
        sig.auth.push(auth_path);
    }

    sig
}

/// Computes a node in a specific FORS tree.
///
/// Helper function that computes a node at position `pos` and height `height`
/// within FORS tree `tree_idx`.
fn fors_tree_node(
    sk_seed: &[u8],
    pk_seed: &[u8],
    tree_idx: usize,
    pos: u32,
    height: u32,
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    // Save keypair before type changes
    let keypair = adrs.keypair();
    let tree_offset = (tree_idx as u32) * (1 << params.a);

    if height == 0 {
        // Leaf: global index is tree_offset + pos
        let global_idx = tree_offset + pos;
        let sk = fors_sk_gen(sk_seed, pk_seed, adrs, global_idx, params);

        adrs.set_type(AddressType::ForsTree);
        adrs.set_keypair(keypair);
        adrs.set_tree_height(0);
        adrs.set_tree_index(global_idx);

        return f(pk_seed, adrs, &sk, params);
    }

    // Internal node
    let left = fors_tree_node(sk_seed, pk_seed, tree_idx, 2 * pos, height - 1, adrs, params);
    let right = fors_tree_node(sk_seed, pk_seed, tree_idx, 2 * pos + 1, height - 1, adrs, params);

    adrs.set_type(AddressType::ForsTree);
    adrs.set_keypair(keypair);
    adrs.set_tree_height(height);
    adrs.set_tree_index(tree_offset + (pos << height));

    h(pk_seed, adrs, &left, &right, params)
}

/// Computes FORS public key from signature.
///
/// fors_pkFromSig(SIG_FORS, md, PK.seed, ADRS) - FIPS 205 Algorithm 15
///
/// Reconstructs the FORS public key from a signature and message digest.
pub fn fors_pk_from_sig(
    sig: &ForsSig,
    md: &[u8],
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    // Save keypair before type changes
    let keypair = adrs.keypair();
    let indices = message_to_indices(md, params);
    let mut roots = Vec::with_capacity(params.k);

    for i in 0..params.k {
        let idx = indices[i];
        let tree_offset = (i as u32) * (1 << params.a);
        let global_idx = tree_offset + idx;

        // Hash the secret key to get the leaf value
        adrs.set_type(AddressType::ForsTree);
        adrs.set_keypair(keypair);
        adrs.set_tree_height(0);
        adrs.set_tree_index(global_idx);
        let mut node = f(pk_seed, adrs, &sig.sk[i], params);

        // Climb the authentication path
        for j in 0..params.a {
            adrs.set_tree_height((j + 1) as u32);

            // Compute the parent node index
            let parent_idx = tree_offset + ((idx >> (j + 1)) << (j + 1));

            if (idx >> j) & 1 == 0 {
                // Current node is left child
                adrs.set_tree_index(parent_idx);
                node = h(pk_seed, adrs, &node, &sig.auth[i][j], params);
            } else {
                // Current node is right child
                adrs.set_tree_index(parent_idx);
                node = h(pk_seed, adrs, &sig.auth[i][j], &node, params);
            }
        }

        roots.push(node);
    }

    // Compress the roots into the FORS public key
    adrs.set_type(AddressType::ForsPk);
    adrs.set_keypair(keypair);

    t_l(pk_seed, adrs, &roots, params)
}

/// Generates the FORS public key.
///
/// Computes the FORS public key by computing all tree roots and compressing them.
pub fn fors_pk_gen(
    sk_seed: &[u8],
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    // Save keypair before type changes
    let keypair = adrs.keypair();
    let mut roots = Vec::with_capacity(params.k);

    for i in 0..params.k {
        // Compute root of tree i
        let root = fors_tree_node(sk_seed, pk_seed, i, 0, params.a as u32, adrs, params);
        roots.push(root);
    }

    // Compress the roots
    adrs.set_type(AddressType::ForsPk);
    adrs.set_keypair(keypair);
    t_l(pk_seed, adrs, &roots, params)
}

/// Converts message digest to FORS indices.
///
/// Extracts k indices of a bits each from the message digest.
fn message_to_indices(md: &[u8], params: &Params) -> Vec<u32> {
    let mut indices = Vec::with_capacity(params.k);
    let mut bits_left = 0u32;
    let mut total = 0u64;
    let mut byte_idx = 0;

    let mask = (1u32 << params.a) - 1;

    for _ in 0..params.k {
        // Load more bits if needed
        while bits_left < params.a as u32 {
            if byte_idx < md.len() {
                total = (total << 8) | (md[byte_idx] as u64);
                byte_idx += 1;
            } else {
                total <<= 8;
            }
            bits_left += 8;
        }

        // Extract a bits
        bits_left -= params.a as u32;
        indices.push(((total >> bits_left) as u32) & mask);
    }

    indices
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::SLH_DSA_SHAKE_128F;

    #[test]
    fn test_message_to_indices() {
        let params = SLH_DSA_SHAKE_128F;
        let md = vec![0xFF; params.md_len()];

        let indices = message_to_indices(&md, &params);

        assert_eq!(indices.len(), params.k);
        for &idx in &indices {
            assert!(idx < (1 << params.a));
        }
    }

    #[test]
    fn test_message_to_indices_deterministic() {
        let params = SLH_DSA_SHAKE_128F;
        let md = vec![0x12, 0x34, 0x56, 0x78, 0x9A, 0xBC, 0xDE, 0xF0];

        let indices1 = message_to_indices(&md, &params);
        let indices2 = message_to_indices(&md, &params);

        assert_eq!(indices1, indices2);
    }

    #[test]
    fn test_fors_sk_gen_deterministic() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];

        let mut adrs1 = Address::new();
        adrs1.set_keypair(0);

        let mut adrs2 = Address::new();
        adrs2.set_keypair(0);

        let sk1 = fors_sk_gen(&sk_seed, &pk_seed, &mut adrs1, 0, &params);
        let sk2 = fors_sk_gen(&sk_seed, &pk_seed, &mut adrs2, 0, &params);

        assert_eq!(sk1, sk2);
        assert_eq!(sk1.len(), params.n);
    }

    #[test]
    fn test_fors_sk_gen_different_indices() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let mut adrs = Address::new();
        adrs.set_keypair(0);

        let sk0 = fors_sk_gen(&sk_seed, &pk_seed, &mut adrs, 0, &params);
        let sk1 = fors_sk_gen(&sk_seed, &pk_seed, &mut adrs, 1, &params);

        assert_ne!(sk0, sk1);
    }

    #[test]
    fn test_fors_sign_structure() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let md = vec![3u8; params.md_len()];

        let mut adrs = Address::new();
        adrs.set_keypair(0);

        let sig = fors_sign(&md, &sk_seed, &pk_seed, &mut adrs, &params);

        // Check structure sizes
        assert_eq!(sig.sk.len(), params.k);
        assert_eq!(sig.auth.len(), params.k);

        for sk in &sig.sk {
            assert_eq!(sk.len(), params.n);
        }

        for auth_path in &sig.auth {
            assert_eq!(auth_path.len(), params.a);
            for node in auth_path {
                assert_eq!(node.len(), params.n);
            }
        }
    }

    #[test]
    fn test_fors_sign_verify_roundtrip() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let md = vec![3u8; params.md_len()];

        // Generate public key
        let mut pk_adrs = Address::new();
        pk_adrs.set_keypair(0);
        let pk = fors_pk_gen(&sk_seed, &pk_seed, &mut pk_adrs, &params);

        // Sign
        let mut sign_adrs = Address::new();
        sign_adrs.set_keypair(0);
        let sig = fors_sign(&md, &sk_seed, &pk_seed, &mut sign_adrs, &params);

        // Verify
        let mut verify_adrs = Address::new();
        verify_adrs.set_keypair(0);
        let pk_recovered = fors_pk_from_sig(&sig, &md, &pk_seed, &mut verify_adrs, &params);

        assert_eq!(pk, pk_recovered);
    }

    #[test]
    fn test_fors_wrong_message() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let md = vec![3u8; params.md_len()];
        let wrong_md = vec![4u8; params.md_len()];

        // Generate public key
        let mut pk_adrs = Address::new();
        pk_adrs.set_keypair(0);
        let pk = fors_pk_gen(&sk_seed, &pk_seed, &mut pk_adrs, &params);

        // Sign correct message
        let mut sign_adrs = Address::new();
        sign_adrs.set_keypair(0);
        let sig = fors_sign(&md, &sk_seed, &pk_seed, &mut sign_adrs, &params);

        // Verify with wrong message
        let mut verify_adrs = Address::new();
        verify_adrs.set_keypair(0);
        let bad_pk = fors_pk_from_sig(&sig, &wrong_md, &pk_seed, &mut verify_adrs, &params);

        assert_ne!(pk, bad_pk);
    }

    #[test]
    fn test_fors_pk_gen_deterministic() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];

        let mut adrs1 = Address::new();
        adrs1.set_keypair(0);

        let mut adrs2 = Address::new();
        adrs2.set_keypair(0);

        let pk1 = fors_pk_gen(&sk_seed, &pk_seed, &mut adrs1, &params);
        let pk2 = fors_pk_gen(&sk_seed, &pk_seed, &mut adrs2, &params);

        assert_eq!(pk1, pk2);
        assert_eq!(pk1.len(), params.n);
    }

    #[test]
    fn test_fors_different_keypairs() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];

        let mut adrs1 = Address::new();
        adrs1.set_keypair(0);

        let mut adrs2 = Address::new();
        adrs2.set_keypair(1);

        let pk1 = fors_pk_gen(&sk_seed, &pk_seed, &mut adrs1, &params);
        let pk2 = fors_pk_gen(&sk_seed, &pk_seed, &mut adrs2, &params);

        // Different keypair addresses should give different public keys
        assert_ne!(pk1, pk2);
    }
}
