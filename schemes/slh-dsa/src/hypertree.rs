//! Hypertree implementation for SLH-DSA.
//!
//! The hypertree provides the multi-layer structure that connects multiple
//! XMSS trees together. This allows for a much larger number of signatures
//! than a single XMSS tree while keeping signatures relatively compact.
//!
//! # Structure
//!
//! The hypertree consists of d layers of XMSS trees:
//! - Layer 0 (bottom): Signs the FORS public key
//! - Layers 1 to d-1: Each XMSS tree signs the root of a tree in the layer below
//! - The root of the layer d-1 tree is the SLH-DSA public key
//!
//! # Addressing
//!
//! Each XMSS tree is identified by:
//! - Its layer (0 to d-1)
//! - Its tree index within the layer

use crate::address::Address;
use crate::params::Params;
use crate::xmss::{xmss_node, xmss_pk_from_sig, xmss_sign, XmssSig};

/// Hypertree signature structure.
#[derive(Clone, Debug)]
pub struct HypertreeSig {
    /// XMSS signatures for each layer (d signatures total)
    pub layers: Vec<XmssSig>,
}

impl HypertreeSig {
    /// Creates a new empty hypertree signature structure.
    pub fn new(params: &Params) -> Self {
        HypertreeSig {
            layers: Vec::with_capacity(params.d),
        }
    }
}

/// Generates a hypertree signature.
///
/// ht_sign(M, SK.seed, PK.seed, idx_tree, idx_leaf) - FIPS 205 Algorithm 10
///
/// Signs message M (typically a FORS public key) using the hypertree.
///
/// # Arguments
/// * `msg` - Message to sign (n bytes, typically FORS public key)
/// * `sk_seed` - Secret seed
/// * `pk_seed` - Public seed
/// * `idx_tree` - Tree index within layer 0
/// * `idx_leaf` - Leaf index within the layer 0 tree
/// * `params` - Parameter set
pub fn ht_sign(
    msg: &[u8],
    sk_seed: &[u8],
    pk_seed: &[u8],
    idx_tree: u64,
    idx_leaf: u32,
    params: &Params,
) -> HypertreeSig {
    let mut sig = HypertreeSig::new(params);

    // Set up address for layer 0
    let mut adrs = Address::new();
    adrs.set_layer(0);
    adrs.set_tree(idx_tree);

    // Sign at layer 0
    let sig_tmp = xmss_sign(msg, sk_seed, idx_leaf, pk_seed, &mut adrs, params);
    sig.layers.push(sig_tmp);

    // Compute root for next layer
    let mut root = xmss_node(sk_seed, pk_seed, 0, params.hp as u32, &mut adrs, params);

    // Sign at layers 1 through d-1
    for j in 1..params.d {
        // Compute indices for this layer
        let idx_leaf_j = extract_leaf_index(idx_tree, j, params);
        let idx_tree_j = extract_tree_index(idx_tree, j, params);

        adrs.set_layer(j as u32);
        adrs.set_tree(idx_tree_j);

        // Sign the root from the previous layer
        let sig_j = xmss_sign(&root, sk_seed, idx_leaf_j, pk_seed, &mut adrs, params);
        sig.layers.push(sig_j);

        // Compute root for next layer (if not at top)
        if j < params.d - 1 {
            root = xmss_node(sk_seed, pk_seed, 0, params.hp as u32, &mut adrs, params);
        }
    }

    sig
}

/// Verifies a hypertree signature.
///
/// ht_verify(M, SIG_HT, PK.seed, idx_tree, idx_leaf, PK.root) - FIPS 205 Algorithm 11
///
/// Verifies the signature and returns true if valid.
///
/// # Arguments
/// * `msg` - Signed message (n bytes)
/// * `sig` - Hypertree signature
/// * `pk_seed` - Public seed
/// * `idx_tree` - Tree index within layer 0
/// * `idx_leaf` - Leaf index within the layer 0 tree
/// * `pk_root` - Expected root (public key)
/// * `params` - Parameter set
pub fn ht_verify(
    msg: &[u8],
    sig: &HypertreeSig,
    pk_seed: &[u8],
    idx_tree: u64,
    idx_leaf: u32,
    pk_root: &[u8],
    params: &Params,
) -> bool {
    // Check signature has correct number of layers
    if sig.layers.len() != params.d {
        return false;
    }

    // Set up address for layer 0
    let mut adrs = Address::new();
    adrs.set_layer(0);
    adrs.set_tree(idx_tree);

    // Recover root at layer 0
    let mut node = xmss_pk_from_sig(idx_leaf, &sig.layers[0], msg, pk_seed, &mut adrs, params);

    // Verify through all layers
    for j in 1..params.d {
        // Compute indices for this layer
        let idx_leaf_j = extract_leaf_index(idx_tree, j, params);
        let idx_tree_j = extract_tree_index(idx_tree, j, params);

        adrs.set_layer(j as u32);
        adrs.set_tree(idx_tree_j);

        // Recover root at layer j
        node = xmss_pk_from_sig(idx_leaf_j, &sig.layers[j], &node, pk_seed, &mut adrs, params);
    }

    // Check if recovered root matches expected public key root
    node == pk_root
}

/// Extracts the leaf index for a given layer from the combined tree index.
///
/// At layer j, the leaf index is determined by bits
/// [(j-1)*h' .. j*h' - 1] of the original tree index.
fn extract_leaf_index(idx_tree: u64, layer: usize, params: &Params) -> u32 {
    let shift = (layer - 1) * params.hp;
    if shift >= 64 {
        return 0;
    }
    let mask = (1u64 << params.hp) - 1;
    ((idx_tree >> shift) & mask) as u32
}

/// Extracts the tree index for a given layer from the combined tree index.
///
/// At layer j, the tree index is idx_tree >> (j * h').
fn extract_tree_index(idx_tree: u64, layer: usize, params: &Params) -> u64 {
    let shift = layer * params.hp;
    if shift >= 64 {
        return 0;
    }
    idx_tree >> shift
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::SLH_DSA_SHAKE_128F;

    #[test]
    fn test_extract_indices() {
        let params = SLH_DSA_SHAKE_128F;

        // For 128f: hp = 3, so each layer uses 3 bits
        let idx_tree: u64 = 0b111_010_001; // Layer 2: 111, Layer 1: 010, Layer 0 leaf: 001

        // At layer 1, leaf index should be bits [0..2] of idx_tree
        let leaf_1 = extract_leaf_index(idx_tree, 1, &params);
        assert_eq!(leaf_1, 0b001);

        // At layer 2, leaf index should be bits [3..5] of idx_tree
        let leaf_2 = extract_leaf_index(idx_tree, 2, &params);
        assert_eq!(leaf_2, 0b010);
    }

    #[test]
    fn test_ht_sign_structure() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        let sig = ht_sign(&msg, &sk_seed, &pk_seed, 0, 0, &params);

        // Should have d layers
        assert_eq!(sig.layers.len(), params.d);

        // Each layer should be a valid XMSS signature
        for layer_sig in &sig.layers {
            assert_eq!(layer_sig.wots_sig.len(), params.wots_len);
            assert_eq!(layer_sig.auth.len(), params.hp);
        }
    }

    #[test]
    fn test_ht_sign_verify_roundtrip() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        // Compute the public key root (top of hypertree)
        let mut root_adrs = Address::new();
        root_adrs.set_layer((params.d - 1) as u32);
        root_adrs.set_tree(0);
        let pk_root = xmss_node(&sk_seed, &pk_seed, 0, params.hp as u32, &mut root_adrs, &params);

        // Sign
        let sig = ht_sign(&msg, &sk_seed, &pk_seed, 0, 0, &params);

        // Verify
        let valid = ht_verify(&msg, &sig, &pk_seed, 0, 0, &pk_root, &params);
        assert!(valid);
    }

    #[test]
    fn test_ht_verify_wrong_message() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];
        let wrong_msg = vec![4u8; params.n];

        // Compute public key root
        let mut root_adrs = Address::new();
        root_adrs.set_layer((params.d - 1) as u32);
        root_adrs.set_tree(0);
        let pk_root = xmss_node(&sk_seed, &pk_seed, 0, params.hp as u32, &mut root_adrs, &params);

        // Sign correct message
        let sig = ht_sign(&msg, &sk_seed, &pk_seed, 0, 0, &params);

        // Verify with wrong message
        let valid = ht_verify(&wrong_msg, &sig, &pk_seed, 0, 0, &pk_root, &params);
        assert!(!valid);
    }

    #[test]
    fn test_ht_verify_wrong_root() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];
        let wrong_root = vec![0u8; params.n];

        // Sign
        let sig = ht_sign(&msg, &sk_seed, &pk_seed, 0, 0, &params);

        // Verify with wrong root
        let valid = ht_verify(&msg, &sig, &pk_seed, 0, 0, &wrong_root, &params);
        assert!(!valid);
    }

    #[test]
    fn test_ht_sign_different_indices() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        // Sign with different leaf indices
        let sig0 = ht_sign(&msg, &sk_seed, &pk_seed, 0, 0, &params);
        let sig1 = ht_sign(&msg, &sk_seed, &pk_seed, 0, 1, &params);

        // Signatures should be different (different authentication paths)
        assert_ne!(sig0.layers[0].wots_sig, sig1.layers[0].wots_sig);
    }

    #[test]
    fn test_ht_sign_different_trees() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        // Sign with different tree indices
        let sig0 = ht_sign(&msg, &sk_seed, &pk_seed, 0, 0, &params);
        let sig1 = ht_sign(&msg, &sk_seed, &pk_seed, 1, 0, &params);

        // Signatures should be different
        // The layer 0 signatures will differ because they're in different trees
        assert_ne!(sig0.layers[0].auth, sig1.layers[0].auth);
    }
}
