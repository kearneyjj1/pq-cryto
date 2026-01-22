//! Hash functions for SLH-DSA using SHAKE256.
//!
//! This module implements all the hash functions required by SLH-DSA (FIPS 205)
//! using SHAKE256 as the underlying XOF (Extendable Output Function).
//!
//! # Functions
//!
//! - `PRF`: Pseudorandom function for deriving secret values
//! - `PRF_msg`: Message-dependent PRF for signature randomization
//! - `H_msg`: Hash function producing message digest and tree indices
//! - `F`: Tweakable hash with single n-byte input (for chain hashing)
//! - `H`: Tweakable hash with 2n-byte input (for tree nodes)
//! - `T_l`: Tweakable hash with variable-length input (for WOTS+ PK compression)

use crate::address::Address;
use crate::params::Params;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

/// PRF function for deriving secret values.
///
/// PRF(PK.seed, SK.seed, ADRS) = SHAKE256(PK.seed || ADRS || SK.seed, 8n)
///
/// Used for deriving WOTS+ and FORS secret key values.
pub fn prf(pk_seed: &[u8], sk_seed: &[u8], adrs: &Address, params: &Params) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(pk_seed);
    hasher.update(adrs.as_bytes());
    hasher.update(sk_seed);

    let mut reader = hasher.finalize_xof();
    let mut output = vec![0u8; params.n];
    reader.read(&mut output);
    output
}

/// PRF_msg function for signature randomization.
///
/// PRF_msg(SK.prf, OptRand, M) = SHAKE256(SK.prf || OptRand || M, 8n)
///
/// Produces the randomizer R for the signature.
pub fn prf_msg(sk_prf: &[u8], opt_rand: &[u8], msg: &[u8], params: &Params) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(sk_prf);
    hasher.update(opt_rand);
    hasher.update(msg);

    let mut reader = hasher.finalize_xof();
    let mut output = vec![0u8; params.n];
    reader.read(&mut output);
    output
}

/// H_msg function for hashing the message.
///
/// H_msg(R, PK.seed, PK.root, M) produces:
/// - Message digest (md) of k*a bits for FORS
/// - Tree index (idx_tree) of h - h' bits
/// - Leaf index (idx_leaf) of h' bits
///
/// Returns (md, idx_tree, idx_leaf)
pub fn h_msg(
    r: &[u8],
    pk_seed: &[u8],
    pk_root: &[u8],
    msg: &[u8],
    params: &Params,
) -> (Vec<u8>, u64, u32) {
    let mut hasher = Shake256::default();
    hasher.update(r);
    hasher.update(pk_seed);
    hasher.update(pk_root);
    hasher.update(msg);

    let mut reader = hasher.finalize_xof();

    // Calculate required output bytes
    let md_len = params.md_len();
    let tree_bits = params.tree_bits();
    let leaf_bits = params.leaf_bits();

    // Total bits needed for indices: tree_bits + leaf_bits = h bits
    let idx_bytes = (tree_bits + leaf_bits + 7) / 8;

    // Read message digest
    let mut md = vec![0u8; md_len];
    reader.read(&mut md);

    // Read index bytes
    let mut idx_raw = vec![0u8; idx_bytes];
    reader.read(&mut idx_raw);

    // Extract tree index and leaf index from the raw bytes
    // Total h bits: upper (h - h') bits are tree index, lower h' bits are leaf index
    let total_bits = tree_bits + leaf_bits;

    // Convert idx_raw to a large integer and extract the indices
    let mut idx_val: u128 = 0;
    for &b in &idx_raw {
        idx_val = (idx_val << 8) | (b as u128);
    }

    // Shift to align to the number of bits we need
    let extra_bits = idx_bytes * 8 - total_bits;
    idx_val >>= extra_bits;

    // Extract leaf index (lower h' bits)
    let leaf_mask = (1u128 << leaf_bits) - 1;
    let idx_leaf = (idx_val & leaf_mask) as u32;

    // Extract tree index (upper h - h' bits)
    let idx_tree = (idx_val >> leaf_bits) as u64;

    // Mask tree index to actual bit width
    let tree_mask = if tree_bits >= 64 {
        u64::MAX
    } else {
        (1u64 << tree_bits) - 1
    };
    let idx_tree = idx_tree & tree_mask;

    (md, idx_tree, idx_leaf)
}

/// F function: tweakable hash with single n-byte input.
///
/// F(PK.seed, ADRS, M1) = SHAKE256(PK.seed || ADRS || M1, 8n)
///
/// Used for hash chain computations in WOTS+ and FORS leaf generation.
pub fn f(pk_seed: &[u8], adrs: &Address, m1: &[u8], params: &Params) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(pk_seed);
    hasher.update(adrs.as_bytes());
    hasher.update(m1);

    let mut reader = hasher.finalize_xof();
    let mut output = vec![0u8; params.n];
    reader.read(&mut output);
    output
}

/// H function: tweakable hash with 2n-byte input.
///
/// H(PK.seed, ADRS, M1 || M2) = SHAKE256(PK.seed || ADRS || M1 || M2, 8n)
///
/// Used for computing internal Merkle tree nodes.
pub fn h(pk_seed: &[u8], adrs: &Address, m1: &[u8], m2: &[u8], params: &Params) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(pk_seed);
    hasher.update(adrs.as_bytes());
    hasher.update(m1);
    hasher.update(m2);

    let mut reader = hasher.finalize_xof();
    let mut output = vec![0u8; params.n];
    reader.read(&mut output);
    output
}

/// T_l function: tweakable hash with variable-length input.
///
/// T_l(PK.seed, ADRS, M) = SHAKE256(PK.seed || ADRS || M, 8n)
///
/// Used for compressing WOTS+ public keys and FORS roots.
pub fn t_l(pk_seed: &[u8], adrs: &Address, messages: &[Vec<u8>], params: &Params) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(pk_seed);
    hasher.update(adrs.as_bytes());
    for m in messages {
        hasher.update(m);
    }

    let mut reader = hasher.finalize_xof();
    let mut output = vec![0u8; params.n];
    reader.read(&mut output);
    output
}

/// T_l function variant that takes a flat byte slice.
///
/// Useful when the input is already concatenated.
pub fn t_l_bytes(pk_seed: &[u8], adrs: &Address, message: &[u8], params: &Params) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(pk_seed);
    hasher.update(adrs.as_bytes());
    hasher.update(message);

    let mut reader = hasher.finalize_xof();
    let mut output = vec![0u8; params.n];
    reader.read(&mut output);
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::address::AddressType;
    use crate::params::SLH_DSA_SHAKE_128S;

    #[test]
    fn test_prf_output_size() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let sk_seed = vec![1u8; params.n];
        let adrs = Address::new();

        let output = prf(&pk_seed, &sk_seed, &adrs, &params);
        assert_eq!(output.len(), params.n);
    }

    #[test]
    fn test_prf_deterministic() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![42u8; params.n];
        let sk_seed = vec![43u8; params.n];
        let adrs = Address::new();

        let output1 = prf(&pk_seed, &sk_seed, &adrs, &params);
        let output2 = prf(&pk_seed, &sk_seed, &adrs, &params);
        assert_eq!(output1, output2);
    }

    #[test]
    fn test_prf_different_inputs() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let sk_seed1 = vec![1u8; params.n];
        let sk_seed2 = vec![2u8; params.n];
        let adrs = Address::new();

        let output1 = prf(&pk_seed, &sk_seed1, &adrs, &params);
        let output2 = prf(&pk_seed, &sk_seed2, &adrs, &params);
        assert_ne!(output1, output2);
    }

    #[test]
    fn test_prf_msg_output_size() {
        let params = SLH_DSA_SHAKE_128S;
        let sk_prf = vec![0u8; params.n];
        let opt_rand = vec![1u8; params.n];
        let msg = b"test message";

        let output = prf_msg(&sk_prf, &opt_rand, msg, &params);
        assert_eq!(output.len(), params.n);
    }

    #[test]
    fn test_h_msg_output_sizes() {
        let params = SLH_DSA_SHAKE_128S;
        let r = vec![0u8; params.n];
        let pk_seed = vec![1u8; params.n];
        let pk_root = vec![2u8; params.n];
        let msg = b"test message";

        let (md, idx_tree, idx_leaf) = h_msg(&r, &pk_seed, &pk_root, msg, &params);

        // Check md size
        assert_eq!(md.len(), params.md_len());

        // Check idx_leaf is within bounds (< 2^h')
        assert!(idx_leaf < (1 << params.hp));

        // Check idx_tree is within bounds (< 2^(h-h'))
        let tree_max = 1u64 << params.tree_bits();
        assert!(idx_tree < tree_max);
    }

    #[test]
    fn test_f_output_size() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let adrs = Address::new();
        let m1 = vec![1u8; params.n];

        let output = f(&pk_seed, &adrs, &m1, &params);
        assert_eq!(output.len(), params.n);
    }

    #[test]
    fn test_f_different_addresses() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let m1 = vec![1u8; params.n];

        let mut adrs1 = Address::new();
        adrs1.set_type(AddressType::WotsHash);
        adrs1.set_hash(0);

        let mut adrs2 = Address::new();
        adrs2.set_type(AddressType::WotsHash);
        adrs2.set_hash(1);

        let output1 = f(&pk_seed, &adrs1, &m1, &params);
        let output2 = f(&pk_seed, &adrs2, &m1, &params);
        assert_ne!(output1, output2);
    }

    #[test]
    fn test_h_output_size() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let adrs = Address::new();
        let m1 = vec![1u8; params.n];
        let m2 = vec![2u8; params.n];

        let output = h(&pk_seed, &adrs, &m1, &m2, &params);
        assert_eq!(output.len(), params.n);
    }

    #[test]
    fn test_t_l_output_size() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let adrs = Address::new();
        let messages: Vec<Vec<u8>> = (0..5).map(|i| vec![i; params.n]).collect();

        let output = t_l(&pk_seed, &adrs, &messages, &params);
        assert_eq!(output.len(), params.n);
    }

    #[test]
    fn test_t_l_bytes_equivalent() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let adrs = Address::new();

        let m1 = vec![1u8; params.n];
        let m2 = vec![2u8; params.n];
        let messages = vec![m1.clone(), m2.clone()];

        let mut concatenated = m1;
        concatenated.extend(&m2);

        let output1 = t_l(&pk_seed, &adrs, &messages, &params);
        let output2 = t_l_bytes(&pk_seed, &adrs, &concatenated, &params);
        assert_eq!(output1, output2);
    }
}
