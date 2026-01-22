//! WOTS+ (Winternitz One-Time Signature) implementation for SLH-DSA.
//!
//! WOTS+ is the fundamental one-time signature building block used in SLH-DSA.
//! It provides a mechanism to sign a fixed-size message using hash chains.
//!
//! # Overview
//!
//! WOTS+ works by:
//! 1. Converting the message to base-w representation
//! 2. Computing a checksum to prevent forgery
//! 3. Using the combined message||checksum to determine chain lengths
//! 4. Signing by revealing appropriate positions in hash chains
//!
//! # Security
//!
//! WOTS+ is a one-time signature scheme - each key pair should only be used once.
//! Multiple uses can allow forgery attacks.

use crate::address::{Address, AddressType};
use crate::hash::{f, prf, t_l};
use crate::params::{Params, LOG_W, W};

/// Computes a single WOTS+ hash chain.
///
/// chain(X, i, s, PK.seed, ADRS) - FIPS 205 Algorithm 2
///
/// Starting from value X, applies the F function s times starting at position i.
///
/// # Arguments
/// * `x` - Starting value (n bytes)
/// * `start` - Starting position in the chain
/// * `steps` - Number of steps to compute
/// * `pk_seed` - Public seed
/// * `adrs` - Address structure (will be modified)
/// * `params` - Parameter set
pub fn chain(
    x: &[u8],
    start: u32,
    steps: u32,
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    if steps == 0 {
        return x.to_vec();
    }

    let mut tmp = x.to_vec();
    for j in start..(start + steps) {
        adrs.set_hash(j);
        tmp = f(pk_seed, adrs, &tmp, params);
    }
    tmp
}

/// Generates a WOTS+ public key.
///
/// wots_PKgen(SK.seed, PK.seed, ADRS) - FIPS 205 Algorithm 4
///
/// Generates the public key by computing the full chain for each secret key element
/// and compressing the result.
pub fn wots_pk_gen(
    sk_seed: &[u8],
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    let mut sk_adrs = adrs.copy();
    sk_adrs.set_type(AddressType::WotsPrf);
    sk_adrs.set_keypair(adrs.keypair());

    let mut tmp = Vec::with_capacity(params.wots_len);

    for i in 0..params.wots_len {
        sk_adrs.set_chain(i as u32);
        let sk_i = prf(pk_seed, sk_seed, &sk_adrs, params);

        adrs.set_chain(i as u32);
        let pk_i = chain(&sk_i, 0, (W - 1) as u32, pk_seed, adrs, params);
        tmp.push(pk_i);
    }

    // Compress the public key using T_l
    let mut wots_pk_adrs = adrs.copy();
    wots_pk_adrs.set_type(AddressType::WotsPk);
    wots_pk_adrs.set_keypair(adrs.keypair());

    t_l(pk_seed, &wots_pk_adrs, &tmp, params)
}

/// Generates a WOTS+ signature.
///
/// wots_sign(M, SK.seed, PK.seed, ADRS) - FIPS 205 Algorithm 5
///
/// Signs an n-byte message by revealing appropriate chain positions.
pub fn wots_sign(
    msg: &[u8],
    sk_seed: &[u8],
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> Vec<Vec<u8>> {
    // Convert message to base-w representation
    let msg_base_w = base_w(msg, params.wots_len1());

    // Compute checksum
    let csum = compute_checksum(&msg_base_w, params);

    // Convert checksum to base-w (len2 elements)
    let csum_base_w = checksum_to_base_w(csum, params.wots_len2());

    // Combine message and checksum
    let mut len_x = msg_base_w;
    len_x.extend(csum_base_w);

    // Generate signature
    let mut sig = Vec::with_capacity(params.wots_len);

    let mut sk_adrs = adrs.copy();
    sk_adrs.set_type(AddressType::WotsPrf);
    sk_adrs.set_keypair(adrs.keypair());

    for i in 0..params.wots_len {
        sk_adrs.set_chain(i as u32);
        let sk_i = prf(pk_seed, sk_seed, &sk_adrs, params);

        adrs.set_chain(i as u32);
        let sig_i = chain(&sk_i, 0, len_x[i] as u32, pk_seed, adrs, params);
        sig.push(sig_i);
    }

    sig
}

/// Computes WOTS+ public key from signature.
///
/// wots_PKFromSig(sig, M, PK.seed, ADRS) - FIPS 205 Algorithm 6
///
/// Recovers the public key by completing the chains from the signature values.
pub fn wots_pk_from_sig(
    sig: &[Vec<u8>],
    msg: &[u8],
    pk_seed: &[u8],
    adrs: &mut Address,
    params: &Params,
) -> Vec<u8> {
    // Convert message to base-w representation
    let msg_base_w = base_w(msg, params.wots_len1());

    // Compute checksum
    let csum = compute_checksum(&msg_base_w, params);

    // Convert checksum to base-w
    let csum_base_w = checksum_to_base_w(csum, params.wots_len2());

    // Combine message and checksum
    let mut len_x = msg_base_w;
    len_x.extend(csum_base_w);

    // Compute public key elements from signature
    let mut tmp = Vec::with_capacity(params.wots_len);

    for i in 0..params.wots_len {
        adrs.set_chain(i as u32);
        let steps = (W - 1) as u32 - len_x[i] as u32;
        let pk_i = chain(&sig[i], len_x[i] as u32, steps, pk_seed, adrs, params);
        tmp.push(pk_i);
    }

    // Compress the public key
    let mut wots_pk_adrs = adrs.copy();
    wots_pk_adrs.set_type(AddressType::WotsPk);
    wots_pk_adrs.set_keypair(adrs.keypair());

    t_l(pk_seed, &wots_pk_adrs, &tmp, params)
}

/// Converts a byte array to base-w representation.
///
/// Each byte is split into LOG_W-bit chunks (nibbles for w=16).
fn base_w(input: &[u8], out_len: usize) -> Vec<u8> {
    let mut output = Vec::with_capacity(out_len);
    let mut bits_left = 0u32;
    let mut total = 0u32;
    let mut idx = 0;

    let w_mask = (W - 1) as u32;

    for _ in 0..out_len {
        // Load more bits if needed
        while bits_left < LOG_W as u32 {
            if idx < input.len() {
                total = (total << 8) | (input[idx] as u32);
                idx += 1;
            } else {
                total <<= 8;
            }
            bits_left += 8;
        }

        // Extract LOG_W bits
        bits_left -= LOG_W as u32;
        output.push(((total >> bits_left) & w_mask) as u8);
    }

    output
}

/// Computes the WOTS+ checksum.
///
/// csum = sum(w - 1 - msg[i]) for all i in 0..len1
fn compute_checksum(msg: &[u8], params: &Params) -> u32 {
    let mut csum = 0u32;
    for &m in msg.iter().take(params.wots_len1()) {
        csum += (W - 1) as u32 - m as u32;
    }
    csum
}

/// Converts checksum to base-w representation.
///
/// The checksum is converted to len2 base-w digits, most significant first.
fn checksum_to_base_w(mut csum: u32, len2: usize) -> Vec<u8> {
    // Shift checksum left to align with LOG_W * len2 bits
    csum <<= (8 - ((len2 * LOG_W) % 8)) % 8;

    let mut output = Vec::with_capacity(len2);
    let w_mask = (W - 1) as u32;

    for i in 0..len2 {
        // Extract digits from most significant to least significant
        let shift = (len2 - 1 - i) * LOG_W;
        output.push(((csum >> shift) & w_mask) as u8);
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::SLH_DSA_SHAKE_128S;

    #[test]
    fn test_base_w_conversion() {
        // For w=16, each byte becomes 2 nibbles
        let input = [0xAB, 0xCD];
        let output = base_w(&input, 4);
        assert_eq!(output, vec![0xA, 0xB, 0xC, 0xD]);
    }

    #[test]
    fn test_base_w_single_byte() {
        let input = [0x12];
        let output = base_w(&input, 2);
        assert_eq!(output, vec![0x1, 0x2]);
    }

    #[test]
    fn test_checksum_all_zeros() {
        let params = SLH_DSA_SHAKE_128S;
        let msg = vec![0u8; params.wots_len1()];
        let csum = compute_checksum(&msg, &params);
        // Each element contributes (w-1) = 15, so total = 15 * len1
        assert_eq!(csum, 15 * params.wots_len1() as u32);
    }

    #[test]
    fn test_checksum_all_max() {
        let params = SLH_DSA_SHAKE_128S;
        let msg = vec![(W - 1) as u8; params.wots_len1()];
        let csum = compute_checksum(&msg, &params);
        // Each element contributes 0, so total = 0
        assert_eq!(csum, 0);
    }

    #[test]
    fn test_chain_zero_steps() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let x = vec![42u8; params.n];
        let mut adrs = Address::new();
        adrs.set_type(AddressType::WotsHash);

        let result = chain(&x, 0, 0, &pk_seed, &mut adrs, &params);
        assert_eq!(result, x);
    }

    #[test]
    fn test_chain_one_step() {
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let x = vec![1u8; params.n];
        let mut adrs = Address::new();
        adrs.set_type(AddressType::WotsHash);

        let result = chain(&x, 0, 1, &pk_seed, &mut adrs, &params);
        assert_eq!(result.len(), params.n);
        assert_ne!(result, x);
    }

    #[test]
    fn test_chain_additivity() {
        // chain(x, i, s1+s2) should equal chain(chain(x, i, s1), i+s1, s2)
        let params = SLH_DSA_SHAKE_128S;
        let pk_seed = vec![0u8; params.n];
        let x = vec![1u8; params.n];

        let mut adrs1 = Address::new();
        adrs1.set_type(AddressType::WotsHash);

        let mut adrs2 = Address::new();
        adrs2.set_type(AddressType::WotsHash);

        let mut adrs3 = Address::new();
        adrs3.set_type(AddressType::WotsHash);

        // Compute chain(x, 0, 5) in one go
        let result1 = chain(&x, 0, 5, &pk_seed, &mut adrs1, &params);

        // Compute chain(x, 0, 3) then chain(result, 3, 2)
        let intermediate = chain(&x, 0, 3, &pk_seed, &mut adrs2, &params);
        let result2 = chain(&intermediate, 3, 2, &pk_seed, &mut adrs3, &params);

        assert_eq!(result1, result2);
    }

    #[test]
    fn test_wots_pk_gen_deterministic() {
        let params = SLH_DSA_SHAKE_128S;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];

        let mut adrs1 = Address::new();
        adrs1.set_type(AddressType::WotsHash);
        adrs1.set_keypair(0);

        let mut adrs2 = Address::new();
        adrs2.set_type(AddressType::WotsHash);
        adrs2.set_keypair(0);

        let pk1 = wots_pk_gen(&sk_seed, &pk_seed, &mut adrs1, &params);
        let pk2 = wots_pk_gen(&sk_seed, &pk_seed, &mut adrs2, &params);

        assert_eq!(pk1, pk2);
        assert_eq!(pk1.len(), params.n);
    }

    #[test]
    fn test_wots_sign_verify_roundtrip() {
        let params = SLH_DSA_SHAKE_128S;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        // Generate public key
        let mut adrs = Address::new();
        adrs.set_type(AddressType::WotsHash);
        adrs.set_keypair(0);
        let pk = wots_pk_gen(&sk_seed, &pk_seed, &mut adrs, &params);

        // Sign
        let mut sign_adrs = Address::new();
        sign_adrs.set_type(AddressType::WotsHash);
        sign_adrs.set_keypair(0);
        let sig = wots_sign(&msg, &sk_seed, &pk_seed, &mut sign_adrs, &params);

        // Verify
        let mut verify_adrs = Address::new();
        verify_adrs.set_type(AddressType::WotsHash);
        verify_adrs.set_keypair(0);
        let pk_recovered = wots_pk_from_sig(&sig, &msg, &pk_seed, &mut verify_adrs, &params);

        assert_eq!(pk, pk_recovered);
    }

    #[test]
    fn test_wots_signature_size() {
        let params = SLH_DSA_SHAKE_128S;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];

        let mut adrs = Address::new();
        adrs.set_type(AddressType::WotsHash);
        adrs.set_keypair(0);

        let sig = wots_sign(&msg, &sk_seed, &pk_seed, &mut adrs, &params);

        assert_eq!(sig.len(), params.wots_len);
        for s in &sig {
            assert_eq!(s.len(), params.n);
        }
    }

    #[test]
    fn test_wots_different_keypairs() {
        let params = SLH_DSA_SHAKE_128S;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];

        let mut adrs1 = Address::new();
        adrs1.set_type(AddressType::WotsHash);
        adrs1.set_keypair(0);

        let mut adrs2 = Address::new();
        adrs2.set_type(AddressType::WotsHash);
        adrs2.set_keypair(1);

        let pk1 = wots_pk_gen(&sk_seed, &pk_seed, &mut adrs1, &params);
        let pk2 = wots_pk_gen(&sk_seed, &pk_seed, &mut adrs2, &params);

        // Different keypair indices should give different public keys
        assert_ne!(pk1, pk2);
    }

    #[test]
    fn test_wots_wrong_message() {
        let params = SLH_DSA_SHAKE_128S;
        let sk_seed = vec![1u8; params.n];
        let pk_seed = vec![2u8; params.n];
        let msg = vec![3u8; params.n];
        let wrong_msg = vec![4u8; params.n];

        // Generate public key
        let mut adrs = Address::new();
        adrs.set_type(AddressType::WotsHash);
        adrs.set_keypair(0);
        let pk = wots_pk_gen(&sk_seed, &pk_seed, &mut adrs, &params);

        // Sign original message
        let mut sign_adrs = Address::new();
        sign_adrs.set_type(AddressType::WotsHash);
        sign_adrs.set_keypair(0);
        let sig = wots_sign(&msg, &sk_seed, &pk_seed, &mut sign_adrs, &params);

        // Try to verify with wrong message
        let mut verify_adrs = Address::new();
        verify_adrs.set_type(AddressType::WotsHash);
        verify_adrs.set_keypair(0);
        let pk_recovered = wots_pk_from_sig(&sig, &wrong_msg, &pk_seed, &mut verify_adrs, &params);

        // Should not match
        assert_ne!(pk, pk_recovered);
    }
}
