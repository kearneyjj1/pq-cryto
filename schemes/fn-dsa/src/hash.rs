//! Hash functions for FALCON.
//!
//! This module implements the hash-to-point function using SHAKE256,
//! which maps a message and nonce to a polynomial in Z_q[X]/(X^n + 1).

use sha3::{Shake256, digest::{ExtendableOutput, Update, XofReader}};
use crate::field::Zq;
use crate::params::{Params, Q, NONCE_SIZE};

/// Hashes a message and nonce to a polynomial in Z_q (constant-time).
///
/// Uses SHAKE256 to generate a uniformly random polynomial c ∈ Z_q[X]/(X^n + 1).
/// The output polynomial has coefficients in [0, q-1].
///
/// Implements constant-time rejection sampling matching the Falcon reference
/// `hash_to_point_ct` to prevent timing side-channel leakage:
/// - Input: SHAKE256(nonce || message)
/// - Reads 2 bytes big-endian as a 16-bit value w
/// - Accepts if w < 5*q = 61445, rejects otherwise
/// - Reduces w mod q to get the coefficient
/// - All samples are processed identically regardless of acceptance
pub fn hash_to_point(message: &[u8], nonce: &[u8], params: &Params) -> Vec<Zq> {
    let n = params.n;
    let mut result = vec![Zq::ZERO; n];

    let mut hasher = Shake256::default();
    hasher.update(nonce);
    hasher.update(message);
    let mut reader = hasher.finalize_xof();

    // Oversample: read enough 2-byte pairs to guarantee >= n accepted values.
    // Acceptance rate is 61445/65536 ≈ 93.8%. For n=512, 840 samples gives
    // a negligible probability of insufficient accepted values.
    let num_samples = n + n / 4 + 200;
    let mut buf = vec![0u8; num_samples * 2];
    reader.read(&mut buf);

    // Pass 1: evaluate all samples, store value + validity
    let mut samples = Vec::with_capacity(num_samples);
    for k in 0..num_samples {
        let w = ((buf[2 * k] as u32) << 8) | (buf[2 * k + 1] as u32);
        let valid = w < 61445;
        let val = (w % (Q as u32)) as i16;
        samples.push((val, valid));
    }

    // Pass 2: constant-time fill — process every sample, no branching on validity
    let mut i: usize = 0;
    for &(val, valid) in &samples {
        // Constant-time: only write if valid AND we still need values
        let need = i < n;
        let accept = valid & need;

        // Branchless conditional write using a mask
        let mask = -(accept as i16); // -1 (0xFFFF) if accept, 0 if not
        let idx = if i < n { i } else { n - 1 };
        let old = result[idx].value();
        result[idx] = Zq::from_i16_unchecked((old & !mask) | (val & mask));

        i += accept as usize;
    }

    assert!(i >= n, "hash_to_point: insufficient accepted samples");

    result
}

/// Generates a random nonce for signing.
pub fn generate_nonce<R: rand::RngCore>(rng: &mut R) -> [u8; NONCE_SIZE] {
    let mut nonce = [0u8; NONCE_SIZE];
    rng.fill_bytes(&mut nonce);
    nonce
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::FALCON_512;

    #[test]
    fn test_hash_to_point_deterministic() {
        let message = b"test message";
        let nonce = [0u8; NONCE_SIZE];

        let c1 = hash_to_point(message, &nonce, &FALCON_512);
        let c2 = hash_to_point(message, &nonce, &FALCON_512);

        assert_eq!(c1.len(), 512);
        assert_eq!(c1, c2);
    }

    #[test]
    fn test_hash_to_point_different_inputs() {
        let message1 = b"message 1";
        let message2 = b"message 2";
        let nonce = [0u8; NONCE_SIZE];

        let c1 = hash_to_point(message1, &nonce, &FALCON_512);
        let c2 = hash_to_point(message2, &nonce, &FALCON_512);

        assert_ne!(c1, c2);
    }

    #[test]
    fn test_hash_to_point_range() {
        let message = b"test";
        let nonce = [42u8; NONCE_SIZE];

        let c = hash_to_point(message, &nonce, &FALCON_512);

        for coeff in &c {
            assert!(coeff.value() >= 0 && coeff.value() < Q as i16);
        }
    }

    #[test]
    fn test_generate_nonce() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;

        let mut rng = StdRng::seed_from_u64(42);
        let nonce1 = generate_nonce(&mut rng);
        let nonce2 = generate_nonce(&mut rng);

        assert_eq!(nonce1.len(), NONCE_SIZE);
        assert_ne!(nonce1, nonce2);
    }
}
