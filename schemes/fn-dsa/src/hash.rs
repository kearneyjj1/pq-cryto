//! Hash functions for FALCON.
//!
//! This module implements the hash-to-point function using SHAKE256,
//! which maps a message and nonce to a polynomial in Z_q[X]/(X^n + 1).

use sha3::{Shake256, digest::{ExtendableOutput, Update, XofReader}};
use crate::field::Zq;
use crate::params::{Params, Q, NONCE_SIZE};

/// Hashes a message and nonce to a polynomial in Z_q.
///
/// Uses SHAKE256 to generate a uniformly random polynomial c ∈ Z_q[X]/(X^n + 1).
/// The output polynomial has coefficients in [0, q-1].
///
/// Matches the Falcon reference implementation (`hash_to_point_ct`):
/// - Input: SHAKE256(nonce || message)
/// - Reads 2 bytes big-endian as a 16-bit value w
/// - Accepts if w < 5*q = 61445, rejects otherwise
/// - Reduces w mod q to get the coefficient
pub fn hash_to_point(message: &[u8], nonce: &[u8], params: &Params) -> Vec<Zq> {
    let n = params.n;
    let mut result = vec![Zq::ZERO; n];

    // Create SHAKE256 instance
    let mut hasher = Shake256::default();
    hasher.update(nonce);
    hasher.update(message);

    let mut reader = hasher.finalize_xof();

    // Generate n coefficients using rejection sampling
    // (matching Falcon reference: big-endian, accept < 5*q, reduce mod q)
    let mut buf = [0u8; 2];
    let mut i = 0;

    while i < n {
        reader.read(&mut buf);

        // Big-endian interpretation (matching Falcon reference implementation)
        let w = ((buf[0] as u32) << 8) | (buf[1] as u32);

        // Accept if w < 5*q = 61445, reject otherwise
        if w < 61445 {
            let val = (w % (Q as u32)) as i16;
            result[i] = Zq::from_i16_unchecked(val);
            i += 1;
        }
    }

    result
}

/// Computes a domain-separated hash for key generation.
pub fn hash_for_keygen(seed: &[u8], counter: u32, output_len: usize) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(b"FALCON-KEYGEN");
    hasher.update(seed);
    hasher.update(&counter.to_le_bytes());

    let mut reader = hasher.finalize_xof();
    let mut output = vec![0u8; output_len];
    reader.read(&mut output);
    output
}

/// Generates a random nonce for signing.
pub fn generate_nonce<R: rand::RngCore>(rng: &mut R) -> [u8; NONCE_SIZE] {
    let mut nonce = [0u8; NONCE_SIZE];
    rng.fill_bytes(&mut nonce);
    nonce
}

/// Computes the hash used in signature verification.
///
/// This takes the public key hash, nonce, and message and produces
/// the challenge polynomial c.
pub fn hash_challenge(
    pk_hash: &[u8],
    nonce: &[u8],
    message: &[u8],
    params: &Params,
) -> Vec<Zq> {
    let mut hasher = Shake256::default();
    hasher.update(pk_hash);
    hasher.update(nonce);
    hasher.update(message);

    let mut reader = hasher.finalize_xof();

    let n = params.n;
    let mut result = vec![Zq::ZERO; n];
    let mut buf = [0u8; 2];
    let mut i = 0;

    while i < n {
        reader.read(&mut buf);

        let w = ((buf[0] as u32) << 8) | (buf[1] as u32);

        if w < 61445 {
            let val = (w % (Q as u32)) as i16;
            result[i] = Zq::from_i16_unchecked(val);
            i += 1;
        }
    }

    result
}

/// Hashes a public key to a fixed-size digest.
pub fn hash_public_key(h: &[i16]) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(b"FALCON-PK");

    for &coeff in h {
        hasher.update(&coeff.to_le_bytes());
    }

    let mut reader = hasher.finalize_xof();
    let mut output = vec![0u8; 32];
    reader.read(&mut output);
    output
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
    fn test_hash_public_key() {
        let h: Vec<i16> = (0..512).map(|i| (i * 7) as i16 % Q as i16).collect();
        let hash1 = hash_public_key(&h);
        let hash2 = hash_public_key(&h);

        assert_eq!(hash1.len(), 32);
        assert_eq!(hash1, hash2);
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
