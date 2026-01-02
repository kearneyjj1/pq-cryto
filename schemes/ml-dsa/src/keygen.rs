//! Key generation for ML-DSA.
//!
//! Implements ML-DSA.KeyGen (FIPS 204 Algorithm 6).

use crate::packing::{pack_eta_vec, pack_t0_vec, pack_t1_vec};
use crate::params::Params;
use crate::polyvec::PolyVec;
use crate::rounding::power2round_vec;
use crate::sampling::{expand_a, expand_s, hash_public_key};
use rand::{CryptoRng, RngCore};

/// ML-DSA public key.
#[derive(Clone, Debug)]
pub struct PublicKey {
    /// Seed for matrix A generation.
    pub rho: [u8; 32],
    /// High bits of t = As1 + s2.
    pub t1: PolyVec,
    /// Parameter set.
    pub params: Params,
}

impl PublicKey {
    /// Serializes the public key to bytes.
    ///
    /// Format: ρ || t1_packed
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.params.public_key_size());
        bytes.extend_from_slice(&self.rho);
        bytes.extend(pack_t1_vec(&self.t1));
        bytes
    }

    /// Returns the public key hash tr = H(pk).
    pub fn hash(&self) -> [u8; 64] {
        hash_public_key(&self.to_bytes())
    }
}

/// ML-DSA secret key.
#[derive(Clone, Debug)]
pub struct SecretKey {
    /// Seed for matrix A generation.
    pub rho: [u8; 32],
    /// Signing key seed.
    pub k: [u8; 32],
    /// Public key hash.
    pub tr: [u8; 64],
    /// Secret vector s1 (l polynomials).
    pub s1: PolyVec,
    /// Secret vector s2 (k polynomials).
    pub s2: PolyVec,
    /// Low bits of t (k polynomials).
    pub t0: PolyVec,
    /// Parameter set.
    pub params: Params,
}

impl SecretKey {
    /// Serializes the secret key to bytes.
    ///
    /// Format: ρ || K || tr || s1_packed || s2_packed || t0_packed
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.params.secret_key_size());
        bytes.extend_from_slice(&self.rho);
        bytes.extend_from_slice(&self.k);
        bytes.extend_from_slice(&self.tr);
        bytes.extend(pack_eta_vec(&self.s1, self.params.eta));
        bytes.extend(pack_eta_vec(&self.s2, self.params.eta));
        bytes.extend(pack_t0_vec(&self.t0));
        bytes
    }
}

/// Generates an ML-DSA key pair.
///
/// This implements ML-DSA.KeyGen (FIPS 204 Algorithm 6).
///
/// # Arguments
/// * `rng` - Cryptographically secure random number generator
/// * `params` - Parameter set (ML_DSA_44, ML_DSA_65, or ML_DSA_87)
///
/// # Returns
/// A tuple of (public_key, secret_key).
pub fn keygen<R: RngCore + CryptoRng>(rng: &mut R, params: Params) -> (PublicKey, SecretKey) {
    // Generate random seeds
    let mut seed = [0u8; 32];
    rng.fill_bytes(&mut seed);

    keygen_internal(&seed, params)
}

/// Internal key generation from seed (for deterministic testing).
pub fn keygen_internal(seed: &[u8; 32], params: Params) -> (PublicKey, SecretKey) {
    use sha3::{
        digest::{ExtendableOutput, Update, XofReader},
        Shake256,
    };

    // Expand seed: (ρ, ρ', K) = H(ξ)
    let mut hasher = Shake256::default();
    hasher.update(seed);
    let mut reader = hasher.finalize_xof();

    let mut rho = [0u8; 32];
    let mut rho_prime = [0u8; 64];
    let mut k = [0u8; 32];

    reader.read(&mut rho);
    reader.read(&mut rho_prime);
    reader.read(&mut k);

    // Generate matrix A from ρ
    let mut a = expand_a(&rho, &params);
    a.ntt();

    // Generate secret vectors s1, s2 from ρ'
    let (mut s1, s2) = expand_s(&rho_prime, &params);

    // Compute t = As1 + s2
    s1.ntt();
    let mut t = a.mul_vec(&s1);
    t.inv_ntt();
    t += &s2;
    t.reduce();

    // Split t into t1 (high bits) and t0 (low bits)
    let (t1, t0) = power2round_vec(&t);

    // Compute public key hash
    let pk = PublicKey {
        rho,
        t1,
        params,
    };
    let tr = pk.hash();

    // Recover s1 from NTT domain for storage
    let mut s1_coeffs = s1.clone();
    s1_coeffs.inv_ntt();

    let sk = SecretKey {
        rho,
        k,
        tr,
        s1: s1_coeffs,
        s2,
        t0,
        params,
    };

    (pk, sk)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ML_DSA_44, ML_DSA_65, ML_DSA_87, D};
    use rand::rngs::OsRng;

    #[test]
    fn test_keygen_ml_dsa_44() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);

        assert_eq!(pk.params.k, 4);
        assert_eq!(pk.params.l, 4);
        assert_eq!(pk.t1.len(), 4);
        assert_eq!(sk.s1.len(), 4);
        assert_eq!(sk.s2.len(), 4);
        assert_eq!(sk.t0.len(), 4);
    }

    #[test]
    fn test_keygen_ml_dsa_65() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_65);

        assert_eq!(pk.params.k, 6);
        assert_eq!(pk.params.l, 5);
        assert_eq!(pk.t1.len(), 6);
        assert_eq!(sk.s1.len(), 5);
        assert_eq!(sk.s2.len(), 6);
        assert_eq!(sk.t0.len(), 6);
    }

    #[test]
    fn test_keygen_ml_dsa_87() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_87);

        assert_eq!(pk.params.k, 8);
        assert_eq!(pk.params.l, 7);
        assert_eq!(pk.t1.len(), 8);
        assert_eq!(sk.s1.len(), 7);
        assert_eq!(sk.s2.len(), 8);
        assert_eq!(sk.t0.len(), 8);
    }

    #[test]
    fn test_keygen_deterministic() {
        let seed = [42u8; 32];

        let (pk1, sk1) = keygen_internal(&seed, ML_DSA_44);
        let (pk2, sk2) = keygen_internal(&seed, ML_DSA_44);

        // Same seed should produce same keys
        assert_eq!(pk1.rho, pk2.rho);
        assert_eq!(sk1.k, sk2.k);
        assert_eq!(pk1.to_bytes(), pk2.to_bytes());
        assert_eq!(sk1.to_bytes(), sk2.to_bytes());
    }

    #[test]
    fn test_keygen_different_seeds() {
        let seed1 = [1u8; 32];
        let seed2 = [2u8; 32];

        let (pk1, _) = keygen_internal(&seed1, ML_DSA_44);
        let (pk2, _) = keygen_internal(&seed2, ML_DSA_44);

        // Different seeds should produce different keys
        assert_ne!(pk1.rho, pk2.rho);
    }

    #[test]
    fn test_public_key_size() {
        let (pk, _) = keygen(&mut OsRng, ML_DSA_44);
        let bytes = pk.to_bytes();

        // pk = ρ (32) + t1 (k * 320)
        assert_eq!(bytes.len(), 32 + 4 * 320);
        assert_eq!(bytes.len(), ML_DSA_44.public_key_size());
    }

    #[test]
    fn test_secret_key_size() {
        let (_, sk) = keygen(&mut OsRng, ML_DSA_44);
        let bytes = sk.to_bytes();

        assert_eq!(bytes.len(), ML_DSA_44.secret_key_size());
    }

    #[test]
    fn test_t1_t0_relationship() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);

        // t1 and t0 should form Power2Round decomposition
        for i in 0..pk.t1.len() {
            for j in 0..256 {
                let t1_coeff = pk.t1.polys[i].coeffs[j];
                let t0_coeff = sk.t0.polys[i].coeffs[j];

                // t1 should be in [0, 2^10)
                assert!(t1_coeff >= 0 && t1_coeff < (1 << 10));

                // t0 should be in [-(2^12), 2^12]
                let half_d = 1 << (D - 1);
                assert!(
                    t0_coeff >= -half_d && t0_coeff <= half_d,
                    "t0[{}][{}] = {} out of bounds",
                    i,
                    j,
                    t0_coeff
                );
            }
        }
    }

    #[test]
    fn test_s1_s2_eta_bounds() {
        use crate::params::Q;

        let (_, sk) = keygen(&mut OsRng, ML_DSA_44);
        let eta = sk.params.eta as i32;

        // s1 and s2 coefficients should be in [-η, η] (mod Q)
        // Negative values are stored as Q - |value|
        for poly in &sk.s1.polys {
            for &c in &poly.coeffs {
                let centered = if c > Q / 2 { c - Q } else { c };
                assert!(
                    centered >= -eta && centered <= eta,
                    "s1 coeff {} (centered: {}) out of bounds",
                    c,
                    centered
                );
            }
        }

        for poly in &sk.s2.polys {
            for &c in &poly.coeffs {
                let centered = if c > Q / 2 { c - Q } else { c };
                assert!(
                    centered >= -eta && centered <= eta,
                    "s2 coeff {} (centered: {}) out of bounds",
                    c,
                    centered
                );
            }
        }
    }

    #[test]
    fn test_tr_is_hash_of_pk() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);

        // tr should equal H(pk)
        let expected_tr = pk.hash();
        assert_eq!(sk.tr, expected_tr);
    }
}
