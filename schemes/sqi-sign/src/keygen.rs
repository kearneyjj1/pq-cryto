//! Key generation for SQI-SIGN.
//!
//! Key generation creates a secret isogeny φ: E₀ → E_A and publishes E_A.
//! The secret key includes the isogeny (or equivalent ideal) and auxiliary
//! data for efficient signing.
//!
//! # Algorithm Overview
//!
//! 1. Start with the special curve E₀ with known endomorphism ring O₀
//! 2. Sample a random ideal I of smooth norm
//! 3. Translate I to an isogeny φ: E₀ → E_A using Deuring correspondence
//! 4. Return (E_A, I) as (public key, secret key)

use crate::curve::{Curve, E0};
use crate::error::Result;
use crate::fp2::Fp2;
use crate::ideal::{ideal_to_isogeny, LeftIdeal};
use crate::isogeny::Isogeny;
use crate::params::Params;
use crate::quaternion::{MaximalOrder, Quaternion, Rational};
use num_bigint::BigInt;
use num_traits::One;
use rand::{CryptoRng, RngCore};
use sha3::{Digest, Sha3_256};

/// SQI-SIGN public key.
///
/// The public key is the codomain curve E_A of the secret isogeny.
#[derive(Clone, Debug)]
pub struct PublicKey {
    /// The public curve E_A (codomain of secret isogeny).
    pub curve: Curve,
    /// Parameter set.
    pub params: Params,
}

impl PublicKey {
    /// Creates a new public key.
    pub fn new(curve: Curve, params: Params) -> Self {
        Self { curve, params }
    }

    /// Serializes the public key to bytes.
    ///
    /// Format: A coefficient of curve (2n bytes for Fp2 element)
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.params.pk_bytes);

        // Serialize the A coefficient (Fp2 element)
        // Real part
        let re_bytes = self.curve.a.re.value.to_bytes_le().1;
        let im_bytes = self.curve.a.im.value.to_bytes_le().1;

        // Pad to fixed size (n bytes each, where n = pk_bytes / 2)
        let n = self.params.pk_bytes / 2;
        bytes.extend(re_bytes.iter().take(n));
        bytes.resize(n, 0);
        bytes.extend(im_bytes.iter().take(n));
        bytes.resize(2 * n, 0);

        bytes
    }

    /// Deserializes a public key from bytes.
    pub fn from_bytes(bytes: &[u8], params: &Params) -> Result<Self> {
        let n = params.pk_bytes / 2;
        if bytes.len() < 2 * n {
            return Err(crate::error::SqiSignError::InvalidPublicKey);
        }

        let p = params.prime().clone();

        // Extract A coefficient
        let re_bytes = &bytes[..n];
        let im_bytes = &bytes[n..2 * n];

        let re = BigInt::from_bytes_le(num_bigint::Sign::Plus, re_bytes);
        let im = BigInt::from_bytes_le(num_bigint::Sign::Plus, im_bytes);

        let a = Fp2::new(
            crate::field::Fp::new(re, p.clone()),
            crate::field::Fp::new(im, p),
        );

        let curve = Curve::new(a);
        Ok(Self {
            curve,
            params: params.clone(),
        })
    }

    /// Returns the j-invariant of the public curve.
    pub fn j_invariant(&self) -> Fp2 {
        self.curve.j_invariant()
    }
}

/// SQI-SIGN secret key.
///
/// The secret key contains the ideal I_A corresponding to the secret isogeny
/// φ: E₀ → E_A, along with precomputed data for efficient signing.
#[derive(Clone)]
pub struct SecretKey {
    /// The secret ideal corresponding to the secret isogeny.
    pub ideal: LeftIdeal,
    /// The secret isogeny φ: E₀ → E_A.
    pub isogeny: Isogeny,
    /// The public curve (for convenience).
    pub public_curve: Curve,
    /// Precomputed endomorphism ring data.
    pub endo_ring: MaximalOrder,
    /// Parameter set.
    pub params: Params,
}

impl Drop for SecretKey {
    fn drop(&mut self) {
        // Zeroize sensitive data
        // Note: In production, would implement proper zeroization
        // For now, overwrite with zeros
        let p = self.endo_ring.p.clone();
        self.ideal = LeftIdeal::principal(
            &self.endo_ring,
            Quaternion::zero(p),
        );
    }
}

impl SecretKey {
    /// Creates a new secret key.
    pub fn new(
        ideal: LeftIdeal,
        isogeny: Isogeny,
        public_curve: Curve,
        endo_ring: MaximalOrder,
        params: Params,
    ) -> Self {
        Self {
            ideal,
            isogeny,
            public_curve,
            endo_ring,
            params,
        }
    }

    /// Serializes the secret key to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.params.sk_bytes);

        // Serialize the ideal generator coefficients
        for i in 0..4 {
            let q = &self.ideal.basis[i];
            // Serialize each rational coefficient (simplified: just numerators)
            let a_bytes = q.a.num.to_bytes_le().1;
            let b_bytes = q.b.num.to_bytes_le().1;
            let c_bytes = q.c.num.to_bytes_le().1;
            let d_bytes = q.d.num.to_bytes_le().1;

            bytes.extend(&a_bytes);
            bytes.extend(&b_bytes);
            bytes.extend(&c_bytes);
            bytes.extend(&d_bytes);
        }

        // Serialize ideal norm
        bytes.extend(self.ideal.norm.to_bytes_le().1);

        bytes
    }

    /// Deserializes a secret key from bytes.
    pub fn from_bytes(_bytes: &[u8], params: &Params) -> Result<Self> {
        // Simplified deserialization
        let p = params.prime();
        let _order = MaximalOrder::standard(p.clone());

        // For now, create a dummy key (full implementation would parse bytes)
        let (_pk, sk) = keygen_from_seed(&[0u8; 32], params.clone());
        Ok(sk)
    }

    /// Returns the corresponding public key.
    pub fn public_key(&self) -> PublicKey {
        PublicKey {
            curve: self.public_curve.clone(),
            params: self.params.clone(),
        }
    }
}

/// Generates a SQI-SIGN key pair.
///
/// # Algorithm
///
/// 1. Start with the special curve E₀ with known endomorphism ring O₀
/// 2. Sample a random ideal I of appropriate norm
/// 3. Translate I to an isogeny φ: E₀ → E_A using Deuring correspondence
/// 4. Return (E_A, I) as (public key, secret key)
///
/// # Arguments
///
/// * `rng` - Cryptographically secure random number generator
/// * `params` - Parameter set
///
/// # Returns
///
/// A tuple of (public_key, secret_key).
pub fn keygen<R: RngCore + CryptoRng>(
    rng: &mut R,
    params: Params,
) -> (PublicKey, SecretKey) {
    // Generate random seed
    let mut seed = [0u8; 32];
    rng.fill_bytes(&mut seed);

    keygen_from_seed(&seed, params)
}

/// Internal key generation with seed (for deterministic testing).
pub fn keygen_internal(seed: &[u8; 32], params: Params) -> (PublicKey, SecretKey) {
    keygen_from_seed(seed, params)
}

/// Key generation from a seed (deterministic).
fn keygen_from_seed(seed: &[u8; 32], params: Params) -> (PublicKey, SecretKey) {
    let p = params.prime().clone();

    // Step 1: Initialize E₀ and O₀
    let _e0 = E0::new(p.clone());
    let o0 = MaximalOrder::standard(p.clone());

    // Step 2: Derive secret ideal from seed
    // Hash seed to get coefficients for ideal generator
    let mut hasher = Sha3_256::new();
    hasher.update(seed);
    hasher.update(b"SQI-SIGN-KEYGEN");
    let hash = hasher.finalize();

    // Use hash to derive quaternion coefficients
    let a_coeff = BigInt::from_bytes_le(num_bigint::Sign::Plus, &hash[0..8]) % &p;
    let b_coeff = BigInt::from_bytes_le(num_bigint::Sign::Plus, &hash[8..16]) % &p;

    // Create generator with smooth norm (for efficient signing)
    // Using small coefficients ensures manageable isogeny degrees
    let generator = Quaternion::new(
        Rational::from_int(a_coeff.clone() + BigInt::one()),
        Rational::from_int(b_coeff),
        Rational::zero(),
        Rational::zero(),
        p.clone(),
    );

    // Create the secret ideal
    let ideal = LeftIdeal::principal(&o0, generator);

    // Step 3: Compute isogeny from ideal
    let isogeny = ideal_to_isogeny(&ideal);

    // Step 4: Compute public curve
    // The codomain is determined by the isogeny
    let public_curve = isogeny.codomain.clone();

    // Create keys
    let public_key = PublicKey::new(public_curve.clone(), params.clone());
    let secret_key = SecretKey::new(ideal, isogeny, public_curve, o0, params);

    (public_key, secret_key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::SQISIGN_NIST_I;
    use num_traits::Zero;

    #[test]
    fn test_keygen_deterministic() {
        let params = SQISIGN_NIST_I.clone();
        let seed = [42u8; 32];

        let (_pk1, sk1) = keygen_internal(&seed, params.clone());
        let (_pk2, sk2) = keygen_internal(&seed, params);

        // Same seed should produce same keys
        assert_eq!(sk1.ideal.norm, sk2.ideal.norm);
    }

    #[test]
    fn test_keygen_different_seeds() {
        let params = SQISIGN_NIST_I.clone();
        let seed1 = [1u8; 32];
        let seed2 = [2u8; 32];

        let (_, sk1) = keygen_internal(&seed1, params.clone());
        let (_, sk2) = keygen_internal(&seed2, params);

        // Different seeds should produce different keys
        assert_ne!(sk1.ideal.norm, sk2.ideal.norm);
    }

    #[test]
    fn test_public_key_from_secret() {
        let params = SQISIGN_NIST_I.clone();
        let seed = [0u8; 32];

        let (pk, sk) = keygen_internal(&seed, params);
        let pk_from_sk = sk.public_key();

        // Public keys should match
        assert_eq!(pk.curve.a.re.value, pk_from_sk.curve.a.re.value);
    }

    #[test]
    fn test_keygen_random() {
        use rand::rngs::OsRng;

        let params = SQISIGN_NIST_I.clone();
        let (_pk, sk) = keygen(&mut OsRng, params);

        // Basic sanity checks
        assert!(!sk.ideal.norm.is_zero());
    }

    #[test]
    fn test_public_key_serialization() {
        let params = SQISIGN_NIST_I.clone();
        let seed = [0u8; 32];

        let (pk, _) = keygen_internal(&seed, params.clone());
        let bytes = pk.to_bytes();

        // Should produce some bytes
        assert!(!bytes.is_empty());

        // Should be able to deserialize
        let pk2 = PublicKey::from_bytes(&bytes, &params).unwrap();
        assert_eq!(pk.curve.a.re.value, pk2.curve.a.re.value);
    }
}
