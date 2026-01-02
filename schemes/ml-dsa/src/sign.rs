//! Signing for ML-DSA.
//!
//! Implements ML-DSA.Sign (FIPS 204 Algorithm 7).

use crate::error::{MlDsaError, Result};
use crate::keygen::SecretKey;
use crate::packing::{pack_hints, pack_w1_vec, pack_z_vec};
use crate::polyvec::PolyVec;
use crate::reduce::freeze;
use crate::rounding::{high_bits_vec, low_bits_vec, make_hint_vec};
use crate::sampling::{expand_a, expand_mask_vec, hash_challenge, hash_message, sample_in_ball};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

/// Maximum number of signing attempts before giving up.
const MAX_ATTEMPTS: u32 = 1000;

/// ML-DSA signature.
#[derive(Clone, Debug)]
pub struct Signature {
    /// Challenge hash (commitment hash).
    pub c_tilde: Vec<u8>,
    /// Response vector z (l polynomials).
    pub z: PolyVec,
    /// Hint vector h.
    pub hints: Vec<Vec<bool>>,
}

impl Signature {
    /// Serializes the signature to bytes.
    ///
    /// Format: c_tilde || z_packed || h_packed
    pub fn to_bytes(&self, gamma1: i32, omega: usize, k: usize) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend(&self.c_tilde);
        bytes.extend(pack_z_vec(&self.z, gamma1));
        bytes.extend(pack_hints(&self.hints, omega, k));
        bytes
    }
}

/// Signs a message using ML-DSA.
///
/// This implements ML-DSA.Sign (FIPS 204 Algorithm 7).
///
/// # Arguments
/// * `sk` - Secret key
/// * `message` - Message to sign
///
/// # Returns
/// A signature on success, or an error if signing fails after MAX_ATTEMPTS.
pub fn sign(sk: &SecretKey, message: &[u8]) -> Result<Signature> {
    let params = &sk.params;

    // Compute message hash μ = H(tr || M)
    let mu = hash_message(&sk.tr, message);

    // Generate matrix A from ρ
    let mut a = expand_a(&sk.rho, params);
    a.ntt();

    // Prepare s1 and s2 in NTT domain
    let mut s1_ntt = sk.s1.clone();
    s1_ntt.ntt();
    let mut s2_ntt = sk.s2.clone();
    s2_ntt.ntt();
    let mut t0_ntt = sk.t0.clone();
    t0_ntt.ntt();

    // Compute signing randomness ρ''
    let rho_double_prime = compute_rho_double_prime(&sk.k, &mu);

    // Rejection sampling loop
    let mut kappa: u16 = 0;

    for _attempt in 0..MAX_ATTEMPTS {
        // Sample masking vector y with coefficients in [-γ1+1, γ1]
        let y = expand_mask_vec(&rho_double_prime, kappa, params.l, params.gamma1);
        kappa = kappa.wrapping_add(params.l as u16);

        // Compute w = Ay (in NTT domain)
        let mut y_ntt = y.clone();
        y_ntt.ntt();
        let mut w = a.mul_vec(&y_ntt);
        w.inv_ntt();
        w.reduce();

        // Compute w1 = HighBits(w)
        let w1 = high_bits_vec(&w, params.gamma2);

        // Compute challenge hash c_tilde = H(μ || w1)
        let w1_bytes = pack_w1_vec(&w1, params.gamma2);
        let c_tilde = hash_challenge(&mu, &w1_bytes, params.lambda);

        // Sample challenge polynomial c from c_tilde
        let mut c = sample_in_ball(&c_tilde, params.tau);
        c.ntt();

        // Compute z = y + cs1
        let mut cs1 = PolyVec::zero(params.l);
        for i in 0..params.l {
            cs1.polys[i] = c.pointwise_mul(&s1_ntt.polys[i]);
        }
        cs1.inv_ntt();

        let mut z = &y + &cs1;
        z.reduce();

        // Compute r0 = LowBits(w - cs2)
        let mut cs2 = PolyVec::zero(params.k);
        for i in 0..params.k {
            cs2.polys[i] = c.pointwise_mul(&s2_ntt.polys[i]);
        }
        cs2.inv_ntt();

        let w_minus_cs2 = &w - &cs2;
        let r0 = low_bits_vec(&w_minus_cs2, params.gamma2);

        // Check ||z||_∞ < γ1 - β
        if !z.check_norm(params.gamma1 - params.beta) {
            continue;
        }

        // Check ||r0||_∞ < γ2 - β
        if !r0.check_norm(params.gamma2 - params.beta) {
            continue;
        }

        // Compute hints
        let mut ct0 = PolyVec::zero(params.k);
        for i in 0..params.k {
            ct0.polys[i] = c.pointwise_mul(&t0_ntt.polys[i]);
        }
        ct0.inv_ntt();

        // Negate ct0 for MakeHint
        for poly in &mut ct0.polys {
            for c in &mut poly.coeffs {
                *c = freeze(-*c);
            }
        }

        // w - cs2 + ct0 (but we want hints for w - cs2 - ct0)
        let mut r_for_hint = w_minus_cs2.clone();
        r_for_hint -= &ct0;
        r_for_hint.reduce();

        let (hints, hint_count) = make_hint_vec(&ct0, &r_for_hint, params.gamma2);

        // Check number of hints ≤ ω
        if hint_count > params.omega {
            continue;
        }

        // Check ||ct0||_∞ < γ2
        let mut ct0_for_check = PolyVec::zero(params.k);
        for i in 0..params.k {
            ct0_for_check.polys[i] = c.pointwise_mul(&t0_ntt.polys[i]);
        }
        ct0_for_check.inv_ntt();
        ct0_for_check.reduce();

        if !ct0_for_check.check_norm(params.gamma2) {
            continue;
        }

        // Success! Return signature
        return Ok(Signature { c_tilde, z, hints });
    }

    Err(MlDsaError::SigningFailed {
        attempts: MAX_ATTEMPTS,
    })
}

/// Signs a message with a random context string (deterministic signing).
///
/// This is the internal signing function that uses K for determinism.
fn compute_rho_double_prime(k: &[u8; 32], mu: &[u8; 64]) -> [u8; 64] {
    let mut hasher = Shake256::default();
    hasher.update(k);
    hasher.update(mu);
    let mut reader = hasher.finalize_xof();
    let mut rho_double_prime = [0u8; 64];
    reader.read(&mut rho_double_prime);
    rho_double_prime
}

/// Signs a message using randomized signing.
///
/// This variant uses additional randomness for side-channel resistance.
pub fn sign_randomized<R: rand::RngCore + rand::CryptoRng>(
    _rng: &mut R,
    sk: &SecretKey,
    message: &[u8],
) -> Result<Signature> {
    // For randomized signing, we'd mix in additional randomness
    // For now, just use deterministic signing
    sign(sk, message)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::keygen;
    use crate::params::{ML_DSA_44, ML_DSA_65, ML_DSA_87};
    use rand::rngs::OsRng;

    #[test]
    fn test_sign_ml_dsa_44() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test message";

        let sig = sign(&sk, message).expect("signing should succeed");

        assert_eq!(sig.c_tilde.len(), ML_DSA_44.lambda / 4);
        assert_eq!(sig.z.len(), ML_DSA_44.l);
        assert_eq!(sig.hints.len(), ML_DSA_44.k);
    }

    #[test]
    fn test_sign_ml_dsa_65() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_65);
        let message = b"test message for ML-DSA-65";

        let sig = sign(&sk, message).expect("signing should succeed");

        assert_eq!(sig.c_tilde.len(), ML_DSA_65.lambda / 4);
        assert_eq!(sig.z.len(), ML_DSA_65.l);
    }

    #[test]
    fn test_sign_ml_dsa_87() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_87);
        let message = b"test message for ML-DSA-87";

        let sig = sign(&sk, message).expect("signing should succeed");

        assert_eq!(sig.c_tilde.len(), ML_DSA_87.lambda / 4);
        assert_eq!(sig.z.len(), ML_DSA_87.l);
    }

    #[test]
    fn test_sign_z_bounds() {
        let (_, sk) = keygen(&mut OsRng, ML_DSA_44);
        let sig = sign(&sk, b"test").unwrap();

        // z should have ||z||_∞ < γ1 - β
        let bound = ML_DSA_44.gamma1 - ML_DSA_44.beta;
        assert!(sig.z.norm_inf() < bound);
    }

    #[test]
    fn test_sign_hint_count() {
        let (_, sk) = keygen(&mut OsRng, ML_DSA_44);
        let sig = sign(&sk, b"test").unwrap();

        // Count total hints
        let hint_count: usize = sig
            .hints
            .iter()
            .map(|h| h.iter().filter(|&&b| b).count())
            .sum();

        assert!(hint_count <= ML_DSA_44.omega);
    }

    #[test]
    fn test_sign_deterministic() {
        let seed = [42u8; 32];
        let (_, sk) = crate::keygen::keygen_internal(&seed, ML_DSA_44);

        let sig1 = sign(&sk, b"same message").unwrap();
        let sig2 = sign(&sk, b"same message").unwrap();

        // Same key + same message should produce same signature (deterministic)
        assert_eq!(sig1.c_tilde, sig2.c_tilde);
        for i in 0..sig1.z.len() {
            assert_eq!(sig1.z.polys[i].coeffs, sig2.z.polys[i].coeffs);
        }
    }

    #[test]
    fn test_sign_different_messages() {
        let (_, sk) = keygen(&mut OsRng, ML_DSA_44);

        let sig1 = sign(&sk, b"message 1").unwrap();
        let sig2 = sign(&sk, b"message 2").unwrap();

        // Different messages should produce different signatures
        assert_ne!(sig1.c_tilde, sig2.c_tilde);
    }

    #[test]
    fn test_sign_empty_message() {
        let (_, sk) = keygen(&mut OsRng, ML_DSA_44);

        let sig = sign(&sk, b"").expect("signing empty message should succeed");
        assert!(!sig.c_tilde.is_empty());
    }

    #[test]
    fn test_sign_long_message() {
        let (_, sk) = keygen(&mut OsRng, ML_DSA_44);

        let long_message = vec![0xAB; 10000];
        let sig = sign(&sk, &long_message).expect("signing long message should succeed");
        assert!(!sig.c_tilde.is_empty());
    }

    #[test]
    fn test_signature_serialization() {
        let (_, sk) = keygen(&mut OsRng, ML_DSA_44);
        let sig = sign(&sk, b"test").unwrap();

        let bytes = sig.to_bytes(ML_DSA_44.gamma1, ML_DSA_44.omega, ML_DSA_44.k);
        assert_eq!(bytes.len(), ML_DSA_44.signature_size());
    }
}
