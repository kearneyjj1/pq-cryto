//! Verification for ML-DSA.
//!
//! Implements ML-DSA.Verify (FIPS 204 Algorithm 8).

use crate::error::{MlDsaError, Result};
use crate::keygen::PublicKey;
use crate::packing::pack_w1_vec;
use crate::params::D;
use crate::polyvec::PolyVec;
use crate::rounding::use_hint_vec;
use crate::sampling::{expand_a, hash_challenge, hash_message, sample_in_ball};
use crate::sign::Signature;

/// Verifies an ML-DSA signature.
///
/// This implements ML-DSA.Verify (FIPS 204 Algorithm 8).
///
/// # Arguments
/// * `pk` - Public key
/// * `message` - Message that was signed
/// * `sig` - Signature to verify
///
/// # Returns
/// Ok(()) if the signature is valid, Err(InvalidSignature) otherwise.
pub fn verify(pk: &PublicKey, message: &[u8], sig: &Signature) -> Result<()> {
    let params = &pk.params;

    // Check signature dimensions
    if sig.z.len() != params.l {
        return Err(MlDsaError::InvalidSignature);
    }
    if sig.hints.len() != params.k {
        return Err(MlDsaError::InvalidSignature);
    }
    if sig.c_tilde.len() != params.lambda / 4 {
        return Err(MlDsaError::InvalidSignature);
    }

    // Check hint counts
    let hint_count: usize = sig
        .hints
        .iter()
        .map(|h| h.iter().filter(|&&b| b).count())
        .sum();
    if hint_count > params.omega {
        return Err(MlDsaError::InvalidHint);
    }

    // Check ||z||_∞ < γ1 - β
    if !sig.z.check_norm(params.gamma1 - params.beta) {
        return Err(MlDsaError::InvalidSignature);
    }

    // Compute message hash μ = H(tr || M)
    let tr = pk.hash();
    let mu = hash_message(&tr, message);

    // Sample challenge polynomial c from c_tilde
    let mut c = sample_in_ball(&sig.c_tilde, params.tau);
    c.ntt();

    // Generate matrix A from ρ
    let mut a = expand_a(&pk.rho, params);
    a.ntt();

    // Compute Az in NTT domain
    let mut z_ntt = sig.z.clone();
    z_ntt.ntt();
    let az = a.mul_vec(&z_ntt);

    // Compute ct1 * 2^d in NTT domain
    let mut ct1 = PolyVec::zero(params.k);
    for i in 0..params.k {
        let mut t1_shifted = pk.t1.polys[i].clone();
        t1_shifted.shift_left(D);
        t1_shifted.ntt();
        ct1.polys[i] = c.pointwise_mul(&t1_shifted);
    }

    // Compute w' = Az - ct1 * 2^d
    let mut w_prime = PolyVec::zero(params.k);
    for i in 0..params.k {
        for j in 0..256 {
            w_prime.polys[i].coeffs[j] = az.polys[i].coeffs[j] - ct1.polys[i].coeffs[j];
        }
    }
    w_prime.inv_ntt();
    w_prime.reduce();

    // Use hints to recover w1'
    let w1_prime = use_hint_vec(&sig.hints, &w_prime, params.gamma2);

    // Compute c_tilde' = H(μ || w1')
    let w1_prime_bytes = pack_w1_vec(&w1_prime, params.gamma2)
        .expect("invalid gamma2 parameter");
    let c_tilde_prime = hash_challenge(&mu, &w1_prime_bytes, params.lambda);

    // Verify c_tilde' == c_tilde
    if c_tilde_prime != sig.c_tilde {
        return Err(MlDsaError::InvalidSignature);
    }

    Ok(())
}

/// Verifies an ML-DSA signature, returning a boolean.
///
/// This is a convenience wrapper around `verify` that returns true for valid
/// signatures and false for invalid ones.
pub fn verify_bool(pk: &PublicKey, message: &[u8], sig: &Signature) -> bool {
    verify(pk, message, sig).is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::keygen;
    use crate::params::{ML_DSA_44, ML_DSA_65, ML_DSA_87, N};
    use crate::rounding::high_bits_vec;
    use crate::sign::sign;
    use rand::rngs::OsRng;

    #[test]
    fn test_verify_debug_comprehensive() {
        // Test with random seeds to ensure robustness
        use rand::Rng;
        for _ in 0..100 {
            let mut seed = [0u8; 32];
            rand::thread_rng().fill(&mut seed);

            let (pk, sk) = crate::keygen::keygen_internal(&seed, ML_DSA_44);
            let message = b"test";
            let sig = sign(&sk, message).expect("signing should succeed");

            assert!(
                verify(&pk, message, &sig).is_ok(),
                "verification should succeed for random seed"
            );
        }
    }

    #[test]
    fn test_verify_debug() {
        // Use deterministic key for debugging
        let seed = [42u8; 32];
        let (pk, sk) = crate::keygen::keygen_internal(&seed, ML_DSA_44);
        let params = &pk.params;
        let message = b"test";

        // Sign
        let sig = sign(&sk, message).expect("signing should succeed");

        // Now manually trace verification
        let tr = pk.hash();
        let mu = hash_message(&tr, message);

        // Regenerate challenge
        let mut c = sample_in_ball(&sig.c_tilde, params.tau);
        c.ntt();

        // Generate A
        let mut a = expand_a(&pk.rho, params);
        a.ntt();

        // Compute Az
        let mut z_ntt = sig.z.clone();
        z_ntt.ntt();
        let az = a.mul_vec(&z_ntt);

        // Compute ct1 * 2^d
        let mut ct1 = PolyVec::zero(params.k);
        for i in 0..params.k {
            let mut t1_shifted = pk.t1.polys[i].clone();
            t1_shifted.shift_left(D);
            t1_shifted.ntt();
            ct1.polys[i] = c.pointwise_mul(&t1_shifted);
        }

        // w' = Az - ct1*2^d
        let mut w_prime = PolyVec::zero(params.k);
        for i in 0..params.k {
            for j in 0..N {
                w_prime.polys[i].coeffs[j] = az.polys[i].coeffs[j] - ct1.polys[i].coeffs[j];
            }
        }
        w_prime.inv_ntt();
        w_prime.reduce();

        // Get high bits with and without hint
        let w1_no_hint = high_bits_vec(&w_prime, params.gamma2);
        let w1_with_hint = use_hint_vec(&sig.hints, &w_prime, params.gamma2);

        // Count differences
        let mut diff_count = 0;
        for i in 0..params.k {
            for j in 0..N {
                if w1_no_hint.polys[i].coeffs[j] != w1_with_hint.polys[i].coeffs[j] {
                    diff_count += 1;
                }
            }
        }
        println!("Differences from hints: {}", diff_count);

        // Compute hashes
        let w1_prime_bytes = pack_w1_vec(&w1_with_hint, params.gamma2)
            .expect("invalid gamma2 parameter");
        let c_tilde_prime = hash_challenge(&mu, &w1_prime_bytes, params.lambda);

        println!("c_tilde     len={}: {:?}", sig.c_tilde.len(), &sig.c_tilde[..8]);
        println!("c_tilde'    len={}: {:?}", c_tilde_prime.len(), &c_tilde_prime[..8]);
        println!("Match: {}", c_tilde_prime == sig.c_tilde);

        // Check hint count
        let hint_count: usize = sig.hints.iter().map(|h| h.iter().filter(|&&b| b).count()).sum();
        println!("Hint count: {}", hint_count);

        assert!(verify(&pk, message, &sig).is_ok(), "verification should succeed");
    }

    #[test]
    fn test_verify_valid_signature_ml_dsa_44() {
        // Run multiple iterations to ensure consistency
        for _ in 0..100 {
            let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
            let message = b"test";
            let sig = sign(&sk, message).expect("signing should succeed");
            assert!(
                verify(&pk, message, &sig).is_ok(),
                "verification should succeed"
            );
        }
    }

    #[test]
    fn test_verify_valid_signature_ml_dsa_65() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_65);
        let message = b"test message for ML-DSA-65";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_verify_valid_signature_ml_dsa_87() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_87);
        let message = b"test message for ML-DSA-87";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_verify_wrong_message() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);

        let sig = sign(&sk, b"correct message").expect("signing should succeed");
        assert!(verify(&pk, b"wrong message", &sig).is_err());
    }

    #[test]
    fn test_verify_wrong_key() {
        let (pk1, sk1) = keygen(&mut OsRng, ML_DSA_44);
        let (pk2, _sk2) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test message";

        let sig = sign(&sk1, message).expect("signing should succeed");

        // Signature should verify with correct key
        assert!(verify(&pk1, message, &sig).is_ok());

        // Signature should not verify with wrong key
        assert!(verify(&pk2, message, &sig).is_err());
    }

    #[test]
    fn test_verify_tampered_signature_z() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test message";

        let mut sig = sign(&sk, message).expect("signing should succeed");

        // Tamper with z
        sig.z.polys[0].coeffs[0] ^= 1;

        assert!(verify(&pk, message, &sig).is_err());
    }

    #[test]
    fn test_verify_tampered_signature_c_tilde() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test message";

        let mut sig = sign(&sk, message).expect("signing should succeed");

        // Tamper with c_tilde
        sig.c_tilde[0] ^= 1;

        assert!(verify(&pk, message, &sig).is_err());
    }

    #[test]
    fn test_verify_tampered_hints() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test message";

        let mut sig = sign(&sk, message).expect("signing should succeed");

        // Flip a hint bit
        sig.hints[0][0] = !sig.hints[0][0];

        assert!(verify(&pk, message, &sig).is_err());
    }

    #[test]
    fn test_verify_empty_message() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);

        let sig = sign(&sk, b"").expect("signing empty message should succeed");
        assert!(verify(&pk, b"", &sig).is_ok());

        // Empty signature should not verify for non-empty message
        assert!(verify(&pk, b"non-empty", &sig).is_err());
    }

    #[test]
    fn test_verify_long_message() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let long_message = vec![0xAB; 10000];

        let sig = sign(&sk, &long_message).expect("signing long message should succeed");
        assert!(verify(&pk, &long_message, &sig).is_ok());
    }

    #[test]
    fn test_verify_bool() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test";

        let sig = sign(&sk, message).unwrap();

        assert!(verify_bool(&pk, message, &sig));
        assert!(!verify_bool(&pk, b"wrong", &sig));
    }

    #[test]
    fn test_verify_multiple_signatures() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);

        for i in 0..10 {
            let message = format!("message {}", i);
            let sig = sign(&sk, message.as_bytes()).expect("signing should succeed");
            assert!(verify(&pk, message.as_bytes(), &sig).is_ok());
        }
    }

    #[test]
    fn test_verify_wrong_z_length() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test";

        let mut sig = sign(&sk, message).unwrap();
        sig.z.polys.pop(); // Remove one polynomial

        assert!(matches!(
            verify(&pk, message, &sig),
            Err(MlDsaError::InvalidSignature)
        ));
    }

    #[test]
    fn test_verify_wrong_hints_length() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test";

        let mut sig = sign(&sk, message).unwrap();
        sig.hints.pop(); // Remove one hint vector

        assert!(matches!(
            verify(&pk, message, &sig),
            Err(MlDsaError::InvalidSignature)
        ));
    }

    #[test]
    fn test_verify_too_many_hints() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"test";

        let mut sig = sign(&sk, message).unwrap();

        // Set all hints to true (exceeds omega)
        for h in &mut sig.hints {
            for b in h.iter_mut() {
                *b = true;
            }
        }

        assert!(matches!(
            verify(&pk, message, &sig),
            Err(MlDsaError::InvalidHint)
        ));
    }
}
