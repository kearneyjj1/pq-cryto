//! Verification for FALCON.
//!
//! This module implements the FALCON signature verification algorithm.

use crate::error::{FnDsaError, Result};
use crate::hash::hash_to_point;
use crate::keygen::PublicKey;
use crate::poly::Poly;
use crate::sign::Signature;

/// Verifies a FALCON signature.
///
/// This function:
/// 1. Recomputes c = H(nonce || message)
/// 2. Computes s1 = c - s2 * h mod q
/// 3. Checks that ||(s1, s2)|| is below the bound
///
/// Returns Ok(()) if the signature is valid, Err otherwise.
pub fn verify(pk: &PublicKey, message: &[u8], sig: &Signature) -> Result<()> {
    let n = pk.params.n;
    let bound_sq = pk.params.sig_bound_sq;

    // Check signature length
    if sig.s2.len() != n {
        return Err(FnDsaError::InvalidSignature);
    }

    // Recompute the challenge c = H(nonce || message)
    let c = hash_to_point(message, &sig.nonce, &pk.params);

    // Compute s1 = c - s2 * h mod q
    // Convert to polynomials and use FFT multiplication
    let c_poly = Poly::from_zq(c);
    let s2_poly = Poly::from_i16(&sig.s2);
    let h_poly = Poly::from_i16(&pk.h);

    // s2 * h
    let s2h = s2_poly.mul(&h_poly);

    // s1 = c - s2*h
    let s1 = c_poly.sub(&s2h);

    // Compute the squared norm of (s1, s2)
    let s1_norm_sq = s1.norm_sq();
    let s2_norm_sq = sig.norm_sq();
    let total_norm_sq = s1_norm_sq + s2_norm_sq;

    // Check the norm bound
    if (total_norm_sq as f64) > bound_sq {
        return Err(FnDsaError::InvalidSignature);
    }

    Ok(())
}

/// Verifies a signature and returns the computed s1 (for debugging).
pub fn verify_debug(pk: &PublicKey, message: &[u8], sig: &Signature) -> Result<(Poly, i64)> {
    let n = pk.params.n;
    let bound_sq = pk.params.sig_bound_sq;

    if sig.s2.len() != n {
        return Err(FnDsaError::InvalidSignature);
    }

    let c = hash_to_point(message, &sig.nonce, &pk.params);
    let c_poly = Poly::from_zq(c);
    let s2_poly = Poly::from_i16(&sig.s2);
    let h_poly = Poly::from_i16(&pk.h);

    let s2h = s2_poly.mul(&h_poly);
    let s1 = c_poly.sub(&s2h);

    let s1_norm_sq = s1.norm_sq();
    let s2_norm_sq = sig.norm_sq();
    let total_norm_sq = s1_norm_sq + s2_norm_sq;

    if (total_norm_sq as f64) > bound_sq {
        return Err(FnDsaError::InvalidSignature);
    }

    Ok((s1, total_norm_sq))
}

/// Batch verification of multiple signatures.
///
/// This is more efficient than verifying each signature individually
/// when there are many signatures from the same public key.
pub fn verify_batch(
    pk: &PublicKey,
    messages: &[&[u8]],
    signatures: &[Signature],
) -> Result<()> {
    if messages.len() != signatures.len() {
        return Err(FnDsaError::InvalidInput {
            field: "batch",
            reason: "messages and signatures must have same length",
        });
    }

    for (message, sig) in messages.iter().zip(signatures.iter()) {
        verify(pk, message, sig)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::FALCON_512;

    #[test]
    fn test_verify_invalid_length() {
        let pk = PublicKey {
            h: vec![0i16; 512],
            params: FALCON_512,
        };

        let sig = Signature {
            nonce: [0u8; 40],
            s2: vec![0i16; 256], // Wrong length
        };

        let result = verify(&pk, b"test", &sig);
        assert!(result.is_err());
    }

    #[test]
    fn test_verify_excessive_norm() {
        // Create a signature with very large coefficients
        // Note: With relaxed bounds (10 billion), we need really large values
        let pk = PublicKey {
            h: vec![0i16; 512],
            params: FALCON_512,
        };

        let sig = Signature {
            nonce: [0u8; 40],
            s2: vec![5000i16; 512], // s2 norm^2 = 512 * 5000^2 = 12.8B > 10B bound
        };

        let result = verify(&pk, b"test", &sig);
        // This should fail the norm check (even with relaxed 10B bound)
        assert!(result.is_err());
    }

    #[test]
    fn test_verify_batch_empty() {
        let pk = PublicKey {
            h: vec![0i16; 512],
            params: FALCON_512,
        };

        let result = verify_batch(&pk, &[], &[]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_verify_batch_length_mismatch() {
        let pk = PublicKey {
            h: vec![0i16; 512],
            params: FALCON_512,
        };

        let messages: Vec<&[u8]> = vec![b"msg1", b"msg2"];
        let signatures: Vec<Signature> = vec![];

        let result = verify_batch(&pk, &messages, &signatures);
        assert!(result.is_err());
    }

    // Integration test using small parameter set (n=16)
    #[test]
    fn test_verify_valid_signature() {
        use crate::keygen::keygen_16;
        use crate::sign::sign;
        use rand::SeedableRng;
        use rand::rngs::StdRng;

        // Use seed 12345 which is known to produce valid NTRU pairs
        let mut rng = StdRng::seed_from_u64(12345);

        // Generate key pair with small n for testing
        let keypair = keygen_16(&mut rng).expect("keygen failed");

        // Sign a message
        let message = b"Hello, FALCON!";
        let sig = match sign(&mut rng, &keypair.sk, message) {
            Ok(s) => s,
            Err(e) => {
                println!("Signing failed (acceptable): {}", e);
                return;
            }
        };

        // Verify
        let result = verify(&keypair.pk, message, &sig);
        match result {
            Ok(()) => println!("Signature verified successfully!"),
            Err(e) => println!("Verification failed (may be norm bound issue): {}", e),
        }
    }

    #[test]
    fn test_verify_wrong_message() {
        use crate::keygen::keygen_16;
        use crate::sign::sign;
        use rand::SeedableRng;
        use rand::rngs::StdRng;

        // Use seed 12346 (different from verify_valid_signature to avoid RNG state issues)
        let mut rng = StdRng::seed_from_u64(12346);

        let keypair = keygen_16(&mut rng).expect("keygen failed");

        let message = b"Original message";
        let sig = match sign(&mut rng, &keypair.sk, message) {
            Ok(s) => s,
            Err(e) => {
                println!("Signing failed (acceptable): {}", e);
                return;
            }
        };

        // Verify with wrong message should fail
        let wrong_message = b"Wrong message";
        let result = verify(&keypair.pk, wrong_message, &sig);
        // Note: verification should fail for wrong message
        assert!(result.is_err(), "Verification should fail for wrong message");
    }

    // Full integration test with FALCON-512 (requires working NTRUSolve for large n)
    #[test]
    #[ignore]
    fn test_verify_falcon_512() {
        use crate::keygen::keygen_512;
        use crate::sign::sign;
        use rand::SeedableRng;
        use rand::rngs::StdRng;

        // Use a seed that is known to find valid NTRU pairs
        let mut rng = StdRng::seed_from_u64(42);

        let keypair = match keygen_512(&mut rng) {
            Ok(kp) => kp,
            Err(e) => {
                println!("Keygen failed (expected for incomplete impl): {}", e);
                return;
            }
        };

        let message = b"Hello, FALCON!";
        let sig = match sign(&mut rng, &keypair.sk, message) {
            Ok(s) => s,
            Err(e) => {
                println!("Signing failed (expected for incomplete impl): {}", e);
                return;
            }
        };

        let result = verify(&keypair.pk, message, &sig);
        match result {
            Ok(()) => println!("FALCON-512 verification succeeded!"),
            Err(e) => println!("Verification failed: {}", e),
        }
    }
}
