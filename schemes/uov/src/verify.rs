//! UOV signature verification.
//!
//! This module implements the verification algorithm for the UOV (Unbalanced Oil
//! and Vinegar) signature scheme. Verification works by:
//!
//! 1. Hashing the message with the signature's salt to get target values
//! 2. Evaluating the public quadratic forms at the signature point
//! 3. Checking that the evaluation matches the target values

use crate::error::{Result, UovError};
use crate::field::F;
use crate::keygen::{hash_to_field, PublicKey, QuadForm, Signature};
use crate::matrix::idx_ut;

/// Evaluates a quadratic form at a given point.
///
/// Computes Q(x) = sum_{i <= j} q_ij * x_i * x_j
fn eval_quad(q: &QuadForm, x: &[F]) -> F {
    let n = x.len();
    let mut acc = F::ZERO;

    for i in 0..n {
        for j in i..n {
            let c = q.coeffs[idx_ut(n, i, j)];
            acc += c * x[i] * x[j];
        }
    }

    acc
}

/// Evaluates all public quadratic forms at a given point.
///
/// Returns a vector of m field elements, one for each form.
fn eval_public(pk: &PublicKey, x: &[F]) -> Vec<F> {
    pk.polys.iter().map(|q| eval_quad(q, x)).collect()
}

/// Verifies a UOV signature.
///
/// This function checks that the signature is valid for the given message
/// under the provided public key.
///
/// # Arguments
///
/// * `pk` - The public key
/// * `msg` - The message that was signed
/// * `sig` - The signature to verify
///
/// # Returns
///
/// `Ok(())` if the signature is valid, or `Err(UovError::InvalidSignature)` if not.
///
/// # Errors
///
/// Returns `UovError::InvalidSignature` if the signature fails verification,
/// including cases where:
/// - The signature length doesn't match the parameters
/// - The salt length doesn't match the parameters
/// - The quadratic form evaluations don't match the hash output
///
/// # Example
///
/// ```
/// use rand::rngs::OsRng;
/// use pqsigs_uov::{keygen::keygen, sign::sign, verify::verify, params::PARAMS_DEMO};
///
/// let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
/// let msg = b"Hello, world!";
/// let sig = sign(&mut OsRng, &pk, &sk, msg).expect("signing should succeed");
///
/// assert!(verify(&pk, msg, &sig).is_ok());
/// ```
pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature) -> Result<()> {
    // Validate signature dimensions
    if sig.x.len() != pk.params.n() {
        return Err(UovError::InvalidInput {
            field: "signature",
            reason: "signature point has wrong length",
        });
    }

    if sig.salt.len() != pk.params.m {
        return Err(UovError::InvalidInput {
            field: "signature",
            reason: "salt has wrong length",
        });
    }

    // Hash the message to get target values
    let target = hash_to_field(msg, &sig.salt, pk.params.m);

    // Evaluate public forms at signature point
    let evaluation = eval_public(pk, &sig.x);

    // Check if evaluation matches target
    if target == evaluation {
        Ok(())
    } else {
        Err(UovError::InvalidSignature)
    }
}

/// Verifies a signature, returning a boolean for convenience.
///
/// This is a convenience wrapper around [`verify`] that returns `true` for
/// valid signatures and `false` otherwise.
///
/// # Arguments
///
/// * `pk` - The public key
/// * `msg` - The message that was signed
/// * `sig` - The signature to verify
///
/// # Returns
///
/// `true` if the signature is valid, `false` otherwise.
pub fn verify_bool(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    verify(pk, msg, sig).is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::keygen;
    use crate::params::PARAMS_DEMO;
    use crate::sign::sign;
    use rand::rngs::OsRng;

    #[test]
    fn test_verify_valid_signature() {
        let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
        let msg = b"test message";
        let sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

        assert!(verify(&pk, msg, &sig).is_ok());
    }

    #[test]
    fn test_verify_wrong_message() {
        let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
        let msg = b"test message";
        let sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

        let result = verify(&pk, b"wrong message", &sig);
        assert!(matches!(result, Err(UovError::InvalidSignature)));
    }

    #[test]
    fn test_verify_wrong_signature_length() {
        let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
        let msg = b"test message";
        let mut sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

        // Truncate signature
        sig.x.pop();

        let result = verify(&pk, msg, &sig);
        assert!(matches!(
            result,
            Err(UovError::InvalidInput {
                field: "signature",
                ..
            })
        ));
    }

    #[test]
    fn test_verify_bool() {
        let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
        let msg = b"test message";
        let sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

        assert!(verify_bool(&pk, msg, &sig));
        assert!(!verify_bool(&pk, b"wrong message", &sig));
    }
}
