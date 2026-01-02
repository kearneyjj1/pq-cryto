//! UOV signature generation.
//!
//! This module implements the signing algorithm for the UOV (Unbalanced Oil and
//! Vinegar) signature scheme. The signing process:
//!
//! 1. Hash the message with a random salt to get target values
//! 2. Choose random vinegar variables
//! 3. Build and solve a linear system in the oil variables
//! 4. Transform the solution back through T^(-1)
//!
//! The linear system may be singular for some vinegar choices, so multiple
//! attempts may be needed.

use rand::{CryptoRng, RngCore};

use crate::error::{Result, UovError};
use crate::field::F;
use crate::keygen::{hash_to_field, PublicKey, SecretKey, Signature};
use crate::matrix::{idx_ut, solve};
use crate::params::Params;

/// Default maximum number of signing attempts before returning an error.
pub const DEFAULT_MAX_ATTEMPTS: u32 = 64;

/// Signs a message using the UOV secret key.
///
/// This function attempts to create a valid signature for the given message.
/// It may need multiple attempts if the randomly chosen vinegar values lead
/// to a singular linear system.
///
/// # Arguments
///
/// * `rng` - A cryptographically secure random number generator
/// * `pk` - The public key (used for parameter validation)
/// * `sk` - The secret key
/// * `msg` - The message to sign
///
/// # Returns
///
/// A `Result` containing the signature on success, or an error if signing fails.
///
/// # Errors
///
/// Returns `UovError::SigningFailed` if a valid signature cannot be found after
/// the maximum number of attempts (64 by default).
///
/// Returns `UovError::InvalidKeyPair` if the public and secret keys have
/// mismatched parameters.
///
/// # Example
///
/// ```
/// use rand::rngs::OsRng;
/// use pqsigs_uov::{keygen::keygen, sign::sign, params::PARAMS_DEMO};
///
/// let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
/// let msg = b"Hello, world!";
/// let sig = sign(&mut OsRng, &pk, &sk, msg).expect("signing should succeed");
/// ```
pub fn sign<R: RngCore + CryptoRng>(
    rng: &mut R,
    pk: &PublicKey,
    sk: &SecretKey,
    msg: &[u8],
) -> Result<Signature> {
    sign_with_attempts(rng, pk, sk, msg, DEFAULT_MAX_ATTEMPTS)
}

/// Signs a message with a configurable maximum number of attempts.
///
/// This is useful for testing or when you need different retry behavior.
///
/// # Arguments
///
/// * `rng` - A cryptographically secure random number generator
/// * `pk` - The public key
/// * `sk` - The secret key
/// * `msg` - The message to sign
/// * `max_attempts` - Maximum number of signing attempts
///
/// # Returns
///
/// A `Result` containing the signature on success, or an error if signing fails.
pub fn sign_with_attempts<R: RngCore + CryptoRng>(
    rng: &mut R,
    pk: &PublicKey,
    sk: &SecretKey,
    msg: &[u8],
    max_attempts: u32,
) -> Result<Signature> {
    // Validate key pair consistency
    if pk.params != sk.params {
        return Err(UovError::InvalidKeyPair {
            reason: "public and secret key parameters do not match",
        });
    }

    let Params { v, m } = sk.params;
    let n = v + m;

    // Generate random salt
    let mut salt = vec![0u8; m];
    rng.fill_bytes(&mut salt);

    // Hash message to get target values
    let target = hash_to_field(msg, &salt, m);

    // Try multiple vinegar assignments
    for attempt in 0..max_attempts {
        if let Some(sig) = try_sign_once(rng, sk, &salt, &target, n, v, m) {
            return Ok(sig);
        }

        // If we failed but have more attempts, we'll try again with the same salt
        // This is acceptable as the vinegar values provide sufficient randomness
        if attempt == max_attempts - 1 {
            return Err(UovError::SigningFailed {
                attempts: max_attempts,
            });
        }
    }

    Err(UovError::SigningFailed {
        attempts: max_attempts,
    })
}

/// Attempts a single signing operation with random vinegar values.
fn try_sign_once<R: RngCore + CryptoRng>(
    rng: &mut R,
    sk: &SecretKey,
    salt: &[u8],
    target: &[F],
    n: usize,
    v: usize,
    m: usize,
) -> Option<Signature> {
    // Choose random vinegar values
    let mut vinegar = vec![F::ZERO; n];
    for i in 0..v {
        vinegar[i] = F(rng.next_u32() as u8);
    }

    // Build the linear system A_oil * y_oil = b
    // where y is the internal representation before applying T^(-1)
    let mut a_mat = vec![vec![F::ZERO; m]; m];
    let mut b_vec = target.to_vec();

    for eq in 0..m {
        let ut = &sk.f_quads[eq];

        // Compute vinegar×vinegar contribution (constant term)
        let mut const_term = F::ZERO;
        for i in 0..v {
            for j in i..v {
                let aij = ut[idx_ut(n, i, j)];
                const_term += aij * vinegar[i] * vinegar[j];
            }
        }

        // Move constant to RHS: b = target - const_term
        b_vec[eq] -= const_term;

        // Compute vinegar×oil coefficients (linear in oil variables)
        for i in 0..v {
            for j in v..n {
                let aij = ut[idx_ut(n, i, j)];
                let coeff = aij * vinegar[i];
                a_mat[eq][j - v] += coeff;
            }
        }

        // Note: oil×oil block is zero by construction in central OV map
    }

    // Solve for oil variables
    let y_oil = solve(a_mat, b_vec)?;

    // Compose y = [vinegar | oil]
    let mut y = vec![F::ZERO; n];
    for i in 0..v {
        y[i] = vinegar[i];
    }
    for j in 0..m {
        y[v + j] = y_oil[j];
    }

    // Compute x = T^(-1) * y
    let x = solve_t_inverse(&sk.t, &y)?;

    Some(Signature {
        x,
        salt: salt.to_vec(),
    })
}

/// Solves T * x = y for x using Gaussian elimination.
fn solve_t_inverse(t: &[Vec<F>], y: &[F]) -> Option<Vec<F>> {
    solve(t.to_vec(), y.to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::keygen;
    use crate::params::PARAMS_DEMO;
    use rand::rngs::OsRng;

    #[test]
    fn test_sign_produces_valid_length() {
        let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
        let msg = b"test message";
        let sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

        assert_eq!(sig.x.len(), pk.params.n());
        assert_eq!(sig.salt.len(), pk.params.m);
    }

    #[test]
    fn test_sign_different_messages_different_signatures() {
        let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
        let sig1 = sign(&mut OsRng, &pk, &sk, b"msg1").unwrap();
        let sig2 = sign(&mut OsRng, &pk, &sk, b"msg2").unwrap();

        // Signatures should differ (with overwhelming probability)
        assert!(sig1.x != sig2.x || sig1.salt != sig2.salt);
    }
}
