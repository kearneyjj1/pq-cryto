//! UOV key generation.
//!
//! This module implements key generation for the UOV (Unbalanced Oil and Vinegar)
//! signature scheme. The key generation process:
//!
//! 1. Samples random central OV quadratic forms F (with no oil×oil terms)
//! 2. Samples a random invertible linear transformation matrix T
//! 3. Computes public polynomials P(x) = F(Tx)
//!
//! The secret key contains T and the central forms F.
//! The public key contains the composed quadratic forms P.

use rand::{CryptoRng, RngCore};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
// Note: We manually zeroize F elements since F doesn't implement Zeroize

use crate::field::F;
use crate::matrix::{idx_ut, sample_invertible};
use crate::params::Params;

/// A quadratic form over GF(256), stored in upper triangular format.
///
/// For n variables, the coefficients are stored as a vector of length n*(n+1)/2,
/// representing the upper triangular part of the symmetric coefficient matrix.
#[derive(Clone, Debug)]
pub struct QuadForm {
    /// Coefficients in upper triangular packed format.
    pub coeffs: Vec<F>,
}

impl QuadForm {
    /// Creates a new quadratic form with the given coefficients.
    pub fn new(coeffs: Vec<F>) -> Self {
        QuadForm { coeffs }
    }

    /// Creates a zero quadratic form for n variables.
    pub fn zero(n: usize) -> Self {
        QuadForm {
            coeffs: vec![F::ZERO; n * (n + 1) / 2],
        }
    }
}

/// Public key for the UOV signature scheme.
///
/// Contains the parameter set and m quadratic polynomials that define
/// the public map.
#[derive(Clone, Debug)]
pub struct PublicKey {
    /// The parameter set used for this key.
    pub params: Params,
    /// The m public quadratic forms.
    pub polys: Vec<QuadForm>,
}

impl PublicKey {
    /// Returns the total number of variables n.
    #[inline]
    pub fn n(&self) -> usize {
        self.params.n()
    }

    /// Returns the number of oil variables / equations m.
    #[inline]
    pub fn m(&self) -> usize {
        self.params.m
    }
}

/// Secret key for the UOV signature scheme.
///
/// Contains the parameter set, the invertible transformation matrix T,
/// and the central OV quadratic forms.
///
/// # Security
///
/// This struct implements `Drop` to zeroize secret key material when dropped.
pub struct SecretKey {
    /// The parameter set used for this key.
    pub params: Params,
    /// The invertible linear map T (n × n matrix).
    pub t: Vec<Vec<F>>,
    /// The central OV quadratic forms (m forms, each in UT packed format).
    pub f_quads: Vec<Vec<F>>,
}

impl Drop for SecretKey {
    fn drop(&mut self) {
        // Zeroize the transformation matrix T
        // We manually zero each element since F doesn't implement Zeroize
        for row in &mut self.t {
            for elem in row.iter_mut() {
                *elem = F::ZERO;
            }
        }

        // Zeroize the central quadratic forms
        for form in &mut self.f_quads {
            for elem in form.iter_mut() {
                *elem = F::ZERO;
            }
        }
    }
}

impl SecretKey {
    /// Returns the total number of variables n.
    #[inline]
    pub fn n(&self) -> usize {
        self.params.n()
    }

    /// Returns the number of oil variables / equations m.
    #[inline]
    pub fn m(&self) -> usize {
        self.params.m
    }
}

/// A signature produced by the UOV scheme.
///
/// Contains the signature point x and the random salt used during signing.
#[derive(Clone, Debug)]
pub struct Signature {
    /// The signature point (n field elements).
    pub x: Vec<F>,
    /// The random salt used in hashing (m bytes).
    pub salt: Vec<u8>,
}

/// Domain separator for the hash function.
const DOMAIN_SEPARATOR: &[u8] = b"pqsigs-uov-v1";

/// Hashes a message with a salt to produce m field elements.
///
/// Uses SHAKE256 with a domain separator for deterministic hashing.
pub fn hash_to_field(msg: &[u8], salt: &[u8], m: usize) -> Vec<F> {
    let mut h = Shake256::default();
    h.update(DOMAIN_SEPARATOR);
    h.update(msg);
    h.update(salt);
    let mut xof = h.finalize_xof();
    let mut out = vec![0u8; m];
    xof.read(&mut out);
    out.into_iter().map(F).collect()
}

/// Generates a UOV key pair.
///
/// # Arguments
///
/// * `rng` - A cryptographically secure random number generator
/// * `params` - The parameter set to use
///
/// # Returns
///
/// A tuple containing the public key and secret key.
///
/// # Example
///
/// ```
/// use rand::rngs::OsRng;
/// use pqsigs_uov::{keygen::keygen, params::PARAMS_DEMO};
///
/// let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
/// ```
pub fn keygen<R: RngCore + CryptoRng>(rng: &mut R, params: Params) -> (PublicKey, SecretKey) {
    let n = params.n();
    let ut_len = params.ut_size();
    let v = params.v;
    let m = params.m;

    // Sample central OV quadratic forms (no oil×oil terms)
    let mut f_quads = Vec::with_capacity(m);
    for _ in 0..m {
        let mut a = vec![F::ZERO; ut_len];
        for i in 0..n {
            for j in i..n {
                let oil_i = i >= v;
                let oil_j = j >= v;
                // Skip oil×oil terms (central OV property)
                if oil_i && oil_j {
                    continue;
                }
                a[idx_ut(n, i, j)] = F(rng.next_u32() as u8);
            }
        }
        f_quads.push(a);
    }

    // Sample random invertible transformation matrix T
    let t = sample_invertible(rng, n);

    // Compose P(x) = F(Tx) using direct coefficient computation
    // For the upper-triangular form Q(x) = sum_{i<=j} c_ij * x_i * x_j,
    // substituting y = Tx gives new coefficients:
    //   P_pp = sum_{i<=j} f_ij * T_ip * T_jp  (diagonal)
    //   P_pq = sum_{i<=j} f_ij * (T_ip * T_jq + T_iq * T_jp)  for p < q
    let mut polys = Vec::with_capacity(m);
    for k in 0..m {
        let f = &f_quads[k];
        let mut p = vec![F::ZERO; ut_len];

        for p_idx in 0..n {
            for q_idx in p_idx..n {
                let mut coeff = F::ZERO;
                for i in 0..n {
                    for j in i..n {
                        let f_ij = f[idx_ut(n, i, j)];
                        if f_ij.is_zero() {
                            continue;
                        }
                        if p_idx == q_idx {
                            // Diagonal term: f_ij * T_ip * T_jp
                            coeff += f_ij * t[i][p_idx] * t[j][p_idx];
                        } else {
                            // Off-diagonal term: f_ij * (T_ip * T_jq + T_iq * T_jp)
                            coeff += f_ij * (t[i][p_idx] * t[j][q_idx] + t[i][q_idx] * t[j][p_idx]);
                        }
                    }
                }
                p[idx_ut(n, p_idx, q_idx)] = coeff;
            }
        }

        polys.push(QuadForm::new(p));
    }

    (
        PublicKey { params, polys },
        SecretKey { params, t, f_quads },
    )
}

// Re-export PARAMS_L1 for backward compatibility
pub use crate::params::PARAMS_DEMO as PARAMS_L1;
