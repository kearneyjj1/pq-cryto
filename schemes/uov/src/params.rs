//! UOV parameter sets including NIST PQC standardization candidates.
//!
//! This module defines parameter sets for the UOV (Unbalanced Oil and Vinegar)
//! signature scheme. Parameters are defined by:
//! - `v`: Number of vinegar variables
//! - `m`: Number of oil variables (also the number of equations)
//! - `n = v + m`: Total number of variables
//!
//! The security level roughly corresponds to NIST security levels:
//! - Level 1: ~128-bit security (comparable to AES-128)
//! - Level 3: ~192-bit security (comparable to AES-192)
//! - Level 5: ~256-bit security (comparable to AES-256)

use crate::error::{Result, UovError};

/// Parameters for the UOV signature scheme.
///
/// The UOV scheme uses `v` vinegar variables and `m` oil variables,
/// with `n = v + m` total variables. The scheme produces `m` quadratic
/// equations in `n` variables.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Params {
    /// Number of vinegar variables.
    pub v: usize,
    /// Number of oil variables (also number of equations).
    pub m: usize,
}

impl Params {
    /// Creates new parameters with validation.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `v` is zero
    /// - `m` is zero
    /// - `v < m` (UOV requires more vinegar than oil variables)
    pub fn new(v: usize, m: usize) -> Result<Self> {
        if v == 0 {
            return Err(UovError::InvalidParams {
                reason: "vinegar count v must be positive",
            });
        }
        if m == 0 {
            return Err(UovError::InvalidParams {
                reason: "oil count m must be positive",
            });
        }
        if v < m {
            return Err(UovError::InvalidParams {
                reason: "UOV requires v >= m for security",
            });
        }
        Ok(Params { v, m })
    }

    /// Creates parameters without validation (for use with known-good constants).
    #[inline]
    pub const fn new_unchecked(v: usize, m: usize) -> Self {
        Params { v, m }
    }

    /// Returns the total number of variables (n = v + m).
    #[inline]
    pub const fn n(&self) -> usize {
        self.v + self.m
    }

    /// Returns the size of an upper triangular matrix representation.
    ///
    /// For n variables, this is n*(n+1)/2 coefficients.
    #[inline]
    pub const fn ut_size(&self) -> usize {
        let n = self.n();
        n * (n + 1) / 2
    }

    /// Returns the public key size in field elements.
    ///
    /// The public key contains m quadratic forms, each with ut_size() coefficients.
    #[inline]
    pub const fn public_key_size(&self) -> usize {
        self.m * self.ut_size()
    }

    /// Returns the signature size in field elements (excluding salt).
    #[inline]
    pub const fn signature_size(&self) -> usize {
        self.n()
    }

    /// Returns the salt size in bytes.
    #[inline]
    pub const fn salt_size(&self) -> usize {
        self.m
    }
}

// ============================================================================
// Demo/Test Parameters (small, for fast testing)
// ============================================================================

/// Demo parameters for testing (NOT for production use).
///
/// These parameters use v=16, m=8 for fast key generation and signing
/// during development and testing. They do NOT provide adequate security.
pub const PARAMS_DEMO: Params = Params::new_unchecked(16, 8);

/// Alias for demo parameters (backward compatibility).
pub const PARAMS_L1_DEMO: Params = PARAMS_DEMO;

// ============================================================================
// NIST PQC Round 2 Parameters
// ============================================================================

/// NIST Level 1 parameters (~128-bit security).
///
/// Uses v=68, m=44, n=112. Provides security roughly equivalent to AES-128.
///
/// Key sizes:
/// - Public key: ~50 KB
/// - Secret key: ~12 KB
/// - Signature: 112 bytes + 44 byte salt
pub const PARAMS_NIST_L1: Params = Params::new_unchecked(68, 44);

/// NIST Level 3 parameters (~192-bit security).
///
/// Uses v=148, m=96, n=244. Provides security roughly equivalent to AES-192.
///
/// Key sizes:
/// - Public key: ~240 KB
/// - Secret key: ~55 KB
/// - Signature: 244 bytes + 96 byte salt
pub const PARAMS_NIST_L3: Params = Params::new_unchecked(148, 96);

/// NIST Level 5 parameters (~256-bit security).
///
/// Uses v=244, m=160, n=404. Provides security roughly equivalent to AES-256.
///
/// Key sizes:
/// - Public key: ~660 KB
/// - Secret key: ~150 KB
/// - Signature: 404 bytes + 160 byte salt
pub const PARAMS_NIST_L5: Params = Params::new_unchecked(244, 160);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_params_new_valid() {
        let p = Params::new(16, 8).unwrap();
        assert_eq!(p.v, 16);
        assert_eq!(p.m, 8);
        assert_eq!(p.n(), 24);
    }

    #[test]
    fn test_params_new_invalid_v_zero() {
        assert!(Params::new(0, 8).is_err());
    }

    #[test]
    fn test_params_new_invalid_m_zero() {
        assert!(Params::new(16, 0).is_err());
    }

    #[test]
    fn test_params_new_invalid_v_less_than_m() {
        assert!(Params::new(4, 8).is_err());
    }

    #[test]
    fn test_ut_size() {
        let p = Params::new_unchecked(16, 8);
        // n = 24, ut_size = 24 * 25 / 2 = 300
        assert_eq!(p.ut_size(), 300);
    }

    #[test]
    fn test_nist_l1_params() {
        assert_eq!(PARAMS_NIST_L1.v, 68);
        assert_eq!(PARAMS_NIST_L1.m, 44);
        assert_eq!(PARAMS_NIST_L1.n(), 112);
    }

    #[test]
    fn test_nist_l3_params() {
        assert_eq!(PARAMS_NIST_L3.v, 148);
        assert_eq!(PARAMS_NIST_L3.m, 96);
        assert_eq!(PARAMS_NIST_L3.n(), 244);
    }

    #[test]
    fn test_nist_l5_params() {
        assert_eq!(PARAMS_NIST_L5.v, 244);
        assert_eq!(PARAMS_NIST_L5.m, 160);
        assert_eq!(PARAMS_NIST_L5.n(), 404);
    }
}
