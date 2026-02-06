//! Parameter sets for FN-DSA (FALCON).
//!
//! This module defines the two FALCON parameter sets:
//! - FALCON-512: NIST Level 1 (~128-bit security)
//! - FALCON-1024: NIST Level 5 (~256-bit security)
//!
//! Both use the same modulus q = 12289 and polynomial ring Z[X]/(X^n + 1).

/// The FALCON modulus q = 12289.
///
/// This is an NTT-friendly prime: q = 12*1024 + 1 = 3*2^12 + 1.
/// The 2n-th roots of unity exist in Z_q for n up to 1024.
pub const Q: i32 = 12289;

/// Log2 of the modulus (for bit operations).
pub const Q_BITS: usize = 14;

/// Parameters for the FN-DSA (FALCON) signature scheme.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Params {
    /// Polynomial degree n (512 or 1024).
    pub n: usize,

    /// Log2 of n (9 for 512, 10 for 1024).
    pub log_n: usize,

    /// Gaussian standard deviation sigma for signing.
    /// This is sigma = 1.17 * sqrt(q) * sqrt(2n / (2n-1)) approximately.
    pub sigma: f64,

    /// Minimum standard deviation at LDL* tree leaves (sigma_min).
    /// Different for FALCON-512 vs FALCON-1024 per the specification.
    pub sigma_min: f64,

    /// Signature bound squared (for norm checking).
    pub sig_bound_sq: f64,

    /// Public key size in bytes.
    pub pk_bytes: usize,

    /// Secret key size in bytes.
    pub sk_bytes: usize,

    /// Maximum signature size in bytes.
    pub sig_bytes_max: usize,

    /// Golomb-Rice parameter k for signature compression.
    /// Low k bits of each coefficient's absolute value are coded in binary;
    /// the remaining high bits are coded in unary. k = 8 for FALCON-512,
    /// k = 9 for FALCON-1024.
    pub rice_k: usize,

    /// NIST security level (1 or 5).
    pub security_level: usize,
}

impl Params {
    /// Returns the polynomial degree n.
    #[inline]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// Returns the number of coefficients in the upper triangular representation.
    #[inline]
    pub const fn ut_size(&self) -> usize {
        self.n * (self.n + 1) / 2
    }

    /// Computes the signature norm bound for verification.
    #[inline]
    pub fn sig_bound(&self) -> f64 {
        self.sig_bound_sq.sqrt()
    }
}

/// FALCON-512 parameters (NIST Level 1, ~128-bit security).
///
/// - Polynomial degree: n = 512
/// - Public key: 897 bytes
/// - Secret key: 1281 bytes
/// - Signature: ~666 bytes (variable, max ~809)
///
/// The Gaussian standard deviation sigma ≈ 165.74 is chosen such that
/// signatures satisfy the norm bound with high probability.
pub const FALCON_512: Params = Params {
    n: 512,
    log_n: 9,
    // sigma = 1.17 * sqrt(q) * sqrt(2*512 / (2*512 - 1)) per reference impl
    sigma: 165.7366171829776,
    sigma_min: 1.2778336969128337,
    // sig_bound^2 per FIPS 206: floor(beta^2 * 2n * sigma^2)
    // = floor(1.1^2 * 2 * 512 * 165.7366...^2) ≈ 34,034,726
    sig_bound_sq: 34_034_726.0,
    pk_bytes: 897,
    sk_bytes: 1281,
    sig_bytes_max: 809, // Worst case, typical is ~666
    rice_k: 8,
    security_level: 1,
};

/// FALCON-1024 parameters (NIST Level 5, ~256-bit security).
///
/// - Polynomial degree: n = 1024
/// - Public key: 1793 bytes
/// - Secret key: 2305 bytes
/// - Signature: ~1280 bytes (variable, max ~1577)
///
/// The Gaussian standard deviation sigma ≈ 168.39 is slightly larger
/// than FALCON-512 to maintain the norm bound.
pub const FALCON_1024: Params = Params {
    n: 1024,
    log_n: 10,
    // sigma = 1.17 * sqrt(q) * sqrt(2*1024 / (2*1024 - 1)) per reference impl
    sigma: 168.38857144654395,
    sigma_min: 1.298280334344292,
    // sig_bound^2 = (1.1 * sigma * sqrt(2n))^2
    sig_bound_sq: 70265242.0,
    pk_bytes: 1793,
    sk_bytes: 2305,
    sig_bytes_max: 1577, // Worst case, typical is ~1280
    rice_k: 9,
    security_level: 5,
};

/// FALCON-16 parameters (for testing only, NOT SECURE).
///
/// - Polynomial degree: n = 16
/// - This is a toy parameter set for unit testing.
/// - DO NOT USE IN PRODUCTION.
///
/// NOTE: NTRUSolve for small n produces F, G with large coefficients,
/// so signatures are much larger than standard FALCON. We use a very
/// relaxed bound that allows valid signatures while still rejecting
/// wrong-message verifications.
#[cfg(test)]
pub const FALCON_16: Params = Params {
    n: 16,
    log_n: 4,
    sigma: 165.7, // Use FALCON-512 sigma for consistency
    sigma_min: 1.2778336969128337, // Use FALCON-512 value
    sig_bound_sq: 6e8, // 600M — relaxed for n=16 (NTRUSolve produces large F,G)
    pk_bytes: 64,
    sk_bytes: 128,
    sig_bytes_max: 64,
    rice_k: 8,
    security_level: 0, // Not secure
};

/// Beta parameter for signature compression.
/// Used in the Golomb-Rice style encoding.
pub const BETA: f64 = 0.5;

/// Maximum number of signing attempts before giving up.
/// With proper ffSampling, signing typically succeeds within a few attempts.
pub const MAX_SIGN_ATTEMPTS: u32 = 100;

/// Nonce size in bytes for signing.
pub const NONCE_SIZE: usize = 40;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_falcon_512_params() {
        assert_eq!(FALCON_512.n, 512);
        assert_eq!(FALCON_512.log_n, 9);
        assert_eq!(FALCON_512.pk_bytes, 897);
        assert_eq!(FALCON_512.sk_bytes, 1281);
        assert_eq!(FALCON_512.security_level, 1);
        assert!(FALCON_512.sigma > 165.0 && FALCON_512.sigma < 166.0);
    }

    #[test]
    fn test_falcon_1024_params() {
        assert_eq!(FALCON_1024.n, 1024);
        assert_eq!(FALCON_1024.log_n, 10);
        assert_eq!(FALCON_1024.pk_bytes, 1793);
        assert_eq!(FALCON_1024.sk_bytes, 2305);
        assert_eq!(FALCON_1024.security_level, 5);
        assert!(FALCON_1024.sigma > 168.0 && FALCON_1024.sigma < 169.0);
    }

    #[test]
    fn test_modulus() {
        // q = 12289 = 3 * 4096 + 1 = 3 * 2^12 + 1
        assert_eq!(Q, 12289);
        assert_eq!(Q, 3 * 4096 + 1);
    }

    #[test]
    fn test_q_is_prime() {
        // Simple primality check for small q
        let q = Q as u64;
        for i in 2..((q as f64).sqrt() as u64 + 1) {
            assert_ne!(q % i, 0, "q should be prime");
        }
    }

    #[test]
    fn test_n_is_power_of_two() {
        assert!(FALCON_512.n.is_power_of_two());
        assert!(FALCON_1024.n.is_power_of_two());
        assert_eq!(1 << FALCON_512.log_n, FALCON_512.n);
        assert_eq!(1 << FALCON_1024.log_n, FALCON_1024.n);
    }

    #[test]
    fn test_sig_bound() {
        // Check that sig_bound is reasonable (should be around sqrt(sig_bound_sq))
        let bound_512 = FALCON_512.sig_bound();
        let bound_1024 = FALCON_1024.sig_bound();

        // FALCON-512: sqrt(34034726) ≈ 5834
        assert!(bound_512 > 5800.0 && bound_512 < 5900.0);

        // FALCON-1024: sqrt(70265242) ≈ 8382
        assert!(bound_1024 > 8300.0 && bound_1024 < 8400.0);
    }
}
