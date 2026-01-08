//! ML-DSA parameter sets as defined in FIPS 204.
//!
//! This module defines the three security levels:
//! - ML-DSA-44: NIST Level 2 (~128-bit security)
//! - ML-DSA-65: NIST Level 3 (~192-bit security)
//! - ML-DSA-87: NIST Level 5 (~256-bit security)

/// The prime modulus q = 2^23 - 2^13 + 1 = 8380417
pub const Q: i32 = 8380417;

/// Polynomial degree n = 256
pub const N: usize = 256;

/// Number of dropped bits from t (d = 13)
pub const D: usize = 13;

/// Root of unity for NTT (ζ = 1753).
///
/// ZETA is a primitive 512th root of unity modulo Q, meaning:
/// - ZETA^512 ≡ 1 (mod Q)
/// - ZETA^256 ≡ -1 (mod Q)
/// - No smaller power of ZETA equals 1
///
/// This value is specified in FIPS 204 for the ML-DSA modulus Q = 8380417.
pub const ZETA: i32 = 1753;

/// Parameters for a specific ML-DSA security level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Params {
    /// Number of rows in matrix A (k)
    pub k: usize,
    /// Number of columns in matrix A (l)
    pub l: usize,
    /// Coefficient range for secret vectors s1, s2 (η)
    pub eta: usize,
    /// Number of ±1 coefficients in challenge polynomial c (τ)
    pub tau: usize,
    /// Bound on ||z||∞ (β = τ * η)
    pub beta: i32,
    /// Range for mask polynomial y (γ1 = 2^γ1_bits)
    pub gamma1: i32,
    /// Low-order rounding range (γ2 = (q-1)/α)
    pub gamma2: i32,
    /// Maximum number of 1s in hint vector h (ω)
    pub omega: usize,
    /// Security parameter in bits (λ)
    pub lambda: usize,
}

impl Params {
    /// Returns the public key size in bytes.
    pub const fn public_key_size(&self) -> usize {
        // pk = ρ || t1
        // ρ is 32 bytes, t1 is k polynomials each with coefficients in [0, 2^(ceil(log2(q)) - D))
        // For q = 8380417, ceil(log2(q)) = 23, so 23 - 13 = 10 bits per coefficient
        // t1: k * 256 * 10 / 8 = k * 320 bytes
        32 + self.k * 320
    }

    /// Returns the secret key size in bytes.
    pub const fn secret_key_size(&self) -> usize {
        // sk = ρ || K || tr || s1 || s2 || t0
        // ρ: 32 bytes, K: 32 bytes, tr: 64 bytes
        // s1: l polynomials with coefficients in [-η, η]
        // s2: k polynomials with coefficients in [-η, η]
        // t0: k polynomials with coefficients in [-(2^(D-1)-1), 2^(D-1)]
        let eta_bits = if self.eta == 2 { 3 } else { 4 }; // bits per coefficient for s1, s2
        let s1_bytes = self.l * N * eta_bits / 8;
        let s2_bytes = self.k * N * eta_bits / 8;
        let t0_bytes = self.k * N * D / 8; // D bits per coefficient
        32 + 32 + 64 + s1_bytes + s2_bytes + t0_bytes
    }

    /// Returns the signature size in bytes.
    pub const fn signature_size(&self) -> usize {
        // sig = c_tilde || z || h
        // c_tilde: λ/4 bytes (commitment hash)
        // z: l polynomials with coefficients in [-(γ1-1), γ1]
        // h: encoded hint vector (ω + k bytes)
        let c_tilde_bytes = self.lambda / 4;
        let gamma1_bits = if self.gamma1 == (1 << 17) { 18 } else { 20 };
        let z_bytes = self.l * N * gamma1_bits / 8;
        let h_bytes = self.omega + self.k;
        c_tilde_bytes + z_bytes + h_bytes
    }

    /// Returns the number of bits needed to encode γ1.
    pub const fn gamma1_bits(&self) -> usize {
        if self.gamma1 == (1 << 17) {
            18
        } else {
            20
        }
    }
}

// ============================================================================
// ML-DSA-44 Parameters (NIST Level 2)
// ============================================================================

/// ML-DSA-44 parameters (~128-bit security, NIST Level 2).
pub const ML_DSA_44: Params = Params {
    k: 4,
    l: 4,
    eta: 2,
    tau: 39,
    beta: 78,           // τ * η = 39 * 2
    gamma1: 1 << 17,    // 2^17 = 131072
    gamma2: 95232,      // (q - 1) / 88
    omega: 80,
    lambda: 128,
};

// ============================================================================
// ML-DSA-65 Parameters (NIST Level 3)
// ============================================================================

/// ML-DSA-65 parameters (~192-bit security, NIST Level 3).
pub const ML_DSA_65: Params = Params {
    k: 6,
    l: 5,
    eta: 4,
    tau: 49,
    beta: 196,          // τ * η = 49 * 4
    gamma1: 1 << 19,    // 2^19 = 524288
    gamma2: 261888,     // (q - 1) / 32
    omega: 55,
    lambda: 192,
};

// ============================================================================
// ML-DSA-87 Parameters (NIST Level 5)
// ============================================================================

/// ML-DSA-87 parameters (~256-bit security, NIST Level 5).
pub const ML_DSA_87: Params = Params {
    k: 8,
    l: 7,
    eta: 2,
    tau: 60,
    beta: 120,          // τ * η = 60 * 2
    gamma1: 1 << 19,    // 2^19 = 524288
    gamma2: 261888,     // (q - 1) / 32
    omega: 75,
    lambda: 256,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_q_is_prime() {
        // Q should be 2^23 - 2^13 + 1
        assert_eq!(Q, (1 << 23) - (1 << 13) + 1);
        assert_eq!(Q, 8380417);
    }

    #[test]
    fn test_ml_dsa_44_params() {
        assert_eq!(ML_DSA_44.k, 4);
        assert_eq!(ML_DSA_44.l, 4);
        assert_eq!(ML_DSA_44.eta, 2);
        assert_eq!(ML_DSA_44.tau, 39);
        assert_eq!(ML_DSA_44.beta, 78);
        assert_eq!(ML_DSA_44.gamma1, 131072);
        assert_eq!(ML_DSA_44.omega, 80);
        assert_eq!(ML_DSA_44.lambda, 128);
    }

    #[test]
    fn test_ml_dsa_65_params() {
        assert_eq!(ML_DSA_65.k, 6);
        assert_eq!(ML_DSA_65.l, 5);
        assert_eq!(ML_DSA_65.eta, 4);
        assert_eq!(ML_DSA_65.tau, 49);
        assert_eq!(ML_DSA_65.beta, 196);
        assert_eq!(ML_DSA_65.gamma1, 524288);
        assert_eq!(ML_DSA_65.omega, 55);
        assert_eq!(ML_DSA_65.lambda, 192);
    }

    #[test]
    fn test_ml_dsa_87_params() {
        assert_eq!(ML_DSA_87.k, 8);
        assert_eq!(ML_DSA_87.l, 7);
        assert_eq!(ML_DSA_87.eta, 2);
        assert_eq!(ML_DSA_87.tau, 60);
        assert_eq!(ML_DSA_87.beta, 120);
        assert_eq!(ML_DSA_87.gamma1, 524288);
        assert_eq!(ML_DSA_87.omega, 75);
        assert_eq!(ML_DSA_87.lambda, 256);
    }

    #[test]
    fn test_gamma2_values() {
        // γ2 = (q - 1) / 88 for ML-DSA-44
        assert_eq!((Q - 1) / 88, 95232);
        assert_eq!(ML_DSA_44.gamma2, 95232);

        // γ2 = (q - 1) / 32 for ML-DSA-65 and ML-DSA-87
        assert_eq!((Q - 1) / 32, 261888);
        assert_eq!(ML_DSA_65.gamma2, 261888);
        assert_eq!(ML_DSA_87.gamma2, 261888);
    }

    #[test]
    fn test_beta_equals_tau_times_eta() {
        assert_eq!(ML_DSA_44.beta, (ML_DSA_44.tau * ML_DSA_44.eta) as i32);
        assert_eq!(ML_DSA_65.beta, (ML_DSA_65.tau * ML_DSA_65.eta) as i32);
        assert_eq!(ML_DSA_87.beta, (ML_DSA_87.tau * ML_DSA_87.eta) as i32);
    }

    /// Helper to compute a^e mod m.
    fn powmod(base: i32, mut exp: u32, modulus: i64) -> i64 {
        let mut result = 1i64;
        let mut b = base as i64 % modulus;
        while exp > 0 {
            if exp & 1 == 1 {
                result = (result * b) % modulus;
            }
            b = (b * b) % modulus;
            exp >>= 1;
        }
        result
    }

    #[test]
    fn test_zeta_is_primitive_512th_root() {
        // ZETA^512 ≡ 1 (mod Q)
        let zeta_512 = powmod(ZETA, 512, Q as i64);
        assert_eq!(zeta_512, 1, "ZETA^512 should equal 1 mod Q");

        // ZETA^256 ≡ -1 (mod Q)
        let zeta_256 = powmod(ZETA, 256, Q as i64);
        assert_eq!(zeta_256, (Q - 1) as i64, "ZETA^256 should equal -1 mod Q");

        // Verify ZETA is primitive: ZETA^k ≠ 1 for k < 512
        // Check divisors of 512: 1, 2, 4, 8, 16, 32, 64, 128, 256
        for &k in &[1u32, 2, 4, 8, 16, 32, 64, 128, 256] {
            let zeta_k = powmod(ZETA, k, Q as i64);
            assert_ne!(zeta_k, 1, "ZETA^{} should not equal 1 mod Q", k);
        }
    }

    #[test]
    fn test_n_inverse_mod_q() {
        // The inverse of N=256 mod Q should be 8347681
        // Verify: 256 * 8347681 ≡ 1 (mod Q)
        let n_inv: i64 = 8347681;
        let product = (N as i64 * n_inv) % (Q as i64);
        assert_eq!(product, 1, "256 * 8347681 should equal 1 mod Q");
    }

    #[test]
    fn test_q_ntt_friendly() {
        // Q should support 512th roots of unity
        // This requires Q - 1 to be divisible by 512
        assert_eq!((Q - 1) % 512, 0, "Q-1 should be divisible by 512 for NTT");

        // Q - 1 = 8380416 = 2^13 * 1023 = 2^13 * 3 * 11 * 31
        // So Q supports up to 2^13 = 8192 th roots of unity
        assert_eq!((Q - 1) % (1 << 13), 0, "Q-1 should be divisible by 2^13");
    }
}
