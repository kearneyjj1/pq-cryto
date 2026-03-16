//! Optimized modular reduction algorithms for ML-DSA.
//!
//! This module provides Montgomery and Barrett reduction for efficient
//! modular arithmetic with q = 8380417.
//!
//! Montgomery reduction transforms values to a "Montgomery domain" where
//! multiplication becomes more efficient by avoiding explicit division.
//!
//! Barrett reduction uses a precomputed reciprocal approximation.

use crate::params::Q;

/// Montgomery parameter: R = 2^32
#[cfg(test)]
const MONTGOMERY_R: i64 = 1 << 32;

/// Montgomery parameter: R^(-1) mod Q
/// Computed as: pow(2^32, Q-2, Q)
pub const R_INV: i32 = 8265825;

/// Montgomery parameter: -Q^(-1) mod R for REDC algorithm
/// -Q^(-1) mod 2^32 = 4236238847
/// Satisfies: Q * Q_INV ≡ -1 (mod 2^32)
pub const Q_INV: i64 = 4236238847;

/// R mod Q = 2^32 mod Q = 4193792
pub const R_MOD_Q: i32 = 4193792;

/// R^2 mod Q = (2^32)^2 mod Q = 2365951
pub const R2_MOD_Q: i32 = 2365951;

/// Barrett reduction constant: floor(2^32 / Q)
/// = floor(4294967296 / 8380417) = 512
const BARRETT_MU: i64 = 512;

/// High-precision Barrett constant for i64 reduction: floor(2^58 / Q)
/// Used by reduce64 to handle products up to (Q-1)^2 without overflow.
const BARRETT_MU_HI: i64 = ((1i128 << 58) / (Q as i128)) as i64;

/// Simple modular reduction of a value in [-Q, 2Q) to [0, Q).
/// Constant-time: uses arithmetic masks instead of branches.
#[inline]
pub fn reduce32(a: i32) -> i32 {
    let mut t = a;
    // If t >= Q, subtract Q (mask is 1 if Q-1-t is negative)
    let mask_sub = ((Q - 1 - t) >> 31) & 1;
    t -= Q * mask_sub;
    // If t < 0, add Q (mask is 1 if t is negative)
    let mask_add = (t >> 31) & 1;
    t += Q * mask_add;
    t
}

/// Conditional reduction: if a >= Q, subtract Q.
/// Constant-time: uses arithmetic mask instead of branch.
#[inline]
pub fn cond_sub_q(a: i32) -> i32 {
    // mask is 1 if a >= Q (i.e., Q-1-a is negative), 0 otherwise
    let mask = ((Q - 1 - a) >> 31) & 1;
    a - Q * mask
}

/// Conditional addition: if a < 0, add Q.
/// Constant-time: uses arithmetic mask instead of branch.
#[inline]
pub fn cond_add_q(a: i32) -> i32 {
    // mask is 1 if a < 0, 0 otherwise
    let mask = (a >> 31) & 1;
    a + Q * mask
}

/// Montgomery reduction: compute a * R^(-1) mod Q.
///
/// Given a in [-Q*R, Q*R), returns t in (-Q, Q) such that t ≡ a * R^(-1) (mod Q).
///
/// The REDC algorithm:
/// 1. m = (a mod R) * (-Q^(-1) mod R) mod R
/// 2. t = (a + m * Q) / R
#[inline]
pub fn montgomery_reduce(a: i64) -> i32 {
    // m = (a * Q_INV) mod R
    // We only need the low 32 bits
    let m = ((a as i64).wrapping_mul(Q_INV as i64)) as i32;

    // t = (a + m * Q) / R
    // Since R = 2^32, division is a right shift
    let t = ((a + (m as i64) * (Q as i64)) >> 32) as i32;

    t
}

/// Converts a value to Montgomery domain: a -> a * R mod Q.
#[inline]
pub fn to_montgomery(a: i32) -> i32 {
    // a * R mod Q = a * R2 * R^(-1) mod Q = montgomery_reduce(a * R2)
    montgomery_reduce((a as i64) * (R2_MOD_Q as i64))
}

/// Converts a value from Montgomery domain: aR -> a mod Q.
#[inline]
pub fn from_montgomery(a: i32) -> i32 {
    montgomery_reduce(a as i64)
}

/// Montgomery multiplication: (aR * bR) -> abR mod Q.
#[inline]
pub fn montgomery_mul(a: i32, b: i32) -> i32 {
    montgomery_reduce((a as i64) * (b as i64))
}

/// Barrett reduction for values in [0, 2^32).
///
/// Uses the approximation: a mod Q ≈ a - floor(a * mu / 2^32) * Q
/// where mu = floor(2^32 / Q).
#[inline]
pub fn barrett_reduce(a: i32) -> i32 {
    // This works for a in [0, Q^2) approximately
    let t = ((a as i64 * BARRETT_MU) >> 32) as i32;
    let r = a - t * Q;

    // Final conditional subtraction
    cond_sub_q(r)
}

/// Reduces a 64-bit product to a value in [0, Q).
/// Constant-time: uses high-precision Barrett reduction instead of `%` operator.
/// Valid for inputs in [-(Q-1)^2, (Q-1)^2] (covers all products of i32 values in (-Q, Q)).
#[inline]
pub fn reduce64(a: i64) -> i32 {
    // High-precision Barrett: t = floor(a * mu / 2^58)
    // Uses i128 intermediate to avoid overflow for large products.
    let t = ((a as i128 * BARRETT_MU_HI as i128) >> 58) as i64;
    let mut r = (a - t * (Q as i64)) as i32;
    // r is in [-Q, 2Q) after Barrett approximation, correct with masks
    let mask_add = (r >> 31) & 1;
    r += Q * mask_add;
    let mask_sub = ((Q - 1 - r) >> 31) & 1;
    r -= Q * mask_sub;
    r
}

/// Freeze: reduce a to the canonical representative in [0, Q).
/// Constant-time: shifts input positive, then uses Barrett reduction.
#[inline]
pub fn freeze(a: i32) -> i32 {
    let a_shifted = ((a as i64) + (2 * Q as i64)) as i32;
    barrett_reduce(a_shifted)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce32() {
        assert_eq!(reduce32(0), 0);
        assert_eq!(reduce32(Q), 0);
        assert_eq!(reduce32(Q + 100), 100);
        assert_eq!(reduce32(-100), Q - 100);
    }

    #[test]
    fn test_cond_sub_q() {
        assert_eq!(cond_sub_q(Q), 0);
        assert_eq!(cond_sub_q(Q + 1), 1);
        assert_eq!(cond_sub_q(100), 100);
    }

    #[test]
    fn test_cond_add_q() {
        assert_eq!(cond_add_q(-1), Q - 1);
        assert_eq!(cond_add_q(-100), Q - 100);
        assert_eq!(cond_add_q(100), 100);
    }

    // Note: Montgomery functions are not currently used in the NTT implementation
    // (we use plain modular arithmetic instead). These tests are ignored until
    // Montgomery form is properly integrated.

    #[test]
    fn test_montgomery_roundtrip() {
        for val in [0, 1, 100, 12345, Q - 1] {
            let mont = to_montgomery(val);
            let back = freeze(from_montgomery(mont));
            assert_eq!(back, val, "roundtrip failed for {}", val);
        }
    }

    #[test]
    fn test_montgomery_mul() {
        let a = 12345;
        let b = 67890;
        let expected = ((a as i64 * b as i64) % (Q as i64)) as i32;

        let a_mont = to_montgomery(a);
        let b_mont = to_montgomery(b);
        let result_mont = montgomery_mul(a_mont, b_mont);
        let result = freeze(from_montgomery(result_mont));

        assert_eq!(result, expected);
    }

    #[test]
    fn test_barrett_reduce() {
        assert_eq!(barrett_reduce(0), 0);
        assert_eq!(barrett_reduce(Q), 0);
        assert_eq!(barrett_reduce(100), 100);
    }

    #[test]
    fn test_reduce64() {
        assert_eq!(reduce64(0), 0);
        assert_eq!(reduce64(Q as i64), 0);
        assert_eq!(reduce64(Q as i64 + 100), 100);
        assert_eq!(reduce64(-100), Q - 100);
    }

    #[test]
    fn test_freeze() {
        assert_eq!(freeze(0), 0);
        assert_eq!(freeze(Q), 0);
        assert_eq!(freeze(2 * Q + 100), 100);
        assert_eq!(freeze(-100), Q - 100);
        assert_eq!(freeze(-Q - 100), Q - 100);
    }

    #[test]
    fn test_constants() {
        // Verify R_MOD_Q = 2^32 mod Q
        assert_eq!((MONTGOMERY_R % (Q as i64)) as i32, R_MOD_Q);

        // Verify R2_MOD_Q = 2^64 mod Q
        let r2 = ((MONTGOMERY_R as u128 * MONTGOMERY_R as u128) % (Q as u128)) as i32;
        assert_eq!(r2, R2_MOD_Q);

        // Verify R_INV: R * R_INV ≡ 1 (mod Q)
        let check = (((MONTGOMERY_R % (Q as i64)) * (R_INV as i64)) % (Q as i64)) as i32;
        assert_eq!(check, 1);
    }

    #[test]
    fn test_montgomery_reduce_range() {
        // Montgomery reduce should handle the full range used in NTT
        let a: i64 = (Q as i64 - 1) * (Q as i64 - 1);
        let result = montgomery_reduce(a);
        // Result should be representable and correct when frozen
        let frozen = freeze(result);
        let expected = ((a % (Q as i64)) as i32 * R_INV) % Q;
        let expected = if expected < 0 { expected + Q } else { expected };
        assert_eq!(frozen, expected);
    }

    #[test]
    fn test_q_inv_correct() {
        // Q_INV should satisfy: Q * Q_INV ≡ -1 (mod 2^32)
        // Using Newton's method (Hensel lifting) to compute -Q^(-1) mod 2^32
        let q = Q as u64;
        // Start: x = Q (Q * Q = Q^2 ≡ ? mod 2, Q is odd so Q^2 is odd, Q * 1 = Q ≡ 1 mod 2)
        // We need x such that Q*x ≡ -1 mod 2^32
        // First find Q^(-1) mod 2^32 using Newton: x_{n+1} = x_n * (2 - Q * x_n)
        let mut x: u64 = 1; // Q * 1 ≡ 1 mod 2 (Q is odd)
        for _ in 0..5 { // 5 iterations: 1 -> 2 -> 4 -> 8 -> 16 -> 32 bits
            x = x.wrapping_mul(2u64.wrapping_sub(q.wrapping_mul(x)));
        }
        let x = x & 0xFFFFFFFF; // mod 2^32
        // Verify: Q * x ≡ 1 (mod 2^32)
        assert_eq!((q.wrapping_mul(x)) & 0xFFFFFFFF, 1, "Q^(-1) computation failed");
        // -Q^(-1) mod 2^32
        let neg_q_inv = (0u64.wrapping_sub(x)) & 0xFFFFFFFF;
        // Verify: Q * neg_q_inv ≡ -1 (mod 2^32)
        assert_eq!((q.wrapping_mul(neg_q_inv)) & 0xFFFFFFFF, 0xFFFFFFFF,
            "Neg Q^(-1) verification failed");
        eprintln!("Correct Q_INV (=-Q^{{-1}} mod 2^32) = {}", neg_q_inv);
        eprintln!("Current Q_INV = {}", Q_INV);
        assert_eq!(Q_INV as u64, neg_q_inv, "Q_INV constant is incorrect");
    }

    #[test]
    fn test_freeze_ct_exhaustive() {
        // Test freeze for a representative range
        for a in (-2 * Q)..=(2 * Q) {
            let expected = {
                let mut t = a % Q;
                if t < 0 { t += Q; }
                t
            };
            assert_eq!(freeze(a), expected, "freeze({}) failed", a);
        }
    }

    #[test]
    fn test_reduce64_ct() {
        // Test reduce64 for products of values in [0, Q)
        for a in [0i64, 1, 100, Q as i64 - 1, Q as i64, -1, -100, -(Q as i64) + 1] {
            for b in [0i64, 1, 100, Q as i64 - 1] {
                let prod = a * b;
                let expected = {
                    let mut r = (prod % Q as i64) as i32;
                    if r < 0 { r += Q; }
                    r
                };
                assert_eq!(reduce64(prod), expected, "reduce64({}) failed", prod);
            }
        }
    }
}
