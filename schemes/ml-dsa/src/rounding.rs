//! Rounding and hint functions for ML-DSA.
//!
//! This module implements the rounding and hint functions defined in FIPS 204:
//! - Power2Round: Decomposes r into (r1, r0) where r = r1 * 2^d + r0
//! - Decompose: Decomposes r into (r1, r0) where r = r1 * α + r0
//! - HighBits: Extracts the high-order bits r1
//! - LowBits: Extracts the low-order bits r0
//! - MakeHint: Computes hint bit for recovering high bits
//! - UseHint: Recovers high bits using hint

use crate::params::{D, Q};
use crate::poly::Poly;
use crate::polyvec::PolyVec;
use crate::reduce::freeze;

/// Power2Round: Decomposes r into (r1, r0) where r = r1 * 2^d + r0.
///
/// r0 is in [-2^(d-1), 2^(d-1)].
/// r1 is in [0, (q-1)/2^d].
///
/// FIPS 204 Algorithm 30
pub fn power2round(r: i32) -> (i32, i32) {
    let r_pos = freeze(r);

    // r0 = r mod 2^d, centered in [-2^(d-1), 2^(d-1))
    let half_d = 1 << (D - 1); // 2^(d-1) = 4096
    let mask_d = (1 << D) - 1; // 2^d - 1 = 8191

    let r0_unsigned = r_pos & mask_d; // r mod 2^d
    let r0 = if r0_unsigned > half_d {
        r0_unsigned - (1 << D) // Center to negative
    } else {
        r0_unsigned
    };

    // r1 = (r - r0) / 2^d
    let r1 = (r_pos - r0) >> D;

    (r1, r0)
}

/// Applies Power2Round to all coefficients of a polynomial.
///
/// Returns (high, low) polynomials.
pub fn power2round_poly(r: &Poly) -> (Poly, Poly) {
    let mut high = Poly::zero();
    let mut low = Poly::zero();

    for i in 0..256 {
        let (r1, r0) = power2round(r.coeffs[i]);
        high.coeffs[i] = r1;
        low.coeffs[i] = r0;
    }

    (high, low)
}

/// Applies Power2Round to all polynomials in a vector.
///
/// Returns (t1, t0) vectors.
pub fn power2round_vec(t: &PolyVec) -> (PolyVec, PolyVec) {
    let len = t.len();
    let mut t1 = PolyVec::zero(len);
    let mut t0 = PolyVec::zero(len);

    for i in 0..len {
        let (high, low) = power2round_poly(&t.polys[i]);
        t1.polys[i] = high;
        t0.polys[i] = low;
    }

    (t1, t0)
}

/// Decompose: Decomposes r into (r1, r0) where r ≡ r1 * α + r0 (mod q).
///
/// For α = 2*γ2:
/// - r1 is the "high" part
/// - r0 is the "low" part in (-γ2, γ2]
///
/// FIPS 204 Algorithm 31
pub fn decompose(r: i32, gamma2: i32) -> (i32, i32) {
    let r_plus = freeze(r); // r in [0, Q)
    let alpha = 2 * gamma2;

    // r0 = r mod± α (centered modular reduction to (-γ2, γ2])
    // First compute r mod α in [0, α)
    let r_mod_alpha = r_plus % alpha;

    // Center to (-γ2, γ2]: if r_mod_alpha > γ2, subtract α
    let r0 = if r_mod_alpha > gamma2 {
        r_mod_alpha - alpha
    } else {
        r_mod_alpha
    };

    // Handle special case: r_plus - r0 = q - 1
    // This happens when r_plus is near q-1 and r0 would make r1 = (q-1)/α
    if r_plus - r0 == Q - 1 {
        return (0, r0 - 1);
    }

    // r1 = (r_plus - r0) / α
    let r1 = (r_plus - r0) / alpha;

    (r1, r0)
}

/// HighBits: Returns the high-order representative r1.
///
/// FIPS 204 Algorithm 32
pub fn high_bits(r: i32, gamma2: i32) -> i32 {
    let (r1, _) = decompose(r, gamma2);
    r1
}

/// LowBits: Returns the low-order representative r0.
///
/// FIPS 204 Algorithm 33
pub fn low_bits(r: i32, gamma2: i32) -> i32 {
    let (_, r0) = decompose(r, gamma2);
    r0
}

/// Applies HighBits to all coefficients of a polynomial.
pub fn high_bits_poly(r: &Poly, gamma2: i32) -> Poly {
    let mut result = Poly::zero();
    for i in 0..256 {
        result.coeffs[i] = high_bits(r.coeffs[i], gamma2);
    }
    result
}

/// Applies LowBits to all coefficients of a polynomial.
pub fn low_bits_poly(r: &Poly, gamma2: i32) -> Poly {
    let mut result = Poly::zero();
    for i in 0..256 {
        result.coeffs[i] = low_bits(r.coeffs[i], gamma2);
    }
    result
}

/// Applies HighBits to all polynomials in a vector.
pub fn high_bits_vec(r: &PolyVec, gamma2: i32) -> PolyVec {
    let mut result = PolyVec::zero(r.len());
    for i in 0..r.len() {
        result.polys[i] = high_bits_poly(&r.polys[i], gamma2);
    }
    result
}

/// Applies LowBits to all polynomials in a vector.
pub fn low_bits_vec(r: &PolyVec, gamma2: i32) -> PolyVec {
    let mut result = PolyVec::zero(r.len());
    for i in 0..r.len() {
        result.polys[i] = low_bits_poly(&r.polys[i], gamma2);
    }
    result
}

/// MakeHint: Computes hint bit h = 1 iff HighBits(r) ≠ HighBits(r + z).
///
/// The hint allows recovering HighBits(r) from HighBits(r + z).
///
/// FIPS 204 Algorithm 34
pub fn make_hint(z: i32, r: i32, gamma2: i32) -> bool {
    let r1 = high_bits(r, gamma2);
    let v1 = high_bits(freeze(r + z), gamma2);
    r1 != v1
}

/// UseHint: Recovers high bits using the hint.
///
/// If h = 0, returns HighBits(r).
/// If h = 1, returns the adjusted high bits.
///
/// FIPS 204 Algorithm 35
pub fn use_hint(h: bool, r: i32, gamma2: i32) -> i32 {
    let (r1, r0) = decompose(r, gamma2);

    if !h {
        return r1;
    }

    // Adjust r1 based on sign of r0
    let alpha = 2 * gamma2;
    let m = (Q - 1) / alpha; // Maximum value of r1

    // Handle special case zone: when r is in (Q-1-alpha/2, Q-1], decompose
    // returns r1=0 due to the special case, but the "true" r1 would be m.
    // We need to use the effective r1 for the adjustment.
    let r_plus = freeze(r);
    let in_special_case = r1 == 0 && r_plus > (Q - 1 - gamma2);
    let effective_r1 = if in_special_case { m } else { r1 };

    let result = if r0 > 0 {
        if effective_r1 == m {
            0
        } else {
            effective_r1 + 1
        }
    } else if effective_r1 == 0 {
        m
    } else {
        effective_r1 - 1
    };

    // HighBits never returns m (it returns 0 due to special case).
    // So if we would return m, return 0 instead.
    if result == m {
        0
    } else {
        result
    }
}

/// Computes hints for a polynomial vector.
///
/// Returns (hint_poly_vec, count) where count is the total number of 1s.
pub fn make_hint_vec(z: &PolyVec, r: &PolyVec, gamma2: i32) -> (Vec<Vec<bool>>, usize) {
    assert_eq!(z.len(), r.len());

    let mut hints = Vec::with_capacity(z.len());
    let mut count = 0;

    for i in 0..z.len() {
        let mut poly_hints = Vec::with_capacity(256);
        for j in 0..256 {
            let h = make_hint(z.polys[i].coeffs[j], r.polys[i].coeffs[j], gamma2);
            if h {
                count += 1;
            }
            poly_hints.push(h);
        }
        hints.push(poly_hints);
    }

    (hints, count)
}

/// Applies UseHint to a polynomial vector.
pub fn use_hint_vec(hints: &[Vec<bool>], r: &PolyVec, gamma2: i32) -> PolyVec {
    assert_eq!(hints.len(), r.len());

    let mut result = PolyVec::zero(r.len());

    for i in 0..r.len() {
        for j in 0..256 {
            result.polys[i].coeffs[j] = use_hint(hints[i][j], r.polys[i].coeffs[j], gamma2);
        }
    }

    result
}

/// Checks if the infinity norm of r0 is less than γ2 - β.
///
/// This is used in the signing loop to reject bad signatures.
pub fn check_low_bits(r: &PolyVec, gamma2: i32, beta: i32) -> bool {
    let bound = gamma2 - beta;

    for poly in &r.polys {
        for &c in &poly.coeffs {
            let r0 = low_bits(c, gamma2);
            if r0.abs() >= bound {
                return false;
            }
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_use_hint_correctness() {
        // Test that use_hint correctly recovers the original high bits
        let gamma2 = 95232; // ML-DSA-44
        let alpha = 2 * gamma2;

        for r1_orig in 0..44 {
            for r0_orig in [-gamma2/2, 0, gamma2/2].iter().copied() {
                let r_orig = r1_orig * alpha + r0_orig;
                if r_orig < 0 || r_orig >= Q { continue; }

                for ct0 in [-1000i32, 0, 1000, 50000, 95000].iter().copied() {
                    let r = freeze(r_orig + ct0);
                    let z = freeze(-ct0);

                    // What make_hint would compute
                    let r_plus_z = ((r as i64 + z as i64) % Q as i64) as i32;
                    let h = high_bits(r, gamma2) != high_bits(freeze(r_plus_z), gamma2);

                    // What use_hint should return
                    let recovered = use_hint(h, r, gamma2);

                    // Expected: high_bits(r + z) = high_bits(r_orig)
                    let expected = high_bits(r_orig, gamma2);

                    if recovered != expected {
                        println!("MISMATCH: r1_orig={}, r0_orig={}, ct0={}", r1_orig, r0_orig, ct0);
                        println!("  r={}, z={}, r+z={}", r, z, freeze(r_plus_z));
                        println!("  h={}, recovered={}, expected={}", h, recovered, expected);
                    }
                    assert_eq!(recovered, expected, "use_hint failed for r1_orig={}, r0_orig={}, ct0={}", r1_orig, r0_orig, ct0);
                }
            }
        }
    }

    #[test]
    fn test_power2round_basic() {
        // Test that r = r1 * 2^d + r0
        for r in [0, 1, 100, 1000, 10000, Q - 1, Q / 2] {
            let (r1, r0) = power2round(r);

            // r0 should be in [-2^(d-1), 2^(d-1)]
            let half_d = 1 << (D - 1);
            assert!(
                r0 >= -half_d && r0 <= half_d,
                "r0 = {} out of bounds for r = {}",
                r0,
                r
            );

            // Verify reconstruction
            let reconstructed = freeze(r1 * (1 << D) + r0);
            assert_eq!(
                reconstructed,
                freeze(r),
                "reconstruction failed for r = {}",
                r
            );
        }
    }

    #[test]
    fn test_power2round_d_value() {
        // D should be 13 for ML-DSA
        assert_eq!(D, 13);
    }

    #[test]
    fn test_decompose_reconstruction() {
        let gamma2 = 95232; // (Q - 1) / 88 for ML-DSA-44
        let alpha = 2 * gamma2;

        for r in [0, 1, 100, 1000, 10000, Q - 1, Q / 2, gamma2, gamma2 + 1] {
            let (r1, r0) = decompose(r, gamma2);

            // r0 should be in [-γ2, γ2]
            assert!(
                r0.abs() <= gamma2,
                "r0 = {} out of bounds for r = {}",
                r0,
                r
            );

            // Verify reconstruction: r ≡ r1 * α + r0 (mod q)
            let reconstructed = freeze(r1 * alpha + r0);
            assert_eq!(
                reconstructed,
                freeze(r),
                "reconstruction failed for r = {}",
                r
            );
        }
    }

    #[test]
    fn test_high_low_bits() {
        let gamma2 = 261888; // (Q - 1) / 32 for ML-DSA-65/87

        for r in [0, 1, 1000, gamma2, 2 * gamma2, Q - 1] {
            let h = high_bits(r, gamma2);
            let l = low_bits(r, gamma2);

            // Low bits should be in [-γ2, γ2]
            assert!(l.abs() <= gamma2);

            // High bits should be non-negative
            assert!(h >= 0);
        }
    }

    #[test]
    fn test_make_use_hint() {
        let gamma2 = 95232;

        // Test cases where hint is needed
        for r in [0, gamma2 - 10, gamma2 + 10, 2 * gamma2 - 10] {
            for z in [-10, 0, 10, 100] {
                let h = make_hint(z, r, gamma2);
                let r_plus_z = freeze(r + z);

                let original_high = high_bits(r_plus_z, gamma2);
                let recovered_high = use_hint(h, r_plus_z, gamma2);

                // UseHint should help recover correct high bits
                if !h {
                    assert_eq!(recovered_high, original_high);
                }
            }
        }
    }

    #[test]
    fn test_hint_preserves_high_bits() {
        let gamma2 = 261888;

        // When hint is false, high bits should be unchanged
        let r = 500000;
        let high = high_bits(r, gamma2);
        let recovered = use_hint(false, r, gamma2);
        assert_eq!(high, recovered);
    }

    #[test]
    fn test_power2round_poly() {
        let mut p = Poly::zero();
        p.coeffs[0] = 12345;
        p.coeffs[1] = Q - 1;
        p.coeffs[100] = Q / 2;

        let (high, low) = power2round_poly(&p);

        // Verify each coefficient
        for i in 0..256 {
            let (h, l) = power2round(p.coeffs[i]);
            assert_eq!(high.coeffs[i], h);
            assert_eq!(low.coeffs[i], l);
        }
    }

    #[test]
    fn test_check_low_bits() {
        let gamma2 = 95232;
        let beta = 78;

        let mut v = PolyVec::zero(4);
        // Set coefficients well within bounds
        for poly in &mut v.polys {
            for c in &mut poly.coeffs {
                *c = gamma2 / 2; // Well within gamma2 - beta
            }
        }

        assert!(check_low_bits(&v, gamma2, beta));

        // Set one coefficient close to bound
        v.polys[0].coeffs[0] = gamma2 - beta + 1; // Just over the bound
        assert!(!check_low_bits(&v, gamma2, beta));
    }
}
