//! Polynomial operations for ML-DSA.
//!
//! Polynomials are elements of the ring R_q = Z_q[X]/(X^256 + 1).
//! They are represented as arrays of 256 coefficients in Z_q.
//!
//! Operations can be performed in either:
//! - Coefficient domain: standard polynomial representation
//! - NTT domain: transformed representation for efficient multiplication

use crate::ntt::{inv_ntt, ntt, pointwise_mul, reduce_coeffs};
use crate::params::{N, Q};
use crate::reduce::{cond_add_q, cond_sub_q, freeze, montgomery_reduce};
use std::ops::{Add, AddAssign, Sub, SubAssign};

/// A polynomial in R_q = Z_q[X]/(X^256 + 1).
///
/// The polynomial is stored as an array of coefficients where
/// coeffs[i] is the coefficient of X^i.
#[derive(Clone, Debug)]
pub struct Poly {
    /// Coefficients of the polynomial.
    pub coeffs: [i32; N],
}

impl Default for Poly {
    fn default() -> Self {
        Self::zero()
    }
}

impl Poly {
    /// Creates the zero polynomial.
    pub fn zero() -> Self {
        Poly { coeffs: [0; N] }
    }

    /// Creates a polynomial from a coefficient array.
    pub fn from_coeffs(coeffs: [i32; N]) -> Self {
        Poly { coeffs }
    }

    /// Returns true if all coefficients are zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c == 0)
    }

    /// Reduces all coefficients to [0, Q).
    pub fn reduce(&mut self) {
        reduce_coeffs(&mut self.coeffs);
    }

    /// Returns a reduced copy of this polynomial.
    pub fn reduced(&self) -> Self {
        let mut p = self.clone();
        p.reduce();
        p
    }

    /// Applies forward NTT in-place.
    ///
    /// After this operation, the polynomial is in NTT domain.
    pub fn ntt(&mut self) {
        // Note: We use plain modular arithmetic NTT (not Montgomery form)
        ntt(&mut self.coeffs);
    }

    /// Applies inverse NTT in-place.
    ///
    /// After this operation, the polynomial is in coefficient domain.
    pub fn inv_ntt(&mut self) {
        inv_ntt(&mut self.coeffs);
        self.reduce();
    }

    /// Pointwise multiplication in NTT domain.
    ///
    /// Both self and other must be in NTT domain.
    /// Result is in NTT domain.
    pub fn pointwise_mul(&self, other: &Poly) -> Poly {
        let mut result = Poly::zero();
        pointwise_mul(&mut result.coeffs, &self.coeffs, &other.coeffs);
        result
    }

    /// Computes self + c * other where c is a scalar.
    pub fn add_scaled(&mut self, other: &Poly, c: i32) {
        for i in 0..N {
            self.coeffs[i] += c * other.coeffs[i];
        }
    }

    /// Computes self - c * other where c is a scalar.
    pub fn sub_scaled(&mut self, other: &Poly, c: i32) {
        for i in 0..N {
            self.coeffs[i] -= c * other.coeffs[i];
        }
    }

    /// Multiplies all coefficients by a scalar.
    pub fn scale(&mut self, c: i32) {
        for coeff in self.coeffs.iter_mut() {
            *coeff = montgomery_reduce((*coeff as i64) * (c as i64));
        }
    }

    /// Shifts coefficients left by d bits (multiply by 2^d).
    pub fn shift_left(&mut self, d: usize) {
        for coeff in self.coeffs.iter_mut() {
            *coeff <<= d;
        }
    }

    /// Computes the infinity norm ||f||_∞.
    ///
    /// Returns the maximum absolute value of any coefficient
    /// when represented in centered form [-Q/2, Q/2].
    pub fn norm_inf(&self) -> i32 {
        let mut max = 0i32;
        let half_q = (Q - 1) / 2;

        for &coeff in self.coeffs.iter() {
            let c = freeze(coeff);
            let centered = if c > half_q { Q - c } else { c };
            if centered > max {
                max = centered;
            }
        }

        max
    }

    /// Checks if ||f||_∞ < bound.
    pub fn check_norm(&self, bound: i32) -> bool {
        self.norm_inf() < bound
    }

    /// Applies conditional subtraction to all coefficients.
    pub fn cond_sub_q(&mut self) {
        for coeff in self.coeffs.iter_mut() {
            *coeff = cond_sub_q(*coeff);
        }
    }

    /// Applies conditional addition to all coefficients.
    pub fn cond_add_q(&mut self) {
        for coeff in self.coeffs.iter_mut() {
            *coeff = cond_add_q(*coeff);
        }
    }

    /// Freezes all coefficients to [0, Q).
    pub fn freeze(&mut self) {
        for coeff in self.coeffs.iter_mut() {
            *coeff = freeze(*coeff);
        }
    }
}

impl Add for Poly {
    type Output = Poly;

    fn add(self, rhs: Poly) -> Poly {
        let mut result = Poly::zero();
        for i in 0..N {
            result.coeffs[i] = self.coeffs[i] + rhs.coeffs[i];
        }
        result
    }
}

impl Add for &Poly {
    type Output = Poly;

    fn add(self, rhs: &Poly) -> Poly {
        let mut result = Poly::zero();
        for i in 0..N {
            result.coeffs[i] = self.coeffs[i] + rhs.coeffs[i];
        }
        result
    }
}

impl AddAssign for Poly {
    fn add_assign(&mut self, rhs: Poly) {
        for i in 0..N {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, rhs: &Poly) {
        for i in 0..N {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl Sub for Poly {
    type Output = Poly;

    fn sub(self, rhs: Poly) -> Poly {
        let mut result = Poly::zero();
        for i in 0..N {
            result.coeffs[i] = self.coeffs[i] - rhs.coeffs[i];
        }
        result
    }
}

impl Sub for &Poly {
    type Output = Poly;

    fn sub(self, rhs: &Poly) -> Poly {
        let mut result = Poly::zero();
        for i in 0..N {
            result.coeffs[i] = self.coeffs[i] - rhs.coeffs[i];
        }
        result
    }
}

impl SubAssign for Poly {
    fn sub_assign(&mut self, rhs: Poly) {
        for i in 0..N {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl SubAssign<&Poly> for Poly {
    fn sub_assign(&mut self, rhs: &Poly) {
        for i in 0..N {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl PartialEq for Poly {
    fn eq(&self, other: &Self) -> bool {
        // Compare in reduced form
        for i in 0..N {
            if freeze(self.coeffs[i]) != freeze(other.coeffs[i]) {
                return false;
            }
        }
        true
    }
}

impl Eq for Poly {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero() {
        let p = Poly::zero();
        assert!(p.is_zero());
        assert_eq!(p.norm_inf(), 0);
    }

    #[test]
    fn test_add_zero() {
        let mut p = Poly::zero();
        p.coeffs[0] = 12345;
        p.coeffs[100] = 67890;

        let result = &p + &Poly::zero();
        assert_eq!(result, p);
    }

    #[test]
    fn test_sub_self_is_zero() {
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = (i * 17) as i32 % Q;
        }

        let result = &p - &p;
        assert!(result.is_zero());
    }

    #[test]
    fn test_add_commutative() {
        let mut a = Poly::zero();
        let mut b = Poly::zero();
        for i in 0..N {
            a.coeffs[i] = (i * 7) as i32 % Q;
            b.coeffs[i] = (i * 11) as i32 % Q;
        }

        let ab = (&a + &b).reduced();
        let ba = (&b + &a).reduced();
        assert_eq!(ab, ba);
    }

    #[test]
    fn test_norm_inf() {
        let mut p = Poly::zero();
        p.coeffs[0] = 100;
        p.coeffs[1] = Q - 200; // -200 in centered form
        assert_eq!(p.norm_inf(), 200);
    }

    #[test]
    fn test_check_norm() {
        let mut p = Poly::zero();
        p.coeffs[0] = 100;
        assert!(p.check_norm(101));
        assert!(!p.check_norm(100));
        assert!(!p.check_norm(50));
    }

    #[test]
    fn test_ntt_roundtrip() {
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = (i * 13) as i32 % Q;
        }
        let original = p.clone();

        p.ntt();
        p.inv_ntt();

        // Should recover original (with possible reduction)
        for i in 0..N {
            let a = freeze(p.coeffs[i]);
            let b = freeze(original.coeffs[i]);
            assert!(
                (a - b).abs() < 2 || (a - b + Q).abs() < 2 || (a - b - Q).abs() < 2,
                "Mismatch at {}: {} vs {}",
                i,
                a,
                b
            );
        }
    }

    #[test]
    fn test_reduce() {
        let mut p = Poly::zero();
        p.coeffs[0] = Q + 100;
        p.coeffs[1] = -100;

        p.reduce();
        assert_eq!(p.coeffs[0], 100);
        assert_eq!(p.coeffs[1], Q - 100);
    }
}
