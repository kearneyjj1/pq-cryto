//! Polynomial vectors and matrices for ML-DSA.
//!
//! This module provides vector and matrix operations over polynomials
//! in R_q = Z_q[X]/(X^256 + 1).
//!
//! In ML-DSA:
//! - s1 is a vector of l polynomials
//! - s2 is a vector of k polynomials
//! - A is a k×l matrix of polynomials
//! - t is a vector of k polynomials

use crate::poly::Poly;
use std::ops::{Add, AddAssign, Sub, SubAssign};

/// A vector of polynomials.
#[derive(Clone, Debug)]
pub struct PolyVec {
    /// The polynomials in this vector.
    pub polys: Vec<Poly>,
}

impl PolyVec {
    /// Creates a zero vector of the given length.
    pub fn zero(len: usize) -> Self {
        PolyVec {
            polys: (0..len).map(|_| Poly::zero()).collect(),
        }
    }

    /// Creates a vector from a Vec of polynomials.
    pub fn from_polys(polys: Vec<Poly>) -> Self {
        PolyVec { polys }
    }

    /// Returns the length of this vector.
    pub fn len(&self) -> usize {
        self.polys.len()
    }

    /// Returns true if the vector is empty.
    pub fn is_empty(&self) -> bool {
        self.polys.is_empty()
    }

    /// Reduces all coefficients in all polynomials.
    pub fn reduce(&mut self) {
        for poly in &mut self.polys {
            poly.reduce();
        }
    }

    /// Returns a reduced copy of this vector.
    pub fn reduced(&self) -> Self {
        let mut v = self.clone();
        v.reduce();
        v
    }

    /// Applies forward NTT to all polynomials.
    pub fn ntt(&mut self) {
        for poly in &mut self.polys {
            poly.ntt();
        }
    }

    /// Applies inverse NTT to all polynomials.
    pub fn inv_ntt(&mut self) {
        for poly in &mut self.polys {
            poly.inv_ntt();
        }
    }

    /// Freezes all coefficients to [0, Q).
    pub fn freeze(&mut self) {
        for poly in &mut self.polys {
            poly.freeze();
        }
    }

    /// Computes the dot product of two vectors (in NTT domain).
    ///
    /// Both vectors must be in NTT domain and have the same length.
    /// Returns a polynomial in NTT domain.
    pub fn dot(&self, other: &PolyVec) -> Poly {
        assert_eq!(self.len(), other.len(), "Vector lengths must match");

        let mut result = Poly::zero();
        for (a, b) in self.polys.iter().zip(other.polys.iter()) {
            let prod = a.pointwise_mul(b);
            result += prod;
        }
        result
    }

    /// Computes the infinity norm of the vector.
    ///
    /// Returns the maximum infinity norm across all polynomials.
    pub fn norm_inf(&self) -> i32 {
        self.polys.iter().map(|p| p.norm_inf()).max().unwrap_or(0)
    }

    /// Checks if ||v||_∞ < bound.
    pub fn check_norm(&self, bound: i32) -> bool {
        self.norm_inf() < bound
    }

    /// Shifts all coefficients left by d bits.
    pub fn shift_left(&mut self, d: usize) {
        for poly in &mut self.polys {
            poly.shift_left(d);
        }
    }
}

impl Add for PolyVec {
    type Output = PolyVec;

    fn add(self, rhs: PolyVec) -> PolyVec {
        assert_eq!(self.len(), rhs.len(), "Vector lengths must match");
        PolyVec {
            polys: self
                .polys
                .into_iter()
                .zip(rhs.polys.into_iter())
                .map(|(a, b)| a + b)
                .collect(),
        }
    }
}

impl Add for &PolyVec {
    type Output = PolyVec;

    fn add(self, rhs: &PolyVec) -> PolyVec {
        assert_eq!(self.len(), rhs.len(), "Vector lengths must match");
        PolyVec {
            polys: self
                .polys
                .iter()
                .zip(rhs.polys.iter())
                .map(|(a, b)| a + b)
                .collect(),
        }
    }
}

impl AddAssign<&PolyVec> for PolyVec {
    fn add_assign(&mut self, rhs: &PolyVec) {
        assert_eq!(self.len(), rhs.len(), "Vector lengths must match");
        for (a, b) in self.polys.iter_mut().zip(rhs.polys.iter()) {
            *a += b;
        }
    }
}

impl Sub for PolyVec {
    type Output = PolyVec;

    fn sub(self, rhs: PolyVec) -> PolyVec {
        assert_eq!(self.len(), rhs.len(), "Vector lengths must match");
        PolyVec {
            polys: self
                .polys
                .into_iter()
                .zip(rhs.polys.into_iter())
                .map(|(a, b)| a - b)
                .collect(),
        }
    }
}

impl Sub for &PolyVec {
    type Output = PolyVec;

    fn sub(self, rhs: &PolyVec) -> PolyVec {
        assert_eq!(self.len(), rhs.len(), "Vector lengths must match");
        PolyVec {
            polys: self
                .polys
                .iter()
                .zip(rhs.polys.iter())
                .map(|(a, b)| a - b)
                .collect(),
        }
    }
}

impl SubAssign<&PolyVec> for PolyVec {
    fn sub_assign(&mut self, rhs: &PolyVec) {
        assert_eq!(self.len(), rhs.len(), "Vector lengths must match");
        for (a, b) in self.polys.iter_mut().zip(rhs.polys.iter()) {
            *a -= b;
        }
    }
}

/// A matrix of polynomials (k rows, l columns).
#[derive(Clone, Debug)]
pub struct PolyMatrix {
    /// Rows of the matrix, where each row is a PolyVec.
    pub rows: Vec<PolyVec>,
    /// Number of rows (k).
    pub nrows: usize,
    /// Number of columns (l).
    pub ncols: usize,
}

impl PolyMatrix {
    /// Creates a zero matrix of dimensions k × l.
    pub fn zero(k: usize, l: usize) -> Self {
        PolyMatrix {
            rows: (0..k).map(|_| PolyVec::zero(l)).collect(),
            nrows: k,
            ncols: l,
        }
    }

    /// Applies forward NTT to all polynomials.
    pub fn ntt(&mut self) {
        for row in &mut self.rows {
            row.ntt();
        }
    }

    /// Applies inverse NTT to all polynomials.
    pub fn inv_ntt(&mut self) {
        for row in &mut self.rows {
            row.inv_ntt();
        }
    }

    /// Reduces all coefficients in all polynomials.
    pub fn reduce(&mut self) {
        for row in &mut self.rows {
            row.reduce();
        }
    }

    /// Matrix-vector multiplication: A * v (in NTT domain).
    ///
    /// Both the matrix and vector must be in NTT domain.
    /// Returns a vector in NTT domain.
    pub fn mul_vec(&self, v: &PolyVec) -> PolyVec {
        assert_eq!(
            self.ncols,
            v.len(),
            "Matrix columns must match vector length"
        );

        let mut result = PolyVec::zero(self.nrows);
        for (i, row) in self.rows.iter().enumerate() {
            result.polys[i] = row.dot(v);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{N, Q};

    #[test]
    fn test_polyvec_zero() {
        let v = PolyVec::zero(4);
        assert_eq!(v.len(), 4);
        assert_eq!(v.norm_inf(), 0);
    }

    #[test]
    fn test_polyvec_add() {
        let mut a = PolyVec::zero(3);
        let mut b = PolyVec::zero(3);

        a.polys[0].coeffs[0] = 100;
        b.polys[0].coeffs[0] = 200;

        let c = &a + &b;
        assert_eq!(c.polys[0].coeffs[0], 300);
    }

    #[test]
    fn test_polyvec_sub_self_is_zero() {
        let mut v = PolyVec::zero(4);
        for i in 0..4 {
            for j in 0..N {
                v.polys[i].coeffs[j] = ((i * N + j) * 17) as i32 % Q;
            }
        }

        let result = &v - &v;
        for poly in &result.polys {
            assert!(poly.is_zero());
        }
    }

    #[test]
    fn test_polyvec_norm_inf() {
        let mut v = PolyVec::zero(2);
        v.polys[0].coeffs[0] = 100;
        v.polys[1].coeffs[50] = Q - 200; // -200 in centered form

        assert_eq!(v.norm_inf(), 200);
    }

    #[test]
    fn test_polymatrix_zero() {
        let m = PolyMatrix::zero(4, 5);
        assert_eq!(m.nrows, 4);
        assert_eq!(m.ncols, 5);
        assert_eq!(m.rows.len(), 4);
        assert_eq!(m.rows[0].len(), 5);
    }

    #[test]
    fn test_dot_product_zero() {
        let a = PolyVec::zero(3);
        let b = PolyVec::zero(3);

        let result = a.dot(&b);
        assert!(result.is_zero());
    }

    #[test]
    fn test_matrix_vec_mul_identity_like() {
        // Create a simple test where first row's first element is "1-like" in NTT
        let mut m = PolyMatrix::zero(2, 2);
        let mut v = PolyVec::zero(2);

        // Set some coefficients
        m.rows[0].polys[0].coeffs[0] = 1;
        v.polys[0].coeffs[0] = 100;

        // In NTT domain
        m.ntt();
        v.ntt();

        let result = m.mul_vec(&v);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_ntt_roundtrip() {
        let mut v = PolyVec::zero(3);
        for i in 0..3 {
            for j in 0..N {
                v.polys[i].coeffs[j] = ((i * N + j) * 13) as i32 % Q;
            }
        }
        let original_norm = v.norm_inf();

        v.ntt();
        v.inv_ntt();

        // Norm should be approximately preserved
        let final_norm = v.norm_inf();
        assert!(
            (original_norm - final_norm).abs() < 10,
            "Norm changed too much: {} -> {}",
            original_norm,
            final_norm
        );
    }
}
