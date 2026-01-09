//! Quadratic extension field Fp².
//!
//! This module implements arithmetic in Fp² = Fp[i]/(i² + 1), the quadratic
//! extension of Fp. Supersingular elliptic curves for SQI-SIGN are defined
//! over Fp².
//!
//! Elements are represented as a + bi where a, b ∈ Fp.

use crate::field::Fp;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::ops::{Add, Mul, Neg, Sub};

/// Element of the quadratic extension field Fp².
///
/// Represented as a + b*i where i² = -1.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fp2 {
    /// Real part (coefficient of 1).
    pub re: Fp,
    /// Imaginary part (coefficient of i).
    pub im: Fp,
}

impl Fp2 {
    /// Creates a new Fp² element.
    pub fn new(re: Fp, im: Fp) -> Self {
        debug_assert_eq!(re.modulus, im.modulus);
        Self { re, im }
    }

    /// Creates the zero element.
    pub fn zero(modulus: BigInt) -> Self {
        Self {
            re: Fp::zero(modulus.clone()),
            im: Fp::zero(modulus),
        }
    }

    /// Creates the one element.
    pub fn one(modulus: BigInt) -> Self {
        Self {
            re: Fp::one(modulus.clone()),
            im: Fp::zero(modulus),
        }
    }

    /// Creates i (the imaginary unit).
    pub fn i(modulus: BigInt) -> Self {
        Self {
            re: Fp::zero(modulus.clone()),
            im: Fp::one(modulus),
        }
    }

    /// Returns true if this is zero.
    pub fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }

    /// Returns the modulus.
    pub fn modulus(&self) -> &BigInt {
        &self.re.modulus
    }

    /// Computes the conjugate: conj(a + bi) = a - bi.
    pub fn conjugate(&self) -> Self {
        Self {
            re: self.re.clone(),
            im: -&self.im,
        }
    }

    /// Computes the norm: N(a + bi) = a² + b² (in Fp).
    pub fn norm(&self) -> Fp {
        &(&self.re * &self.re) + &(&self.im * &self.im)
    }

    /// Computes the multiplicative inverse.
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        // (a + bi)^(-1) = (a - bi) / (a² + b²)
        let n = self.norm();
        let n_inv = n.inverse()?;
        let conj = self.conjugate();
        Some(Self {
            re: &conj.re * &n_inv,
            im: &conj.im * &n_inv,
        })
    }

    /// Computes self^exp.
    pub fn pow(&self, exp: &BigInt) -> Self {
        if exp.is_zero() {
            return Self::one(self.modulus().clone());
        }

        let mut result = Self::one(self.modulus().clone());
        let mut base = self.clone();
        let mut e = exp.clone();

        while e > BigInt::zero() {
            if &e % 2 == BigInt::one() {
                result = &result * &base;
            }
            base = &base * &base;
            e >>= 1;
        }

        result
    }

    /// Computes a square root in Fp² if one exists.
    pub fn sqrt(&self) -> Option<Self> {
        // TODO: Implement square root in Fp²
        unimplemented!("Fp² square root not yet implemented")
    }
}

impl Add for &Fp2 {
    type Output = Fp2;

    fn add(self, other: &Fp2) -> Fp2 {
        Fp2 {
            re: &self.re + &other.re,
            im: &self.im + &other.im,
        }
    }
}

impl Sub for &Fp2 {
    type Output = Fp2;

    fn sub(self, other: &Fp2) -> Fp2 {
        Fp2 {
            re: &self.re - &other.re,
            im: &self.im - &other.im,
        }
    }
}

impl Mul for &Fp2 {
    type Output = Fp2;

    fn mul(self, other: &Fp2) -> Fp2 {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        Fp2 {
            re: &(&self.re * &other.re) - &(&self.im * &other.im),
            im: &(&self.re * &other.im) + &(&self.im * &other.re),
        }
    }
}

impl Neg for &Fp2 {
    type Output = Fp2;

    fn neg(self) -> Fp2 {
        Fp2 {
            re: -&self.re,
            im: -&self.im,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_modulus() -> BigInt {
        BigInt::from(101)
    }

    #[test]
    fn test_conjugate_norm() {
        let p = test_modulus();
        let a = Fp2::new(
            Fp::new(BigInt::from(3), p.clone()),
            Fp::new(BigInt::from(4), p.clone()),
        );
        let norm = a.norm();
        // 3² + 4² = 25
        assert_eq!(norm.value, BigInt::from(25));
    }

    #[test]
    fn test_inverse() {
        let p = test_modulus();
        let a = Fp2::new(
            Fp::new(BigInt::from(3), p.clone()),
            Fp::new(BigInt::from(4), p.clone()),
        );
        let a_inv = a.inverse().unwrap();
        let product = &a * &a_inv;
        assert!((&product.re.value - BigInt::one()).is_zero());
        assert!(product.im.value.is_zero());
    }

    #[test]
    fn test_i_squared() {
        let p = test_modulus();
        let i = Fp2::i(p.clone());
        let i_sq = &i * &i;
        // i² = -1
        assert_eq!(i_sq.re.value, p - BigInt::one());
        assert!(i_sq.im.is_zero());
    }
}
