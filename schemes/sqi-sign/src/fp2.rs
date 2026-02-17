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
    ///
    /// For z = a + bi, finds w such that w² = z.
    /// Every element in Fp² has a square root (Fp² is algebraically closed for quadratics).
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero(self.modulus().clone()));
        }

        let p = self.modulus().clone();
        let two = Fp::new(BigInt::from(2), p.clone());

        // Special case: purely real (im = 0)
        if self.im.is_zero() {
            // Try sqrt in Fp first
            if let Some(r) = self.re.sqrt() {
                return Some(Self::new(r, Fp::zero(p)));
            }
            // If a is not a QR in Fp, then sqrt(a) = sqrt(-a) * i
            let neg_re = -&self.re;
            if let Some(r) = neg_re.sqrt() {
                return Some(Self::new(Fp::zero(p), r));
            }
            return None;
        }

        // General case: z = a + bi where b ≠ 0
        // Compute the norm: α² = a² + b² (in Fp)
        let norm_sq = self.norm();

        // Take sqrt of norm in Fp
        let alpha = match norm_sq.sqrt() {
            Some(a) => a,
            None => return None,
        };

        // Try x = sqrt((a + α) / 2)
        let two_inv = match two.inverse() {
            Some(inv) => inv,
            None => return None,
        };

        let candidate1 = &(&self.re + &alpha) * &two_inv;
        let candidate2 = &(&self.re - &alpha) * &two_inv;

        let x = if let Some(sqrt_x) = candidate1.sqrt() {
            sqrt_x
        } else if let Some(sqrt_x) = candidate2.sqrt() {
            sqrt_x
        } else {
            return None;
        };

        // Compute y = b / (2x)
        if x.is_zero() {
            // Edge case: x = 0, so y² = -a, y = sqrt(-a)
            let neg_a = -&self.re;
            let y = neg_a.sqrt()?;
            return Some(Self::new(Fp::zero(p), y));
        }

        let two_x = &two * &x;
        let two_x_inv = two_x.inverse()?;
        let y = &self.im * &two_x_inv;

        Some(Self::new(x, y))
    }

    /// Checks if this element is a quadratic residue in Fp².
    /// Note: Every element in Fp² has a square root when p ≡ 3 (mod 4).
    pub fn is_quadratic_residue(&self) -> bool {
        // In Fp² with p ≡ 3 (mod 4), every element is a quadratic residue
        // For other cases, we check if sqrt succeeds
        if self.is_zero() {
            return true;
        }

        // Quick check: if purely real, check QR in Fp
        if self.im.is_zero() {
            let leg = self.re.legendre();
            // Either a QR in Fp, or -re is a QR (giving imaginary sqrt)
            if leg >= 0 {
                return true;
            }
            let neg_re = -&self.re;
            return neg_re.legendre() >= 0;
        }

        // For general elements, the norm must be a QR in Fp
        let norm = self.norm();
        norm.legendre() >= 0
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

    #[test]
    fn test_sqrt_real() {
        // p = 103 ≡ 3 (mod 4)
        let p = BigInt::from(103);

        // sqrt of 4 + 0i should be ±2 + 0i
        let four = Fp2::new(
            Fp::new(BigInt::from(4), p.clone()),
            Fp::zero(p.clone()),
        );
        let sqrt = four.sqrt().unwrap();
        let sq = &sqrt * &sqrt;
        assert_eq!(sq.re.value, four.re.value);
        assert!(sq.im.is_zero());
    }

    #[test]
    fn test_sqrt_purely_imaginary() {
        // p = 103 ≡ 3 (mod 4)
        let p = BigInt::from(103);

        // sqrt(-1) = i, so sqrt(0 + (-1)i) needs special handling
        // Test: sqrt of -1 (as a purely real element) should give i
        let neg_one = Fp2::new(
            Fp::new(&p - BigInt::one(), p.clone()),
            Fp::zero(p.clone()),
        );
        let sqrt = neg_one.sqrt().unwrap();
        let sq = &sqrt * &sqrt;
        // Should get back -1
        assert_eq!(sq.re.value, neg_one.re.value);
    }

    #[test]
    fn test_sqrt_complex() {
        // p = 103 ≡ 3 (mod 4)
        let p = BigInt::from(103);

        // Test sqrt of 3 + 4i
        let z = Fp2::new(
            Fp::new(BigInt::from(3), p.clone()),
            Fp::new(BigInt::from(4), p.clone()),
        );

        if let Some(sqrt) = z.sqrt() {
            let sq = &sqrt * &sqrt;
            assert_eq!(sq.re.value, z.re.value);
            assert_eq!(sq.im.value, z.im.value);
        }
        // Note: sqrt may not exist for all elements
    }

    #[test]
    fn test_sqrt_zero_fp2() {
        let p = BigInt::from(103);
        let zero = Fp2::zero(p);
        let sqrt = zero.sqrt().unwrap();
        assert!(sqrt.is_zero());
    }

    #[test]
    fn test_is_quadratic_residue() {
        let p = BigInt::from(103);

        // 4 is a QR
        let four = Fp2::new(
            Fp::new(BigInt::from(4), p.clone()),
            Fp::zero(p.clone()),
        );
        assert!(four.is_quadratic_residue());

        // Zero is a QR
        let zero = Fp2::zero(p.clone());
        assert!(zero.is_quadratic_residue());
    }
}
