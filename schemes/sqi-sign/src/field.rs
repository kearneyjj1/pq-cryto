//! Prime field arithmetic Fp.
//!
//! This module implements arithmetic in the prime field Fp where p is the
//! SQI-SIGN prime. Operations include addition, subtraction, multiplication,
//! inversion, and square root computation.
//!
//! # Implementation Notes
//!
//! For educational clarity, this uses num-bigint for arbitrary precision.
//! A production implementation would use optimized multi-limb arithmetic
//! with Montgomery reduction.

use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::ops::{Add, Mul, Neg, Sub};

/// Element of the prime field Fp.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fp {
    /// The value, always reduced modulo p.
    pub value: BigInt,
    /// The modulus p.
    pub modulus: BigInt,
}

impl Fp {
    /// Creates a new field element.
    pub fn new(value: BigInt, modulus: BigInt) -> Self {
        let value = value % &modulus;
        let value = if value < BigInt::zero() {
            value + &modulus
        } else {
            value
        };
        Self { value, modulus }
    }

    /// Returns the zero element.
    pub fn zero(modulus: BigInt) -> Self {
        Self {
            value: BigInt::zero(),
            modulus,
        }
    }

    /// Returns the one element.
    pub fn one(modulus: BigInt) -> Self {
        Self {
            value: BigInt::one(),
            modulus,
        }
    }

    /// Returns true if this is zero.
    pub fn is_zero(&self) -> bool {
        self.value.is_zero()
    }

    /// Computes the multiplicative inverse using extended Euclidean algorithm.
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        // Use Fermat's little theorem: a^(-1) = a^(p-2) mod p
        let exp = &self.modulus - BigInt::from(2);
        Some(self.pow(&exp))
    }

    /// Computes self^exp mod p.
    pub fn pow(&self, exp: &BigInt) -> Self {
        if exp.is_zero() {
            return Self::one(self.modulus.clone());
        }

        let mut result = Self::one(self.modulus.clone());
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

    /// Computes the Legendre symbol (self/p).
    /// Returns 1 if self is a quadratic residue, -1 if not, 0 if self is 0.
    pub fn legendre(&self) -> i32 {
        if self.is_zero() {
            return 0;
        }
        let exp = (&self.modulus - BigInt::one()) >> 1;
        let result = self.pow(&exp);
        if result.value == BigInt::one() {
            1
        } else {
            -1
        }
    }

    /// Computes a square root if one exists (Tonelli-Shanks algorithm).
    pub fn sqrt(&self) -> Option<Self> {
        // TODO: Implement Tonelli-Shanks
        unimplemented!("Square root computation not yet implemented")
    }
}

impl Add for &Fp {
    type Output = Fp;

    fn add(self, other: &Fp) -> Fp {
        debug_assert_eq!(self.modulus, other.modulus);
        Fp::new(&self.value + &other.value, self.modulus.clone())
    }
}

impl Sub for &Fp {
    type Output = Fp;

    fn sub(self, other: &Fp) -> Fp {
        debug_assert_eq!(self.modulus, other.modulus);
        Fp::new(&self.value - &other.value, self.modulus.clone())
    }
}

impl Mul for &Fp {
    type Output = Fp;

    fn mul(self, other: &Fp) -> Fp {
        debug_assert_eq!(self.modulus, other.modulus);
        Fp::new(&self.value * &other.value, self.modulus.clone())
    }
}

impl Neg for &Fp {
    type Output = Fp;

    fn neg(self) -> Fp {
        Fp::new(-&self.value, self.modulus.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_modulus() -> BigInt {
        BigInt::from(101) // Small prime for testing
    }

    #[test]
    fn test_add() {
        let p = test_modulus();
        let a = Fp::new(BigInt::from(50), p.clone());
        let b = Fp::new(BigInt::from(60), p.clone());
        let c = &a + &b;
        assert_eq!(c.value, BigInt::from(9)); // 110 mod 101 = 9
    }

    #[test]
    fn test_mul() {
        let p = test_modulus();
        let a = Fp::new(BigInt::from(10), p.clone());
        let b = Fp::new(BigInt::from(11), p.clone());
        let c = &a * &b;
        assert_eq!(c.value, BigInt::from(9)); // 110 mod 101 = 9
    }

    #[test]
    fn test_inverse() {
        let p = test_modulus();
        let a = Fp::new(BigInt::from(5), p.clone());
        let a_inv = a.inverse().unwrap();
        let product = &a * &a_inv;
        assert_eq!(product.value, BigInt::one());
    }

    #[test]
    fn test_legendre() {
        let p = test_modulus();
        // 1 is always a QR
        let one = Fp::new(BigInt::one(), p.clone());
        assert_eq!(one.legendre(), 1);
    }
}
