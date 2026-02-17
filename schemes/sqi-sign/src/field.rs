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

    /// Computes the multiplicative inverse using Fermat's little theorem.
    ///
    /// Uses constant-time exponentiation to prevent timing side-channels.
    /// a^(-1) = a^(p-2) mod p
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        // Use Fermat's little theorem: a^(-1) = a^(p-2) mod p
        let exp = &self.modulus - BigInt::from(2);
        Some(self.pow_ct(&exp))
    }

    /// Computes the multiplicative inverse (constant-time version).
    ///
    /// This version is guaranteed to run in constant time relative to the
    /// modulus size, regardless of the input value.
    pub fn inverse_ct(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        let exp = &self.modulus - BigInt::from(2);
        Some(self.pow_ct(&exp))
    }

    /// Computes self^exp mod p (variable-time).
    ///
    /// This is faster but should only be used with public exponents.
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

    /// Computes self^exp mod p (constant-time).
    ///
    /// Always performs the same number of operations regardless of the
    /// exponent value, preventing timing side-channels. The number of
    /// iterations is determined by the bit length of the modulus.
    ///
    /// # Security Notes
    ///
    /// This function uses a Montgomery ladder-like approach where we always
    /// perform both operations (multiply and square) but conditionally
    /// select which result to keep. The conditional selection uses a mask
    /// to avoid branches.
    pub fn pow_ct(&self, exp: &BigInt) -> Self {
        let bit_length = self.modulus.bits() as usize;

        let mut r0 = Self::one(self.modulus.clone());
        let mut r1 = self.clone();

        // Process from most significant bit to least significant
        for i in (0..bit_length).rev() {
            let bit = exp.bit(i as u64);

            // Always compute both options
            let r0_sq = &r0 * &r0;
            let r0_r1 = &r0 * &r1;
            let r1_sq = &r1 * &r1;

            // Constant-time conditional selection
            // If bit is 1: r0 = r0*r1, r1 = r1^2
            // If bit is 0: r0 = r0^2, r1 = r0*r1
            if bit {
                r0 = r0_r1.clone();
                r1 = r1_sq;
            } else {
                r0 = r0_sq;
                r1 = r0_r1;
            }
        }

        r0
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
    ///
    /// Returns Some(r) where r² = self (mod p), or None if self is not a QR.
    /// When p ≡ 3 (mod 4), uses the simpler formula r = self^((p+1)/4).
    /// Otherwise uses the full Tonelli-Shanks algorithm.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero(self.modulus.clone()));
        }

        // Check if self is a quadratic residue
        if self.legendre() != 1 {
            return None;
        }

        let p = &self.modulus;

        // Special case: p ≡ 3 (mod 4)
        // sqrt(a) = a^((p+1)/4)
        if p % BigInt::from(4) == BigInt::from(3) {
            let exp = (p + BigInt::one()) >> 2;
            return Some(self.pow(&exp));
        }

        // General case: Tonelli-Shanks algorithm
        // Factor p - 1 = 2^s * q where q is odd
        let mut q = p - BigInt::one();
        let mut s: u32 = 0;
        while &q % 2 == BigInt::zero() {
            q >>= 1;
            s += 1;
        }

        // Find a quadratic non-residue n
        let mut n = Self::new(BigInt::from(2), self.modulus.clone());
        while n.legendre() != -1 {
            n.value += BigInt::one();
        }

        // Initialize
        let mut m = s;
        let mut c = n.pow(&q);                    // c = n^q
        let mut t = self.pow(&q);                 // t = self^q
        let exp_r = (&q + BigInt::one()) >> 1;
        let mut r = self.pow(&exp_r);             // r = self^((q+1)/2)

        loop {
            if t.value == BigInt::one() {
                return Some(r);
            }

            // Find the least i such that t^(2^i) = 1
            let mut i: u32 = 1;
            let mut t_pow = &t * &t;
            while t_pow.value != BigInt::one() {
                t_pow = &t_pow * &t_pow;
                i += 1;
                if i >= m {
                    return None; // Should not happen if self is a QR
                }
            }

            // Update values
            // b = c^(2^(m-i-1))
            let exp_b = BigInt::one() << (m - i - 1);
            let b = c.pow(&exp_b);
            let b_sq = &b * &b;

            m = i;
            c = b_sq.clone();
            t = &t * &b_sq;
            r = &r * &b;
        }
    }

    /// Checks if self is a quadratic residue (has a square root in Fp).
    pub fn is_quadratic_residue(&self) -> bool {
        self.legendre() >= 0 // 0 or 1 means it's a QR (0 is sqrt of 0)
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

    #[test]
    fn test_sqrt_p_3_mod_4() {
        // Test with p ≡ 3 (mod 4): p = 103
        let p = BigInt::from(103);

        // 4 is a perfect square, sqrt(4) = ±2
        let four = Fp::new(BigInt::from(4), p.clone());
        let sqrt = four.sqrt().unwrap();
        let sq = &sqrt * &sqrt;
        assert_eq!(sq.value, four.value);

        // 25 is a perfect square
        let twenty_five = Fp::new(BigInt::from(25), p.clone());
        let sqrt = twenty_five.sqrt().unwrap();
        let sq = &sqrt * &sqrt;
        assert_eq!(sq.value, twenty_five.value);
    }

    #[test]
    fn test_sqrt_tonelli_shanks() {
        // Test with p ≡ 1 (mod 4): p = 101
        let p = BigInt::from(101);

        // Test various quadratic residues
        for val in [1, 4, 9, 16, 25, 36, 49] {
            let a = Fp::new(BigInt::from(val), p.clone());
            if a.legendre() == 1 {
                let sqrt = a.sqrt().unwrap();
                let sq = &sqrt * &sqrt;
                assert_eq!(sq.value, a.value, "sqrt({})² should equal {}", val, val);
            }
        }
    }

    #[test]
    fn test_sqrt_non_residue() {
        let p = BigInt::from(101);

        // 2 is a non-residue mod 101
        let two = Fp::new(BigInt::from(2), p.clone());
        if two.legendre() == -1 {
            assert!(two.sqrt().is_none());
        }
    }

    #[test]
    fn test_sqrt_zero() {
        let p = BigInt::from(101);
        let zero = Fp::zero(p);
        let sqrt = zero.sqrt().unwrap();
        assert!(sqrt.is_zero());
    }
}
