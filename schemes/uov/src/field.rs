//! GF(2^8) finite field arithmetic using the AES polynomial.
//!
//! This module implements arithmetic in GF(256) with the irreducible polynomial
//! x^8 + x^4 + x^3 + x + 1 (0x11B), which is the same polynomial used in AES.

use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

/// An element of the finite field GF(2^8).
///
/// Elements are represented as bytes, with arithmetic defined by the
/// AES irreducible polynomial x^8 + x^4 + x^3 + x + 1.
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, Hash)]
pub struct F(pub u8);

impl F {
    /// The additive identity (zero element).
    pub const ZERO: F = F(0);

    /// The multiplicative identity (one element).
    pub const ONE: F = F(1);

    /// Creates a new field element from a byte.
    #[inline]
    pub const fn new(val: u8) -> Self {
        F(val)
    }

    /// Returns the underlying byte value.
    #[inline]
    pub const fn value(self) -> u8 {
        self.0
    }

    /// Returns true if this is the zero element.
    #[inline]
    pub const fn is_zero(self) -> bool {
        self.0 == 0
    }

    /// Computes the square of this element.
    #[inline]
    pub fn square(self) -> F {
        self * self
    }

    /// Computes the multiplicative inverse using Fermat's little theorem.
    ///
    /// In GF(2^8), a^255 = 1 for all non-zero a, so a^(-1) = a^254.
    /// We compute this as a^254 = a^(128+64+32+16+8+4+2).
    ///
    /// Returns F(0) for zero input (undefined, but safe).
    #[inline]
    pub fn inverse(self) -> F {
        // a^2, a^4, a^8, a^16, a^32, a^64, a^128
        let x2 = self.square();
        let x4 = x2.square();
        let x8 = x4.square();
        let x16 = x8.square();
        let x32 = x16.square();
        let x64 = x32.square();
        let x128 = x64.square();

        // a^254 = a^128 * a^64 * a^32 * a^16 * a^8 * a^4 * a^2
        x128 * x64 * x32 * x16 * x8 * x4 * x2
    }
}

/// Addition in GF(2^8) is XOR.
#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for F {
    type Output = F;

    #[inline]
    fn add(self, rhs: F) -> F {
        F(self.0 ^ rhs.0)
    }
}

#[allow(clippy::suspicious_op_assign_impl)]
impl AddAssign for F {
    #[inline]
    fn add_assign(&mut self, rhs: F) {
        self.0 ^= rhs.0;
    }
}

/// Subtraction in GF(2^8) is the same as addition (XOR).
#[allow(clippy::suspicious_arithmetic_impl)]
impl Sub for F {
    type Output = F;

    #[inline]
    fn sub(self, rhs: F) -> F {
        F(self.0 ^ rhs.0)
    }
}

#[allow(clippy::suspicious_op_assign_impl)]
impl SubAssign for F {
    #[inline]
    fn sub_assign(&mut self, rhs: F) {
        self.0 ^= rhs.0;
    }
}

/// Multiplication in GF(2^8) using the AES polynomial for reduction.
///
/// # Security
///
/// This implementation is constant-time: it always performs exactly 8 iterations
/// with the same operations, using masking instead of data-dependent branches.
impl Mul for F {
    type Output = F;

    #[inline]
    fn mul(self, rhs: F) -> F {
        let mut result = 0u8;
        let mut a = self.0;
        let mut b = rhs.0;

        // Fixed 8 iterations for constant-time execution
        for _ in 0..8 {
            // Mask is 0xFF if (b & 1) == 1, else 0x00
            let mask = 0u8.wrapping_sub(b & 1);
            result ^= a & mask;

            // Mask for high bit reduction
            let high_bit_mask = 0u8.wrapping_sub((a >> 7) & 1);
            a = (a << 1) ^ (0x1b & high_bit_mask);

            b >>= 1;
        }

        F(result)
    }
}

impl MulAssign for F {
    #[inline]
    fn mul_assign(&mut self, rhs: F) {
        *self = *self * rhs;
    }
}

impl From<u8> for F {
    #[inline]
    fn from(val: u8) -> Self {
        F(val)
    }
}

impl From<F> for u8 {
    #[inline]
    fn from(f: F) -> u8 {
        f.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_identity() {
        let a = F(0x57);
        assert_eq!(a + F::ZERO, a);
    }

    #[test]
    fn test_add_self_is_zero() {
        let a = F(0x57);
        assert_eq!(a + a, F::ZERO);
    }

    #[test]
    fn test_mul_identity() {
        let a = F(0x57);
        assert_eq!(a * F::ONE, a);
    }

    #[test]
    fn test_mul_zero() {
        let a = F(0x57);
        assert_eq!(a * F::ZERO, F::ZERO);
    }

    #[test]
    fn test_mul_known_value() {
        // Known AES multiplication: 0x57 * 0x83 = 0xC1
        assert_eq!(F(0x57) * F(0x83), F(0xC1));
    }

    #[test]
    fn test_inverse() {
        let a = F(0x57);
        let a_inv = a.inverse();
        assert_eq!(a * a_inv, F::ONE);
    }

    #[test]
    fn test_inverse_of_one() {
        assert_eq!(F::ONE.inverse(), F::ONE);
    }

    #[test]
    fn test_sub_equals_add() {
        let a = F(0x57);
        let b = F(0x83);
        assert_eq!(a - b, a + b);
    }
}
