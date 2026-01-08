//! Field arithmetic in Z_q where q = 12289.
//!
//! This module provides modular arithmetic operations for the FALCON
//! signature scheme. Since q = 12289 < 2^14, we can use i16 for storage
//! and i32 for intermediate computations without overflow.

use crate::params::Q;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// An element of Z_q where q = 12289.
///
/// Values are stored in the range [0, q-1].
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, Hash)]
pub struct Zq(pub i16);

impl Zq {
    /// The additive identity (zero).
    pub const ZERO: Zq = Zq(0);

    /// The multiplicative identity (one).
    pub const ONE: Zq = Zq(1);

    /// Creates a new field element, reducing modulo q.
    #[inline]
    pub fn new(val: i32) -> Self {
        Zq(reduce(val))
    }

    /// Creates a new field element from an i16 without reduction.
    /// The caller must ensure the value is already in [0, q-1].
    #[inline]
    pub const fn from_i16_unchecked(val: i16) -> Self {
        Zq(val)
    }

    /// Returns the value as an i16.
    #[inline]
    pub const fn value(self) -> i16 {
        self.0
    }

    /// Returns the value as an i32.
    #[inline]
    pub const fn as_i32(self) -> i32 {
        self.0 as i32
    }

    /// Returns true if this is zero.
    #[inline]
    pub const fn is_zero(self) -> bool {
        self.0 == 0
    }

    /// Computes the square of this element.
    #[inline]
    pub fn square(self) -> Self {
        let v = self.0 as i32;
        Zq::new(v * v)
    }

    /// Computes self^exp using binary exponentiation.
    pub fn pow(self, mut exp: u32) -> Self {
        if exp == 0 {
            return Zq::ONE;
        }

        let mut base = self;
        let mut result = Zq::ONE;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base;
            }
            base = base.square();
            exp >>= 1;
        }

        result
    }

    /// Computes the multiplicative inverse using Fermat's little theorem.
    ///
    /// For prime q, a^(-1) = a^(q-2) mod q.
    /// Returns 0 for input 0 (undefined, but safe).
    #[inline]
    pub fn inverse(self) -> Self {
        self.pow((Q - 2) as u32)
    }

    /// Returns the centered representation in [-(q-1)/2, (q-1)/2].
    #[inline]
    pub fn centered(self) -> i16 {
        let v = self.0;
        let half_q = (Q / 2) as i16;
        if v > half_q {
            v - Q as i16
        } else {
            v
        }
    }
}

/// Reduces a value modulo q into [0, q-1].
#[inline]
pub fn reduce(val: i32) -> i16 {
    let q = Q;
    let mut r = val % q;
    if r < 0 {
        r += q;
    }
    r as i16
}

/// Reduces a 64-bit value modulo q.
#[inline]
pub fn reduce64(val: i64) -> i16 {
    let q = Q as i64;
    let mut r = val % q;
    if r < 0 {
        r += q;
    }
    r as i16
}

impl Add for Zq {
    type Output = Zq;

    #[inline]
    fn add(self, rhs: Zq) -> Zq {
        let sum = self.0 as i32 + rhs.0 as i32;
        if sum >= Q {
            Zq((sum - Q) as i16)
        } else {
            Zq(sum as i16)
        }
    }
}

impl AddAssign for Zq {
    #[inline]
    fn add_assign(&mut self, rhs: Zq) {
        *self = *self + rhs;
    }
}

impl Sub for Zq {
    type Output = Zq;

    #[inline]
    fn sub(self, rhs: Zq) -> Zq {
        let diff = self.0 as i32 - rhs.0 as i32;
        if diff < 0 {
            Zq((diff + Q) as i16)
        } else {
            Zq(diff as i16)
        }
    }
}

impl SubAssign for Zq {
    #[inline]
    fn sub_assign(&mut self, rhs: Zq) {
        *self = *self - rhs;
    }
}

impl Mul for Zq {
    type Output = Zq;

    #[inline]
    fn mul(self, rhs: Zq) -> Zq {
        let prod = self.0 as i32 * rhs.0 as i32;
        Zq(reduce(prod))
    }
}

impl MulAssign for Zq {
    #[inline]
    fn mul_assign(&mut self, rhs: Zq) {
        *self = *self * rhs;
    }
}

impl Neg for Zq {
    type Output = Zq;

    #[inline]
    fn neg(self) -> Zq {
        if self.0 == 0 {
            Zq::ZERO
        } else {
            Zq((Q as i16) - self.0)
        }
    }
}

impl From<i16> for Zq {
    #[inline]
    fn from(val: i16) -> Self {
        Zq::new(val as i32)
    }
}

impl From<i32> for Zq {
    #[inline]
    fn from(val: i32) -> Self {
        Zq::new(val)
    }
}

impl From<u16> for Zq {
    #[inline]
    fn from(val: u16) -> Self {
        Zq::new(val as i32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_identity() {
        let a = Zq::new(1234);
        assert_eq!(a + Zq::ZERO, a);
    }

    #[test]
    fn test_add_commutative() {
        let a = Zq::new(1234);
        let b = Zq::new(5678);
        assert_eq!(a + b, b + a);
    }

    #[test]
    fn test_add_wraps() {
        let a = Zq::new(Q - 1);
        let b = Zq::new(2);
        assert_eq!(a + b, Zq::new(1));
    }

    #[test]
    fn test_sub_self_is_zero() {
        let a = Zq::new(1234);
        assert_eq!(a - a, Zq::ZERO);
    }

    #[test]
    fn test_sub_wraps() {
        let a = Zq::new(1);
        let b = Zq::new(3);
        assert_eq!(a - b, Zq::new(Q - 2));
    }

    #[test]
    fn test_mul_identity() {
        let a = Zq::new(1234);
        assert_eq!(a * Zq::ONE, a);
    }

    #[test]
    fn test_mul_zero() {
        let a = Zq::new(1234);
        assert_eq!(a * Zq::ZERO, Zq::ZERO);
    }

    #[test]
    fn test_mul_commutative() {
        let a = Zq::new(1234);
        let b = Zq::new(5678);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn test_neg() {
        let a = Zq::new(1234);
        let neg_a = -a;
        assert_eq!(a + neg_a, Zq::ZERO);
    }

    #[test]
    fn test_neg_zero() {
        assert_eq!(-Zq::ZERO, Zq::ZERO);
    }

    #[test]
    fn test_square() {
        let a = Zq::new(123);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn test_pow() {
        let a = Zq::new(7);
        assert_eq!(a.pow(0), Zq::ONE);
        assert_eq!(a.pow(1), a);
        assert_eq!(a.pow(2), a * a);
        assert_eq!(a.pow(3), a * a * a);
    }

    #[test]
    fn test_inverse() {
        for val in [1, 2, 7, 123, 1000, 12288] {
            let a = Zq::new(val);
            let a_inv = a.inverse();
            assert_eq!(
                a * a_inv,
                Zq::ONE,
                "inverse failed for {}",
                val
            );
        }
    }

    #[test]
    fn test_inverse_one() {
        assert_eq!(Zq::ONE.inverse(), Zq::ONE);
    }

    #[test]
    fn test_centered() {
        assert_eq!(Zq::new(0).centered(), 0);
        assert_eq!(Zq::new(1).centered(), 1);
        assert_eq!(Zq::new(Q - 1).centered(), -1);
        assert_eq!(Zq::new(Q / 2).centered(), (Q / 2) as i16);
        assert_eq!(Zq::new(Q / 2 + 1).centered(), -((Q / 2) as i16));
    }

    #[test]
    fn test_reduce_negative() {
        assert_eq!(reduce(-1), (Q - 1) as i16);
        assert_eq!(reduce(-Q), 0);
        assert_eq!(reduce(-Q - 1), (Q - 1) as i16);
    }

    #[test]
    fn test_reduce_large() {
        assert_eq!(reduce(Q), 0);
        assert_eq!(reduce(Q + 1), 1);
        assert_eq!(reduce(2 * Q), 0);
        assert_eq!(reduce(100 * Q + 42), 42);
    }
}
