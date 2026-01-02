//! Field arithmetic for Zq where q = 8380417.
//!
//! This module implements modular arithmetic in the ring Z_q used throughout ML-DSA.
//! The prime q = 2^23 - 2^13 + 1 = 8380417 has special structure that enables
//! efficient NTT operations.

use crate::params::Q;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// An element of the ring Z_q where q = 8380417.
///
/// Elements are stored as signed 32-bit integers in the range [0, q).
/// The signed representation simplifies reduction and comparison operations.
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, Hash)]
pub struct Zq(pub i32);

impl Zq {
    /// The additive identity (zero element).
    pub const ZERO: Zq = Zq(0);

    /// The multiplicative identity (one element).
    pub const ONE: Zq = Zq(1);

    /// Creates a new field element, reducing modulo q.
    #[inline]
    pub fn new(val: i32) -> Self {
        Zq(reduce(val))
    }

    /// Creates a field element from an i32 without reduction.
    ///
    /// # Safety
    /// The caller must ensure that val is in [0, Q).
    #[inline]
    pub const fn from_unreduced(val: i32) -> Self {
        Zq(val)
    }

    /// Returns the value as an i32 in [0, q).
    #[inline]
    pub const fn value(self) -> i32 {
        self.0
    }

    /// Returns true if this is the zero element.
    #[inline]
    pub const fn is_zero(self) -> bool {
        self.0 == 0
    }

    /// Reduces the value to be in [0, q).
    #[inline]
    pub fn reduce(self) -> Self {
        Zq(reduce(self.0))
    }

    /// Computes the negation: -a mod q.
    #[inline]
    pub fn neg(self) -> Self {
        if self.0 == 0 {
            Zq::ZERO
        } else {
            Zq(Q - self.0)
        }
    }

    /// Computes the square of this element.
    #[inline]
    pub fn square(self) -> Self {
        // Use i64 to avoid overflow
        let product = (self.0 as i64) * (self.0 as i64);
        Zq(reduce64(product))
    }

    /// Computes self^exp mod q using square-and-multiply.
    pub fn pow(self, mut exp: u32) -> Self {
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
    /// In Z_q, a^(q-1) = 1 for all non-zero a, so a^(-1) = a^(q-2).
    ///
    /// Returns Zq(0) for zero input (undefined, but safe).
    pub fn inverse(self) -> Self {
        self.pow((Q - 2) as u32)
    }

    /// Converts to centered representation in [-(q-1)/2, (q-1)/2].
    #[inline]
    pub fn to_centered(self) -> i32 {
        let half_q = (Q - 1) / 2;
        if self.0 > half_q {
            self.0 - Q
        } else {
            self.0
        }
    }

    /// Creates from centered representation.
    #[inline]
    pub fn from_centered(val: i32) -> Self {
        if val < 0 {
            Zq(val + Q)
        } else {
            Zq(val)
        }
    }
}

/// Reduces a value to [0, q).
#[inline]
fn reduce(mut a: i32) -> i32 {
    // Handle negative values
    a = a % Q;
    if a < 0 {
        a += Q;
    }
    a
}

/// Reduces an i64 product to i32 in [0, q).
#[inline]
fn reduce64(a: i64) -> i32 {
    let r = (a % (Q as i64)) as i32;
    if r < 0 {
        r + Q
    } else {
        r
    }
}

impl Add for Zq {
    type Output = Zq;

    #[inline]
    fn add(self, rhs: Zq) -> Zq {
        let sum = self.0 + rhs.0;
        // Conditional subtraction for reduction
        if sum >= Q {
            Zq(sum - Q)
        } else {
            Zq(sum)
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
        let diff = self.0 - rhs.0;
        if diff < 0 {
            Zq(diff + Q)
        } else {
            Zq(diff)
        }
    }
}

impl SubAssign for Zq {
    #[inline]
    fn sub_assign(&mut self, rhs: Zq) {
        *self = *self - rhs;
    }
}

impl Neg for Zq {
    type Output = Zq;

    #[inline]
    fn neg(self) -> Zq {
        self.neg()
    }
}

impl Mul for Zq {
    type Output = Zq;

    #[inline]
    fn mul(self, rhs: Zq) -> Zq {
        // Use i64 to avoid overflow
        let product = (self.0 as i64) * (rhs.0 as i64);
        Zq(reduce64(product))
    }
}

impl MulAssign for Zq {
    #[inline]
    fn mul_assign(&mut self, rhs: Zq) {
        *self = *self * rhs;
    }
}

impl From<i32> for Zq {
    #[inline]
    fn from(val: i32) -> Self {
        Zq::new(val)
    }
}

impl From<u32> for Zq {
    #[inline]
    fn from(val: u32) -> Self {
        Zq::new((val % (Q as u32)) as i32)
    }
}

impl From<Zq> for i32 {
    #[inline]
    fn from(z: Zq) -> i32 {
        z.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_identity() {
        let a = Zq::new(12345);
        assert_eq!(a + Zq::ZERO, a);
    }

    #[test]
    fn test_add_commutative() {
        let a = Zq::new(12345);
        let b = Zq::new(67890);
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
        let a = Zq::new(12345);
        assert_eq!(a - a, Zq::ZERO);
    }

    #[test]
    fn test_sub_wraps() {
        let a = Zq::new(1);
        let b = Zq::new(2);
        assert_eq!(a - b, Zq::new(Q - 1));
    }

    #[test]
    fn test_neg() {
        let a = Zq::new(12345);
        assert_eq!(a + (-a), Zq::ZERO);
    }

    #[test]
    fn test_neg_zero() {
        assert_eq!(-Zq::ZERO, Zq::ZERO);
    }

    #[test]
    fn test_mul_identity() {
        let a = Zq::new(12345);
        assert_eq!(a * Zq::ONE, a);
    }

    #[test]
    fn test_mul_zero() {
        let a = Zq::new(12345);
        assert_eq!(a * Zq::ZERO, Zq::ZERO);
    }

    #[test]
    fn test_mul_commutative() {
        let a = Zq::new(12345);
        let b = Zq::new(67890);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn test_mul_overflow() {
        // Test with values that would overflow i32 when multiplied
        let a = Zq::new(Q - 1);
        let b = Zq::new(Q - 1);
        let result = a * b;
        // (Q-1)^2 mod Q = 1
        assert_eq!(result, Zq::ONE);
    }

    #[test]
    fn test_square() {
        let a = Zq::new(12345);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn test_inverse() {
        let a = Zq::new(12345);
        let a_inv = a.inverse();
        assert_eq!(a * a_inv, Zq::ONE);
    }

    #[test]
    fn test_inverse_of_one() {
        assert_eq!(Zq::ONE.inverse(), Zq::ONE);
    }

    #[test]
    fn test_pow() {
        let a = Zq::new(2);
        assert_eq!(a.pow(10), Zq::new(1024));
    }

    #[test]
    fn test_pow_zero() {
        let a = Zq::new(12345);
        assert_eq!(a.pow(0), Zq::ONE);
    }

    #[test]
    fn test_centered_representation() {
        // Small positive stays positive
        let a = Zq::new(100);
        assert_eq!(a.to_centered(), 100);

        // Large value becomes negative
        let b = Zq::new(Q - 100);
        assert_eq!(b.to_centered(), -100);
    }

    #[test]
    fn test_from_centered() {
        let a = Zq::from_centered(-100);
        assert_eq!(a, Zq::new(Q - 100));

        let b = Zq::from_centered(100);
        assert_eq!(b, Zq::new(100));
    }

    #[test]
    fn test_reduce_negative() {
        let a = Zq::new(-1);
        assert_eq!(a, Zq::new(Q - 1));
    }

    #[test]
    fn test_reduce_large() {
        let a = Zq::new(Q + 100);
        assert_eq!(a, Zq::new(100));
    }

    #[test]
    fn test_q_value() {
        // Verify q = 2^23 - 2^13 + 1
        assert_eq!(Q, (1 << 23) - (1 << 13) + 1);
        assert_eq!(Q, 8380417);
    }
}
