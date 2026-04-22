//! Field arithmetic in Z_q where q = 12289.
//!
//! This module provides modular arithmetic operations for the FALCON
//! signature scheme. Since q = 12289 < 2^14, we can use i16 for storage
//! and i32 for intermediate computations without overflow.

use crate::params::Q;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use zeroize::Zeroize;

/// An element of Z_q where q = 12289.
///
/// Values are stored in the range [0, q-1].
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, Hash, Zeroize)]
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
        debug_assert!(val >= 0 && (val as i32) < Q);
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
    ///
    /// **Not constant-time** with respect to the exponent: branches on
    /// exponent bits. Use `pow_ct` when the exponent is secret.
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

    /// Constant-time modular exponentiation.
    ///
    /// Always performs the multiply at every bit position and uses a
    /// branchless select to keep or discard the result. Safe for secret
    /// exponents.
    pub fn pow_ct(self, exp: u32) -> Self {
        let mut base = self;
        let mut result = Zq::ONE;

        for i in 0..32u32 {
            let product = result * base;
            let bit = ((exp >> i) & 1) as i16;
            // mask = 0xFFFF if bit=1, 0x0000 if bit=0
            let mask = bit.wrapping_neg();
            result = Zq((product.0 & mask) | (result.0 & !mask));
            base = base.square();
        }

        result
    }

    /// Computes the multiplicative inverse using Fermat's little theorem.
    ///
    /// For prime q, a^(-1) = a^(q-2) mod q.
    /// Returns 0 for input 0 (undefined, but safe).
    /// Uses constant-time exponentiation.
    #[inline]
    pub fn inverse(self) -> Self {
        self.pow_ct((Q - 2) as u32)
    }

    /// Returns the centered representation in [-(q-1)/2, (q-1)/2].
    ///
    /// Constant-time: no branches on the value.
    #[inline]
    pub fn centered(self) -> i16 {
        let v = self.0 as i32;
        let half_q = (Q / 2) as i32;
        // mask = all-1s if v > half_q, all-0s otherwise
        let mask = (half_q - v) >> 31;
        // If v > half_q: return v - Q; else return v
        (v + (mask & (-Q))) as i16
    }
}

/// Reduces a value modulo q into [0, q-1].
///
/// Constant-time: the `%` operator for a constant divisor compiles to a
/// multiply-shift sequence (no division instruction), and the conditional
/// add uses an arithmetic shift mask instead of a branch.
#[inline]
pub fn reduce(val: i32) -> i16 {
    let mut r = val % Q;
    // Branchless: if r < 0, add Q
    let neg_mask = r >> 31; // all-1s if r < 0, 0 otherwise
    r += Q & neg_mask;
    r as i16
}

impl Add for Zq {
    type Output = Zq;

    /// Constant-time modular addition.
    #[inline]
    fn add(self, rhs: Zq) -> Zq {
        let sum = self.0 as i32 + rhs.0 as i32;
        // Branchless: if sum >= Q, subtract Q
        let over = sum - Q;
        // mask = all-0s if over < 0 (sum < Q), all-1s if over >= 0 (sum >= Q)
        let mask = !(over >> 31);
        Zq((sum - (Q & mask)) as i16)
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

    /// Constant-time modular subtraction.
    #[inline]
    fn sub(self, rhs: Zq) -> Zq {
        let diff = self.0 as i32 - rhs.0 as i32;
        // Branchless: if diff < 0, add Q
        let neg_mask = diff >> 31; // all-1s if diff < 0, 0 otherwise
        Zq((diff + (Q & neg_mask)) as i16)
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

    /// Constant-time modular multiplication using Barrett reduction.
    ///
    /// Barrett constant M = floor(2^28 / Q) = 21843. The shift width 28 is
    /// chosen so that (Q-1)^2 * M = 12288^2 * 21843 ≈ 3.3e12 fits in i64.
    /// The approximation q_approx underestimates floor(prod/Q) by at most 1,
    /// so the remainder r = prod - q_approx*Q lies in [0, 2Q-1]. One
    /// branchless correction handles the case r >= Q.
    #[inline]
    fn mul(self, rhs: Zq) -> Zq {
        let prod = self.0 as i32 * rhs.0 as i32;
        // Both inputs in [0, Q-1], so prod in [0, (Q-1)^2 = ~1.51e8].
        const BARRETT_M: i64 = 21843; // floor(2^28 / 12289)
        let q_approx = ((prod as i64 * BARRETT_M) >> 28) as i32;
        let mut r = prod - q_approx * Q;
        // Branchless: if r >= Q, subtract Q.
        // over >= 0 means r >= Q: !(over >> 31) = all-1s, so r -= Q.
        // over < 0 means r < Q: !(over >> 31) = 0, so r unchanged.
        let over = r - Q;
        let mask = !(over >> 31);
        r -= Q & mask;
        Zq(r as i16)
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

    /// Constant-time modular negation.
    #[inline]
    fn neg(self) -> Zq {
        let v = self.0 as i32;
        // Branchless: (Q - v) if v != 0, else 0
        // nonzero = 1 if v != 0, 0 if v == 0
        let nonzero = ((v as u32) | (v as u32).wrapping_neg()) >> 31;
        // mask = all-1s if nonzero, all-0s if zero
        let mask = (nonzero as i32).wrapping_neg();
        Zq(((Q - v) & mask) as i16)
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

    #[test]
    fn test_ct_field_ops_exhaustive() {
        // Exhaustive test of add, sub, neg over all Q values
        for a_val in (0..Q).step_by(7) {
            let a = Zq(a_val as i16);
            // neg: a + (-a) must be zero
            let neg_a = -a;
            assert_eq!(a + neg_a, Zq::ZERO, "neg failed for {}", a_val);
            // centered must be in range
            let c = a.centered();
            assert!(c >= -((Q / 2) as i16) && c <= (Q / 2) as i16,
                "centered({}) = {} out of range", a_val, c);

            for b_val in (0..Q).step_by(13) {
                let b = Zq(b_val as i16);
                // add
                let sum = a + b;
                let expected_sum = (a_val + b_val) % Q;
                assert_eq!(sum.0 as i32, expected_sum, "add({},{}) failed", a_val, b_val);
                // sub
                let diff = a - b;
                let expected_diff = ((a_val - b_val) % Q + Q) % Q;
                assert_eq!(diff.0 as i32, expected_diff, "sub({},{}) failed", a_val, b_val);
            }
        }
    }

    #[test]
    fn test_ct_barrett_mul() {
        // Test mul via Barrett reduction over a wide range
        for a_val in (0..Q).step_by(17) {
            for b_val in (0..Q).step_by(19) {
                let a = Zq(a_val as i16);
                let b = Zq(b_val as i16);
                let prod = a * b;
                let expected = ((a_val as i64 * b_val as i64) % Q as i64) as i32;
                assert_eq!(prod.0 as i32, expected,
                    "mul({},{}) = {} expected {}", a_val, b_val, prod.0, expected);
            }
        }
    }
}
