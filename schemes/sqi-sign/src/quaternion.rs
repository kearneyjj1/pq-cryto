//! Quaternion algebra computations.
//!
//! This module implements arithmetic in the quaternion algebra B_{p,∞}
//! ramified at p and infinity. This algebra is central to SQI-SIGN as
//! the endomorphism ring of supersingular curves embeds into it.
//!
//! # Structure
//!
//! B_{p,∞} = Q + Qi + Qj + Qk where:
//! - i² = -1
//! - j² = -p
//! - k = ij = -ji
//!
//! For computational purposes, we work with orders (lattices) in B_{p,∞}.

use num_bigint::BigInt;
use num_traits::Zero;

/// An element of the quaternion algebra B_{p,∞}.
///
/// Represented as α = a + bi + cj + dk where a, b, c, d ∈ Q.
/// We store these as rationals (numerator/denominator pairs).
#[derive(Clone, Debug)]
pub struct Quaternion {
    /// Coefficient of 1.
    pub a: Rational,
    /// Coefficient of i.
    pub b: Rational,
    /// Coefficient of j.
    pub c: Rational,
    /// Coefficient of k.
    pub d: Rational,
    /// The prime p (j² = -p).
    pub p: BigInt,
}

/// A rational number a/b.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rational {
    /// Numerator.
    pub num: BigInt,
    /// Denominator (always positive).
    pub den: BigInt,
}

impl Rational {
    /// Creates a new rational, reducing to lowest terms.
    pub fn new(num: BigInt, den: BigInt) -> Self {
        // TODO: Reduce to lowest terms using GCD
        Self { num, den }
    }

    /// Creates a rational from an integer.
    pub fn from_int(n: BigInt) -> Self {
        Self {
            num: n,
            den: BigInt::from(1),
        }
    }

    /// Returns zero.
    pub fn zero() -> Self {
        Self {
            num: BigInt::zero(),
            den: BigInt::from(1),
        }
    }

    /// Returns true if this is zero.
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

impl Quaternion {
    /// Creates a new quaternion.
    pub fn new(a: Rational, b: Rational, c: Rational, d: Rational, p: BigInt) -> Self {
        Self { a, b, c, d, p }
    }

    /// Creates the zero quaternion.
    pub fn zero(p: BigInt) -> Self {
        Self {
            a: Rational::zero(),
            b: Rational::zero(),
            c: Rational::zero(),
            d: Rational::zero(),
            p,
        }
    }

    /// Creates the identity quaternion (1).
    pub fn one(p: BigInt) -> Self {
        Self {
            a: Rational::from_int(BigInt::from(1)),
            b: Rational::zero(),
            c: Rational::zero(),
            d: Rational::zero(),
            p,
        }
    }

    /// Computes the conjugate: conj(a + bi + cj + dk) = a - bi - cj - dk.
    pub fn conjugate(&self) -> Self {
        // TODO: Implement conjugation
        unimplemented!("Quaternion conjugation not yet implemented")
    }

    /// Computes the reduced norm: nrd(α) = α * conj(α).
    pub fn reduced_norm(&self) -> Rational {
        // nrd(a + bi + cj + dk) = a² + b² + pc² + pd²
        // TODO: Implement reduced norm
        unimplemented!("Reduced norm not yet implemented")
    }

    /// Computes the reduced trace: trd(α) = α + conj(α) = 2a.
    pub fn reduced_trace(&self) -> Rational {
        // TODO: Implement reduced trace
        unimplemented!("Reduced trace not yet implemented")
    }

    /// Multiplies two quaternions.
    pub fn mul(&self, other: &Quaternion) -> Quaternion {
        // TODO: Implement quaternion multiplication
        // Remember: ij = k, ji = -k, i² = -1, j² = -p
        unimplemented!("Quaternion multiplication not yet implemented")
    }

    /// Computes the inverse (if it exists).
    pub fn inverse(&self) -> Option<Quaternion> {
        // α^(-1) = conj(α) / nrd(α)
        // TODO: Implement quaternion inversion
        unimplemented!("Quaternion inversion not yet implemented")
    }
}

/// A maximal order in B_{p,∞}.
///
/// For SQI-SIGN, we work with the maximal order O₀ corresponding to End(E₀).
#[derive(Clone, Debug)]
pub struct MaximalOrder {
    /// Basis for the order: O = Z⟨b₁, b₂, b₃, b₄⟩.
    pub basis: [Quaternion; 4],
    /// The prime p.
    pub p: BigInt,
}

impl MaximalOrder {
    /// Creates the standard maximal order O₀ for the given prime.
    ///
    /// This is the order corresponding to End(E₀) where E₀: y² = x³ + x.
    pub fn standard(p: BigInt) -> Self {
        // TODO: Compute the standard maximal order
        unimplemented!("Standard maximal order not yet implemented")
    }

    /// Tests if a quaternion is in this order.
    pub fn contains(&self, q: &Quaternion) -> bool {
        // TODO: Check if q can be expressed as integer combination of basis
        unimplemented!("Order membership test not yet implemented")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rational_zero() {
        let z = Rational::zero();
        assert!(z.is_zero());
    }

    #[test]
    fn test_quaternion_creation() {
        let p = BigInt::from(101);
        let q = Quaternion::one(p.clone());
        assert!(!q.a.is_zero());
        assert!(q.b.is_zero());
    }
}
