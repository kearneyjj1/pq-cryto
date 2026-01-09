//! Supersingular elliptic curve operations.
//!
//! This module implements elliptic curve arithmetic for supersingular curves
//! over Fp². Curves are in Montgomery form: By² = x³ + Ax² + x.
//!
//! # Supersingular Curves
//!
//! A curve E/Fp² is supersingular if and only if its trace of Frobenius is
//! divisible by p. For SQI-SIGN, we use curves with j-invariant in Fp.

use crate::fp2::Fp2;
use num_bigint::BigInt;

/// A point on a Montgomery curve.
///
/// Uses projective coordinates (X : Z) where x = X/Z.
/// The point at infinity is represented as (1 : 0).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point {
    /// X coordinate.
    pub x: Fp2,
    /// Z coordinate.
    pub z: Fp2,
}

impl Point {
    /// Creates a new point from projective coordinates.
    pub fn new(x: Fp2, z: Fp2) -> Self {
        Self { x, z }
    }

    /// Returns the point at infinity.
    pub fn infinity(modulus: BigInt) -> Self {
        Self {
            x: Fp2::one(modulus.clone()),
            z: Fp2::zero(modulus),
        }
    }

    /// Returns true if this is the point at infinity.
    pub fn is_infinity(&self) -> bool {
        self.z.is_zero()
    }

    /// Normalizes the point so that Z = 1 (if not infinity).
    pub fn normalize(&self) -> Option<Self> {
        if self.is_infinity() {
            return None;
        }
        let z_inv = self.z.inverse()?;
        Some(Self {
            x: &self.x * &z_inv,
            z: Fp2::one(self.x.modulus().clone()),
        })
    }
}

/// A supersingular elliptic curve in Montgomery form.
///
/// By² = x³ + Ax² + x over Fp².
#[derive(Clone, Debug)]
pub struct Curve {
    /// The A coefficient in Montgomery form.
    pub a: Fp2,
    /// The modulus p.
    pub modulus: BigInt,
}

impl Curve {
    /// Creates a new curve with the given A coefficient.
    pub fn new(a: Fp2) -> Self {
        let modulus = a.modulus().clone();
        Self { a, modulus }
    }

    /// Computes the j-invariant of the curve.
    ///
    /// j(E) = 256(A² - 3)³ / (A² - 4)
    pub fn j_invariant(&self) -> Fp2 {
        // TODO: Implement j-invariant computation
        unimplemented!("j-invariant computation not yet implemented")
    }

    /// Montgomery differential addition.
    ///
    /// Given P, Q, and P-Q, computes P+Q.
    pub fn xadd(&self, p: &Point, q: &Point, diff: &Point) -> Point {
        // TODO: Implement Montgomery ladder differential addition
        unimplemented!("xADD not yet implemented")
    }

    /// Montgomery point doubling.
    pub fn xdbl(&self, p: &Point) -> Point {
        // TODO: Implement Montgomery doubling
        unimplemented!("xDBL not yet implemented")
    }

    /// Scalar multiplication using Montgomery ladder.
    pub fn scalar_mul(&self, k: &BigInt, p: &Point) -> Point {
        // TODO: Implement Montgomery ladder
        unimplemented!("Scalar multiplication not yet implemented")
    }

    /// Computes the order of a point (for testing).
    pub fn point_order(&self, p: &Point, max_order: &BigInt) -> Option<BigInt> {
        // TODO: Implement point order computation
        unimplemented!("Point order computation not yet implemented")
    }
}

/// The special starting curve E₀ with known endomorphism ring.
///
/// For SQI-SIGN, E₀ is typically the curve y² = x³ + x (j = 1728)
/// which has End(E₀) isomorphic to a maximal order in the quaternion
/// algebra ramified at p and ∞.
pub struct E0 {
    /// The curve.
    pub curve: Curve,
}

impl E0 {
    /// Creates E₀ for the given prime.
    pub fn new(modulus: BigInt) -> Self {
        // E₀: y² = x³ + x has A = 0 in Montgomery form
        let a = Fp2::zero(modulus.clone());
        Self {
            curve: Curve::new(a),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_infinity() {
        let p = BigInt::from(101);
        let inf = Point::infinity(p);
        assert!(inf.is_infinity());
    }

    #[test]
    fn test_e0_creation() {
        let p = BigInt::from(101);
        let e0 = E0::new(p);
        assert!(e0.curve.a.is_zero());
    }
}
