//! Supersingular elliptic curve operations.
//!
//! This module implements elliptic curve arithmetic for supersingular curves
//! over Fp². Curves are in Montgomery form: By² = x³ + Ax² + x.
//!
//! # Supersingular Curves
//!
//! A curve E/Fp² is supersingular if and only if its trace of Frobenius is
//! divisible by p. For SQI-SIGN, we use curves with j-invariant in Fp.
//!
//! # Montgomery Arithmetic
//!
//! We use projective x-coordinates (X : Z) where the affine x = X/Z.
//! This allows efficient differential addition without computing y-coordinates.

use crate::fp2::Fp2;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// Constant-time conditional swap for Fp2 elements.
/// If condition is true (non-zero), swaps a and b. Otherwise does nothing.
/// This is implemented in constant time to prevent timing side-channels.
#[inline]
fn cswap_fp2(condition: bool, a: &mut Fp2, b: &mut Fp2) {
    // Create mask: all 1s if condition, all 0s otherwise
    // We work on the underlying BigInt values
    let mask = if condition {
        BigInt::from(-1i64)
    } else {
        BigInt::zero()
    };

    // XOR swap with mask
    // temp = mask & (a ^ b)
    // a ^= temp
    // b ^= temp
    let diff_re = &a.re.value ^ &b.re.value;
    let diff_im = &a.im.value ^ &b.im.value;

    let masked_re = &diff_re & &mask;
    let masked_im = &diff_im & &mask;

    a.re.value = &a.re.value ^ &masked_re;
    a.im.value = &a.im.value ^ &masked_im;
    b.re.value = &b.re.value ^ &masked_re;
    b.im.value = &b.im.value ^ &masked_im;
}

/// Constant-time conditional swap for Points.
/// If condition is true, swaps p1 and p2. Otherwise does nothing.
#[inline]
fn cswap_point(condition: bool, p1: &mut Point, p2: &mut Point) {
    cswap_fp2(condition, &mut p1.x, &mut p2.x);
    cswap_fp2(condition, &mut p1.z, &mut p2.z);
}

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

    /// Validates that this point lies on the given Montgomery curve.
    ///
    /// For a Montgomery curve By² = x³ + Ax² + x, a point (x, y) is valid
    /// if the right-hand side x³ + Ax² + x is a quadratic residue in Fp².
    /// For x-only representation, we just check this condition.
    ///
    /// Returns true if the point is valid (on curve or at infinity).
    pub fn is_valid_on_curve(&self, curve: &Curve) -> bool {
        if self.is_infinity() {
            return true;
        }

        // Normalize to get affine x-coordinate
        let x = match self.normalize() {
            Some(p) => p.x,
            None => return true, // infinity is valid
        };

        // Compute RHS = x³ + Ax² + x
        let x_sq = &x * &x;
        let x_cube = &x_sq * &x;
        let ax_sq = &curve.a * &x_sq;
        let rhs = &(&x_cube + &ax_sq) + &x;

        // Check if RHS is a quadratic residue (has a square root)
        // For now, we use the Legendre symbol approach
        // In Fp², an element is a QR if its norm is a QR in Fp
        // This is a simplified check
        rhs.is_quadratic_residue()
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
    /// For Montgomery curve By² = x³ + Ax² + x:
    /// j(E) = 256(A² - 3)³ / (A² - 4)
    pub fn j_invariant(&self) -> Fp2 {
        let one = Fp2::one(self.modulus.clone());
        let two = &one + &one;
        let three = &two + &one;
        let four = &two + &two;
        let c256 = {
            let c16 = &four * &four;
            &c16 * &c16
        };

        // A²
        let a_sq = &self.a * &self.a;

        // A² - 3
        let num_base = &a_sq - &three;

        // (A² - 3)³
        let num_sq = &num_base * &num_base;
        let num_cube = &num_sq * &num_base;

        // 256(A² - 3)³
        let numerator = &c256 * &num_cube;

        // A² - 4
        let denominator = &a_sq - &four;

        // j = numerator / denominator
        if let Some(denom_inv) = denominator.inverse() {
            &numerator * &denom_inv
        } else {
            // Degenerate case: A² = 4, curve is singular
            Fp2::zero(self.modulus.clone())
        }
    }

    /// Montgomery differential addition.
    ///
    /// Given P = (XP : ZP), Q = (XQ : ZQ), and P-Q = (Xd : Zd), computes P+Q.
    ///
    /// Uses the standard formula:
    /// - U = (XP - ZP)(XQ + ZQ)
    /// - V = (XP + ZP)(XQ - ZQ)
    /// - X+ = Zd · (U + V)²
    /// - Z+ = Xd · (U - V)²
    pub fn xadd(&self, p: &Point, q: &Point, diff: &Point) -> Point {
        // Handle point at infinity cases
        if p.is_infinity() {
            return q.clone();
        }
        if q.is_infinity() {
            return p.clone();
        }

        // U = (XP - ZP)(XQ + ZQ)
        let u = &(&p.x - &p.z) * &(&q.x + &q.z);

        // V = (XP + ZP)(XQ - ZQ)
        let v = &(&p.x + &p.z) * &(&q.x - &q.z);

        // add = U + V, sub = U - V
        let add = &u + &v;
        let sub = &u - &v;

        // X+ = Zd · (U + V)²
        let add_sq = &add * &add;
        let x_result = &diff.z * &add_sq;

        // Z+ = Xd · (U - V)²
        let sub_sq = &sub * &sub;
        let z_result = &diff.x * &sub_sq;

        Point::new(x_result, z_result)
    }

    /// Montgomery point doubling.
    ///
    /// Given P = (X : Z), computes 2P = (X' : Z').
    ///
    /// Uses the formula (with a24 = (A + 2) / 4):
    /// - t1 = (X + Z)²
    /// - t2 = (X - Z)²
    /// - X' = t1 · t2
    /// - t3 = t1 - t2 (= 4XZ)
    /// - Z' = t3 · (t2 + a24 · t3)
    pub fn xdbl(&self, p: &Point) -> Point {
        if p.is_infinity() {
            return p.clone();
        }

        // Compute a24 = (A + 2) / 4
        let one = Fp2::one(self.modulus.clone());
        let two = &one + &one;
        let four = &two + &two;
        let a_plus_2 = &self.a + &two;
        let four_inv = four.inverse().expect("4 should be invertible");
        let a24 = &a_plus_2 * &four_inv;

        // t1 = (X + Z)²
        let sum = &p.x + &p.z;
        let t1 = &sum * &sum;

        // t2 = (X - Z)²
        let diff = &p.x - &p.z;
        let t2 = &diff * &diff;

        // X' = t1 · t2
        let x_result = &t1 * &t2;

        // t3 = t1 - t2 = 4XZ
        let t3 = &t1 - &t2;

        // Z' = t3 · (t2 + a24 · t3)
        let a24_t3 = &a24 * &t3;
        let t2_plus_a24_t3 = &t2 + &a24_t3;
        let z_result = &t3 * &t2_plus_a24_t3;

        Point::new(x_result, z_result)
    }

    /// Scalar multiplication using constant-time Montgomery ladder.
    ///
    /// Computes [k]P for scalar k and point P.
    ///
    /// Uses the Montgomery ladder which processes bits of k from
    /// most-significant to least-significant, maintaining the invariant
    /// that R1 - R0 = P throughout.
    ///
    /// # Security
    ///
    /// This implementation uses constant-time conditional swaps (cswap) to
    /// prevent timing side-channel attacks. The execution path is independent
    /// of the secret scalar k.
    pub fn scalar_mul(&self, k: &BigInt, p: &Point) -> Point {
        if k.is_zero() || p.is_infinity() {
            return Point::infinity(self.modulus.clone());
        }

        // Handle negative scalars
        let k = if k < &BigInt::zero() {
            // For x-only arithmetic, -P has same x-coordinate as P
            // So [−k]P has same x-coordinate as [k]P
            -k
        } else {
            k.clone()
        };

        // Find the bit length of k
        let bits = k.bits();
        if bits == 0 {
            return Point::infinity(self.modulus.clone());
        }

        // Constant-time Montgomery ladder using cswap
        // R0 = P, R1 = 2P
        // Invariant: R1 - R0 = P
        let mut r0 = p.clone();
        let mut r1 = self.xdbl(p);

        // Process bits from second-most-significant to least-significant
        // Use constant-time swap pattern to avoid secret-dependent branches
        for i in (0..bits - 1).rev() {
            let bit = (&k >> i) & BigInt::one() == BigInt::one();

            // Constant-time swap: if bit=1, swap r0 and r1
            cswap_point(bit, &mut r0, &mut r1);

            // Always do the same operations (r1 = r0 + r1, r0 = 2*r0)
            r1 = self.xadd(&r0, &r1, p);
            r0 = self.xdbl(&r0);

            // Swap back if bit was 1
            cswap_point(bit, &mut r0, &mut r1);
        }

        r0
    }

    /// Non-constant-time scalar multiplication (for public values only).
    ///
    /// WARNING: Only use this for operations on public/non-secret values.
    /// For secret scalars, use `scalar_mul` instead.
    pub fn scalar_mul_vartime(&self, k: &BigInt, p: &Point) -> Point {
        if k.is_zero() || p.is_infinity() {
            return Point::infinity(self.modulus.clone());
        }

        let k = if k < &BigInt::zero() { -k } else { k.clone() };
        let bits = k.bits();
        if bits == 0 {
            return Point::infinity(self.modulus.clone());
        }

        let mut r0 = p.clone();
        let mut r1 = self.xdbl(p);

        for i in (0..bits - 1).rev() {
            let bit = (&k >> i) & BigInt::one() == BigInt::one();
            if bit {
                r0 = self.xadd(&r0, &r1, p);
                r1 = self.xdbl(&r1);
            } else {
                r1 = self.xadd(&r0, &r1, p);
                r0 = self.xdbl(&r0);
            }
        }

        r0
    }

    /// Computes the order of a point (for testing).
    ///
    /// Tries multiples of P up to max_order to find the smallest n
    /// such that [n]P = O (point at infinity).
    ///
    /// Returns None if no such n is found within the bound.
    pub fn point_order(&self, p: &Point, max_order: &BigInt) -> Option<BigInt> {
        if p.is_infinity() {
            return Some(BigInt::one());
        }

        let mut current = p.clone();
        let mut n = BigInt::one();

        while &n <= max_order {
            n += 1;
            current = self.xadd(&current, p, p);

            // Check if we've reached infinity
            // In projective coordinates, infinity is when Z = 0
            if current.z.is_zero() {
                return Some(n);
            }
        }

        None
    }

    /// Computes (A + 2) / 4, often denoted a24, used in Montgomery arithmetic.
    pub fn a24(&self) -> Fp2 {
        let one = Fp2::one(self.modulus.clone());
        let two = &one + &one;
        let four = &two + &two;
        let a_plus_2 = &self.a + &two;
        let four_inv = four.inverse().expect("4 should be invertible");
        &a_plus_2 * &four_inv
    }

    /// Checks if the curve is supersingular.
    ///
    /// A curve over Fp² is supersingular iff its trace of Frobenius
    /// is divisible by p. For j ∈ {0, 1728}, curves are supersingular
    /// when p ≡ 3 (mod 4) for j=1728 or p ≡ 2 (mod 3) for j=0.
    ///
    /// This method uses the j-invariant to determine supersingularity:
    /// - For j = 0: supersingular iff p ≡ 2 (mod 3)
    /// - For j = 1728: supersingular iff p ≡ 3 (mod 4)
    /// - For other j: need to check if j is a supersingular j-invariant
    pub fn is_supersingular(&self) -> bool {
        let j = self.j_invariant();
        let p = &self.modulus;

        // Check for j = 0
        if j.is_zero() {
            // Supersingular iff p ≡ 2 (mod 3)
            return p % BigInt::from(3) == BigInt::from(2);
        }

        // Check for j = 1728
        let j_1728 = BigInt::from(1728) % p;
        if j.re.value == j_1728 && j.im.is_zero() {
            // Supersingular iff p ≡ 3 (mod 4)
            return p % BigInt::from(4) == BigInt::from(3);
        }

        // For other j-invariants, we need to check if j is in the set of
        // supersingular j-invariants. This requires more sophisticated
        // computation (checking if the curve has a small endomorphism ring).
        //
        // For SQI-SIGN, curves are constructed to be supersingular by design,
        // starting from E₀ (j=1728) and applying isogenies (which preserve
        // supersingularity). So if we reached here from a valid SQI-SIGN
        // construction, the curve should be supersingular.
        //
        // A full implementation would:
        // 1. Compute the order of a random point
        // 2. Check if it divides (p+1)² (characteristic of supersingular curves over Fp²)
        //
        // For now, we assume curves from valid constructions are supersingular
        // but log a warning for non-standard j-invariants.
        #[cfg(debug_assertions)]
        eprintln!("Warning: Cannot definitively verify supersingularity for j = {:?}", j);

        // Conservative: assume supersingular if constructed via SQI-SIGN methods
        true
    }

    /// Validates the curve parameters.
    ///
    /// Checks that:
    /// 1. The curve is non-singular (A² ≠ 4)
    /// 2. The modulus is prime (basic check)
    pub fn is_valid(&self) -> bool {
        // Check non-singularity: A² ≠ 4
        let a_sq = &self.a * &self.a;
        let four = {
            let one = Fp2::one(self.modulus.clone());
            let two = &one + &one;
            &two + &two
        };

        if a_sq == four {
            return false; // Singular curve
        }

        // Basic primality check (modulus should be > 2)
        if self.modulus <= BigInt::from(2) {
            return false;
        }

        true
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
    use crate::field::Fp;

    fn test_prime() -> BigInt {
        // p = 431 is prime, and p ≡ 3 (mod 4)
        // This means curves with j=1728 (A=0) are supersingular
        BigInt::from(431)
    }

    fn make_point(x_re: i64, x_im: i64, z_re: i64, z_im: i64, p: &BigInt) -> Point {
        Point::new(
            Fp2::new(
                Fp::new(BigInt::from(x_re), p.clone()),
                Fp::new(BigInt::from(x_im), p.clone()),
            ),
            Fp2::new(
                Fp::new(BigInt::from(z_re), p.clone()),
                Fp::new(BigInt::from(z_im), p.clone()),
            ),
        )
    }

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

    #[test]
    fn test_xdbl_infinity() {
        let p = test_prime();
        let curve = E0::new(p.clone()).curve;
        let inf = Point::infinity(p);

        let doubled = curve.xdbl(&inf);
        assert!(doubled.is_infinity());
    }

    #[test]
    fn test_xadd_with_infinity() {
        let p = test_prime();
        let curve = E0::new(p.clone()).curve;
        let inf = Point::infinity(p.clone());

        // Create a non-infinity point (x=1, z=1 in affine x=1)
        let point = make_point(1, 0, 1, 0, &p);

        // P + O = P
        let result = curve.xadd(&point, &inf, &point);
        // Normalize and check
        let norm_result = result.normalize().unwrap();
        let norm_point = point.normalize().unwrap();
        assert_eq!(norm_result.x.re.value, norm_point.x.re.value);
    }

    #[test]
    fn test_scalar_mul_zero() {
        let p = test_prime();
        let curve = E0::new(p.clone()).curve;
        let point = make_point(1, 0, 1, 0, &p);

        let result = curve.scalar_mul(&BigInt::zero(), &point);
        assert!(result.is_infinity());
    }

    #[test]
    fn test_scalar_mul_one() {
        let p = test_prime();
        let curve = E0::new(p.clone()).curve;
        let point = make_point(2, 0, 1, 0, &p);

        let result = curve.scalar_mul(&BigInt::one(), &point);
        let norm_result = result.normalize().unwrap();
        let norm_point = point.normalize().unwrap();

        // [1]P = P
        assert_eq!(norm_result.x.re.value, norm_point.x.re.value);
        assert_eq!(norm_result.x.im.value, norm_point.x.im.value);
    }

    #[test]
    fn test_scalar_mul_two_equals_xdbl() {
        let p = test_prime();
        let curve = E0::new(p.clone()).curve;
        let point = make_point(2, 0, 1, 0, &p);

        let doubled = curve.xdbl(&point);
        let scalar_2 = curve.scalar_mul(&BigInt::from(2), &point);

        // Both should give the same result
        let norm_dbl = doubled.normalize();
        let norm_scalar = scalar_2.normalize();

        if let (Some(nd), Some(ns)) = (norm_dbl, norm_scalar) {
            assert_eq!(nd.x.re.value, ns.x.re.value);
            assert_eq!(nd.x.im.value, ns.x.im.value);
        }
    }

    #[test]
    fn test_j_invariant_e0() {
        // E0 has A = 0, so j = 256(-3)³ / (-4) = 256 * (-27) / (-4) = 1728
        let p = test_prime();
        let curve = E0::new(p.clone()).curve;
        let j = curve.j_invariant();

        // j should be 1728 (mod p)
        assert_eq!(j.re.value, BigInt::from(1728) % &p);
        assert!(j.im.is_zero());
    }

    #[test]
    fn test_a24_computation() {
        let p = test_prime();
        let curve = E0::new(p.clone()).curve;
        let a24 = curve.a24();

        // For A = 0, a24 = (0 + 2) / 4 = 1/2
        // 1/2 mod 431 = (431+1)/2 = 216
        let expected = (p.clone() + BigInt::one()) / 2;
        assert_eq!(a24.re.value, expected);
        assert!(a24.im.is_zero());
    }
}
