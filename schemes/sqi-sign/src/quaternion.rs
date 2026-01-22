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
use num_traits::{One, Zero};
use std::ops::{Add, Mul, Neg, Sub};

/// Computes the greatest common divisor of two integers.
fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    let mut a = if a < &BigInt::zero() { -a } else { a.clone() };
    let mut b = if b < &BigInt::zero() { -b } else { b.clone() };

    while !b.is_zero() {
        let t = b.clone();
        b = &a % &b;
        a = t;
    }
    a
}

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
        if den.is_zero() {
            panic!("Rational denominator cannot be zero");
        }

        // Ensure denominator is positive
        let (num, den) = if den < BigInt::zero() {
            (-num, -den)
        } else {
            (num, den)
        };

        // Reduce to lowest terms using GCD
        let g = gcd(&num, &den);
        if g == BigInt::one() {
            Self { num, den }
        } else {
            Self {
                num: &num / &g,
                den: &den / &g,
            }
        }
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

    /// Returns one.
    pub fn one() -> Self {
        Self {
            num: BigInt::one(),
            den: BigInt::one(),
        }
    }

    /// Computes the multiplicative inverse.
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        Some(Self::new(self.den.clone(), self.num.clone()))
    }
}

impl Add for &Rational {
    type Output = Rational;

    fn add(self, other: &Rational) -> Rational {
        // a/b + c/d = (ad + bc) / bd
        let num = &self.num * &other.den + &self.den * &other.num;
        let den = &self.den * &other.den;
        Rational::new(num, den)
    }
}

impl Sub for &Rational {
    type Output = Rational;

    fn sub(self, other: &Rational) -> Rational {
        // a/b - c/d = (ad - bc) / bd
        let num = &self.num * &other.den - &self.den * &other.num;
        let den = &self.den * &other.den;
        Rational::new(num, den)
    }
}

impl Mul for &Rational {
    type Output = Rational;

    fn mul(self, other: &Rational) -> Rational {
        // (a/b) * (c/d) = (ac) / (bd)
        Rational::new(&self.num * &other.num, &self.den * &other.den)
    }
}

impl Neg for &Rational {
    type Output = Rational;

    fn neg(self) -> Rational {
        Rational::new(-&self.num, self.den.clone())
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
        Self {
            a: self.a.clone(),
            b: -&self.b,
            c: -&self.c,
            d: -&self.d,
            p: self.p.clone(),
        }
    }

    /// Computes the reduced norm: nrd(α) = α * conj(α).
    ///
    /// For α = a + bi + cj + dk in B_{p,∞}:
    /// nrd(α) = a² + b² + p·c² + p·d²
    pub fn reduced_norm(&self) -> Rational {
        let p_rat = Rational::from_int(self.p.clone());

        // a²
        let a_sq = &self.a * &self.a;
        // b²
        let b_sq = &self.b * &self.b;
        // c²
        let c_sq = &self.c * &self.c;
        // d²
        let d_sq = &self.d * &self.d;

        // p·c² + p·d²
        let p_c_sq = &p_rat * &c_sq;
        let p_d_sq = &p_rat * &d_sq;

        // a² + b² + p·c² + p·d²
        let sum1 = &a_sq + &b_sq;
        let sum2 = &p_c_sq + &p_d_sq;
        &sum1 + &sum2
    }

    /// Computes the reduced trace: trd(α) = α + conj(α) = 2a.
    pub fn reduced_trace(&self) -> Rational {
        let two = Rational::from_int(BigInt::from(2));
        &two * &self.a
    }

    /// Multiplies two quaternions.
    ///
    /// For B_{p,∞} with i² = -1, j² = -p, k = ij:
    /// - ij = k, ji = -k
    /// - jk = pi, kj = -pi
    /// - ik = -j, ki = j
    ///
    /// (a₁ + b₁i + c₁j + d₁k)(a₂ + b₂i + c₂j + d₂k) yields:
    /// - Coeff of 1: a₁a₂ - b₁b₂ - p·c₁c₂ - p·d₁d₂
    /// - Coeff of i: a₁b₂ + b₁a₂ + p·c₁d₂ - p·d₁c₂
    /// - Coeff of j: a₁c₂ - b₁d₂ + c₁a₂ + d₁b₂
    /// - Coeff of k: a₁d₂ + b₁c₂ - c₁b₂ + d₁a₂
    pub fn mul(&self, other: &Quaternion) -> Quaternion {
        debug_assert_eq!(self.p, other.p, "Quaternions must have same p");

        let p_rat = Rational::from_int(self.p.clone());

        // Coefficient of 1: a₁a₂ - b₁b₂ - p·c₁c₂ - p·d₁d₂
        let a1a2 = &self.a * &other.a;
        let b1b2 = &self.b * &other.b;
        let c1c2 = &self.c * &other.c;
        let d1d2 = &self.d * &other.d;
        let new_a = &(&a1a2 - &b1b2) - &(&p_rat * &(&c1c2 + &d1d2));

        // Coefficient of i: a₁b₂ + b₁a₂ + p·c₁d₂ - p·d₁c₂
        let a1b2 = &self.a * &other.b;
        let b1a2 = &self.b * &other.a;
        let c1d2 = &self.c * &other.d;
        let d1c2 = &self.d * &other.c;
        let new_b = &(&a1b2 + &b1a2) + &(&p_rat * &(&c1d2 - &d1c2));

        // Coefficient of j: a₁c₂ - b₁d₂ + c₁a₂ + d₁b₂
        let a1c2 = &self.a * &other.c;
        let b1d2 = &self.b * &other.d;
        let c1a2 = &self.c * &other.a;
        let d1b2 = &self.d * &other.b;
        let new_c = &(&a1c2 - &b1d2) + &(&c1a2 + &d1b2);

        // Coefficient of k: a₁d₂ + b₁c₂ - c₁b₂ + d₁a₂
        let a1d2 = &self.a * &other.d;
        let b1c2 = &self.b * &other.c;
        let c1b2 = &self.c * &other.b;
        let d1a2 = &self.d * &other.a;
        let new_d = &(&a1d2 + &b1c2) - &(&c1b2 - &d1a2);

        Quaternion {
            a: new_a,
            b: new_b,
            c: new_c,
            d: new_d,
            p: self.p.clone(),
        }
    }

    /// Computes the inverse (if it exists).
    ///
    /// α⁻¹ = conj(α) / nrd(α)
    pub fn inverse(&self) -> Option<Quaternion> {
        let norm = self.reduced_norm();
        if norm.is_zero() {
            return None;
        }

        let norm_inv = norm.inverse()?;
        let conj = self.conjugate();

        // Multiply each coefficient by 1/nrd(α)
        Some(Quaternion {
            a: &conj.a * &norm_inv,
            b: &conj.b * &norm_inv,
            c: &conj.c * &norm_inv,
            d: &conj.d * &norm_inv,
            p: self.p.clone(),
        })
    }

    /// Adds two quaternions.
    pub fn add(&self, other: &Quaternion) -> Quaternion {
        debug_assert_eq!(self.p, other.p, "Quaternions must have same p");
        Quaternion {
            a: &self.a + &other.a,
            b: &self.b + &other.b,
            c: &self.c + &other.c,
            d: &self.d + &other.d,
            p: self.p.clone(),
        }
    }

    /// Subtracts two quaternions.
    pub fn sub(&self, other: &Quaternion) -> Quaternion {
        debug_assert_eq!(self.p, other.p, "Quaternions must have same p");
        Quaternion {
            a: &self.a - &other.a,
            b: &self.b - &other.b,
            c: &self.c - &other.c,
            d: &self.d - &other.d,
            p: self.p.clone(),
        }
    }

    /// Negates a quaternion.
    pub fn neg(&self) -> Quaternion {
        Quaternion {
            a: -&self.a,
            b: -&self.b,
            c: -&self.c,
            d: -&self.d,
            p: self.p.clone(),
        }
    }

    /// Checks if this quaternion is zero.
    pub fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero() && self.c.is_zero() && self.d.is_zero()
    }

    /// Creates a quaternion from just the real part.
    pub fn from_rational(a: Rational, p: BigInt) -> Self {
        Self {
            a,
            b: Rational::zero(),
            c: Rational::zero(),
            d: Rational::zero(),
            p,
        }
    }

    /// Creates the basis element i.
    pub fn i(p: BigInt) -> Self {
        Self {
            a: Rational::zero(),
            b: Rational::one(),
            c: Rational::zero(),
            d: Rational::zero(),
            p,
        }
    }

    /// Creates the basis element j.
    pub fn j(p: BigInt) -> Self {
        Self {
            a: Rational::zero(),
            b: Rational::zero(),
            c: Rational::one(),
            d: Rational::zero(),
            p,
        }
    }

    /// Creates the basis element k.
    pub fn k(p: BigInt) -> Self {
        Self {
            a: Rational::zero(),
            b: Rational::zero(),
            c: Rational::zero(),
            d: Rational::one(),
            p,
        }
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
    /// For p ≡ 3 (mod 4), the standard maximal order is:
    /// O₀ = Z⟨1, i, (1+j)/2, (i+k)/2⟩
    ///
    /// This is the order corresponding to End(E₀) where E₀: y² = x³ + x.
    pub fn standard(p: BigInt) -> Self {
        // Basis element 1
        let b1 = Quaternion::one(p.clone());

        // Basis element i
        let b2 = Quaternion::i(p.clone());

        // Basis element (1+j)/2
        let one_plus_j = Quaternion::new(
            Rational::one(),
            Rational::zero(),
            Rational::one(),
            Rational::zero(),
            p.clone(),
        );
        let half = Rational::new(BigInt::one(), BigInt::from(2));
        let b3 = Quaternion::new(
            &half * &one_plus_j.a,
            &half * &one_plus_j.b,
            &half * &one_plus_j.c,
            &half * &one_plus_j.d,
            p.clone(),
        );

        // Basis element (i+k)/2
        let i_plus_k = Quaternion::new(
            Rational::zero(),
            Rational::one(),
            Rational::zero(),
            Rational::one(),
            p.clone(),
        );
        let b4 = Quaternion::new(
            &half * &i_plus_k.a,
            &half * &i_plus_k.b,
            &half * &i_plus_k.c,
            &half * &i_plus_k.d,
            p.clone(),
        );

        Self {
            basis: [b1, b2, b3, b4],
            p,
        }
    }

    /// Tests if a quaternion is in this order.
    ///
    /// A quaternion q is in the order if it can be written as an integer
    /// linear combination of the basis elements.
    pub fn contains(&self, q: &Quaternion) -> bool {
        // For the standard order O₀ = Z⟨1, i, (1+j)/2, (i+k)/2⟩,
        // q = a + bi + cj + dk is in O₀ iff:
        // - a, b ∈ Z or a, b ∈ Z + 1/2 (with matching parity)
        // - c, d ∈ Z or c, d ∈ Z + 1/2 (with matching parity)
        // More precisely: 2a, 2b, 2c, 2d ∈ Z and a+c, b+d ∈ Z

        // Check if 2a, 2b, 2c, 2d are integers
        let two = BigInt::from(2);
        let two_a = &q.a.num * &two;
        let two_b = &q.b.num * &two;
        let two_c = &q.c.num * &two;
        let two_d = &q.d.num * &two;

        // Check divisibility by denominators
        if &two_a % &q.a.den != BigInt::zero() {
            return false;
        }
        if &two_b % &q.b.den != BigInt::zero() {
            return false;
        }
        if &two_c % &q.c.den != BigInt::zero() {
            return false;
        }
        if &two_d % &q.d.den != BigInt::zero() {
            return false;
        }

        // Check that a+c and b+d are integers
        let a_plus_c = &q.a + &q.c;
        let b_plus_d = &q.b + &q.d;

        a_plus_c.den == BigInt::one() && b_plus_d.den == BigInt::one()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_prime() -> BigInt {
        // p = 431 ≡ 3 (mod 4)
        BigInt::from(431)
    }

    #[test]
    fn test_rational_zero() {
        let z = Rational::zero();
        assert!(z.is_zero());
    }

    #[test]
    fn test_rational_arithmetic() {
        let a = Rational::new(BigInt::from(1), BigInt::from(2)); // 1/2
        let b = Rational::new(BigInt::from(1), BigInt::from(3)); // 1/3

        // 1/2 + 1/3 = 5/6
        let sum = &a + &b;
        assert_eq!(sum.num, BigInt::from(5));
        assert_eq!(sum.den, BigInt::from(6));

        // 1/2 * 1/3 = 1/6
        let prod = &a * &b;
        assert_eq!(prod.num, BigInt::from(1));
        assert_eq!(prod.den, BigInt::from(6));
    }

    #[test]
    fn test_rational_gcd_reduction() {
        // 4/6 should reduce to 2/3
        let r = Rational::new(BigInt::from(4), BigInt::from(6));
        assert_eq!(r.num, BigInt::from(2));
        assert_eq!(r.den, BigInt::from(3));
    }

    #[test]
    fn test_quaternion_creation() {
        let p = BigInt::from(101);
        let q = Quaternion::one(p.clone());
        assert!(!q.a.is_zero());
        assert!(q.b.is_zero());
    }

    #[test]
    fn test_quaternion_conjugate() {
        let p = test_prime();
        let q = Quaternion::new(
            Rational::from_int(BigInt::from(1)),
            Rational::from_int(BigInt::from(2)),
            Rational::from_int(BigInt::from(3)),
            Rational::from_int(BigInt::from(4)),
            p,
        );

        let conj = q.conjugate();
        assert_eq!(conj.a.num, BigInt::from(1));
        assert_eq!(conj.b.num, BigInt::from(-2));
        assert_eq!(conj.c.num, BigInt::from(-3));
        assert_eq!(conj.d.num, BigInt::from(-4));
    }

    #[test]
    fn test_quaternion_reduced_norm() {
        let p = test_prime();

        // For q = 1 + 2i + 3j + 4k with p = 431:
        // nrd(q) = 1² + 2² + 431·3² + 431·4² = 1 + 4 + 3879 + 6896 = 10780
        let q = Quaternion::new(
            Rational::from_int(BigInt::from(1)),
            Rational::from_int(BigInt::from(2)),
            Rational::from_int(BigInt::from(3)),
            Rational::from_int(BigInt::from(4)),
            p.clone(),
        );

        let norm = q.reduced_norm();
        let expected = 1 + 4 + 431 * 9 + 431 * 16;
        assert_eq!(norm.num, BigInt::from(expected));
        assert_eq!(norm.den, BigInt::one());
    }

    #[test]
    fn test_quaternion_reduced_trace() {
        let p = test_prime();
        let q = Quaternion::new(
            Rational::from_int(BigInt::from(5)),
            Rational::from_int(BigInt::from(2)),
            Rational::from_int(BigInt::from(3)),
            Rational::from_int(BigInt::from(4)),
            p,
        );

        let trace = q.reduced_trace();
        // trd(q) = 2a = 10
        assert_eq!(trace.num, BigInt::from(10));
        assert_eq!(trace.den, BigInt::one());
    }

    #[test]
    fn test_quaternion_i_squared() {
        let p = test_prime();
        let i = Quaternion::i(p);

        // i² = -1
        let i_sq = i.mul(&i);
        assert_eq!(i_sq.a.num, BigInt::from(-1));
        assert!(i_sq.b.is_zero());
        assert!(i_sq.c.is_zero());
        assert!(i_sq.d.is_zero());
    }

    #[test]
    fn test_quaternion_j_squared() {
        let p = test_prime();
        let j = Quaternion::j(p.clone());

        // j² = -p
        let j_sq = j.mul(&j);
        assert_eq!(j_sq.a.num, -p);
        assert!(j_sq.b.is_zero());
        assert!(j_sq.c.is_zero());
        assert!(j_sq.d.is_zero());
    }

    #[test]
    fn test_quaternion_ij_equals_k() {
        let p = test_prime();
        let i = Quaternion::i(p.clone());
        let j = Quaternion::j(p.clone());
        let _k = Quaternion::k(p); // Reference for comparison

        // ij = k
        let ij = i.mul(&j);
        assert!(ij.a.is_zero());
        assert!(ij.b.is_zero());
        assert!(ij.c.is_zero());
        assert_eq!(ij.d.num, BigInt::one());
    }

    #[test]
    fn test_quaternion_ji_equals_neg_k() {
        let p = test_prime();
        let i = Quaternion::i(p.clone());
        let j = Quaternion::j(p);

        // ji = -k
        let ji = j.mul(&i);
        assert!(ji.a.is_zero());
        assert!(ji.b.is_zero());
        assert!(ji.c.is_zero());
        assert_eq!(ji.d.num, BigInt::from(-1));
    }

    #[test]
    fn test_quaternion_inverse() {
        let p = test_prime();
        let q = Quaternion::new(
            Rational::from_int(BigInt::from(1)),
            Rational::from_int(BigInt::from(1)),
            Rational::from_int(BigInt::from(0)),
            Rational::from_int(BigInt::from(0)),
            p.clone(),
        );

        // q = 1 + i, nrd(q) = 1 + 1 = 2
        // q^(-1) = (1 - i) / 2 = 1/2 - i/2
        let q_inv = q.inverse().unwrap();

        // Check q * q^(-1) = 1
        let product = q.mul(&q_inv);
        assert_eq!(product.a.num, BigInt::one());
        assert_eq!(product.a.den, BigInt::one());
        assert!(product.b.is_zero());
        assert!(product.c.is_zero());
        assert!(product.d.is_zero());
    }

    #[test]
    fn test_maximal_order_standard() {
        let p = test_prime();
        let order = MaximalOrder::standard(p.clone());

        // Check that basis elements are in the order
        for basis_elem in &order.basis {
            assert!(order.contains(basis_elem));
        }
    }

    #[test]
    fn test_maximal_order_contains_integers() {
        let p = test_prime();
        let order = MaximalOrder::standard(p.clone());

        // Integer quaternions should be in the order
        let q = Quaternion::new(
            Rational::from_int(BigInt::from(5)),
            Rational::from_int(BigInt::from(3)),
            Rational::from_int(BigInt::from(0)),
            Rational::from_int(BigInt::from(0)),
            p,
        );
        assert!(order.contains(&q));
    }

    #[test]
    fn test_maximal_order_half_basis() {
        let p = test_prime();
        let order = MaximalOrder::standard(p.clone());

        // (1+j)/2 should be in the order
        let half_one_j = Quaternion::new(
            Rational::new(BigInt::from(1), BigInt::from(2)),
            Rational::zero(),
            Rational::new(BigInt::from(1), BigInt::from(2)),
            Rational::zero(),
            p,
        );
        assert!(order.contains(&half_one_j));
    }
}
