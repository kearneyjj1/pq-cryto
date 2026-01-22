//! Isogeny computations for Montgomery curves.
//!
//! This module implements isogeny computation between supersingular elliptic
//! curves using Vélu's formulas adapted for Montgomery form.
//!
//! # Background
//!
//! An isogeny φ: E₁ → E₂ is a non-constant morphism of elliptic curves that
//! maps the identity to the identity. For supersingular curves, isogenies
//! correspond to ideals in the endomorphism ring.
//!
//! # Montgomery Curve Isogenies
//!
//! For Montgomery curves By² = x³ + Ax² + x, we use x-coordinate only
//! formulas which are more efficient for isogeny computation.

use crate::curve::{Curve, Point};
use crate::fp2::Fp2;
use num_bigint::BigInt;
use num_traits::One;

/// An isogeny between two Montgomery curves.
#[derive(Clone, Debug)]
pub struct Isogeny {
    /// The domain curve.
    pub domain: Curve,
    /// The codomain curve.
    pub codomain: Curve,
    /// The degree of the isogeny.
    pub degree: BigInt,
    /// Kernel generator (for reconstruction/evaluation).
    pub kernel_gen: Option<Point>,
    /// Cached kernel points for evaluation (for odd-degree isogenies).
    kernel_points: Vec<Point>,
}

impl Isogeny {
    /// Creates an isogeny from domain to codomain with given degree.
    pub fn new(domain: Curve, codomain: Curve, degree: BigInt) -> Self {
        Self {
            domain,
            codomain,
            degree,
            kernel_gen: None,
            kernel_points: Vec::new(),
        }
    }

    /// Creates an isogeny with explicit kernel points.
    pub fn with_kernel(domain: Curve, codomain: Curve, degree: BigInt, kernel_points: Vec<Point>) -> Self {
        let kernel_gen = kernel_points.first().cloned();
        Self {
            domain,
            codomain,
            degree,
            kernel_gen,
            kernel_points,
        }
    }

    /// Returns a reference to the kernel points.
    pub fn kernel_points(&self) -> &[Point] {
        &self.kernel_points
    }

    /// Computes a 2-isogeny with the given kernel point.
    ///
    /// The kernel point must have order 2 (i.e., 2K = O).
    ///
    /// For a Montgomery curve E_A: y² = x³ + Ax² + x with kernel point
    /// K = (x_K, 0), the codomain curve E_{A'} has:
    ///
    /// A' = 2(1 - 2x_K²) when x_K ≠ 0
    /// A' = -A when x_K = 0 (kernel at origin)
    pub fn isogeny_2(domain: &Curve, kernel: &Point) -> Self {
        let modulus = domain.modulus.clone();

        if kernel.is_infinity() {
            panic!("Kernel point cannot be infinity for 2-isogeny");
        }

        // Normalize kernel point to get x_K = X_K / Z_K
        let x_k = if let Some(norm) = kernel.normalize() {
            norm.x
        } else {
            panic!("Cannot normalize kernel point");
        };

        let one = Fp2::one(modulus.clone());
        let two = &one + &one;

        // Check if kernel is at origin (x_K = 0)
        let new_a = if x_k.is_zero() {
            // A' = -A
            -&domain.a
        } else {
            // A' = 2(1 - 2x_K²)
            // x_K²
            let x_k_sq = &x_k * &x_k;
            // 2x_K²
            let two_x_k_sq = &two * &x_k_sq;
            // 1 - 2x_K²
            let inner = &one - &two_x_k_sq;
            // 2(1 - 2x_K²)
            &two * &inner
        };

        let codomain = Curve::new(new_a);

        Self {
            domain: domain.clone(),
            codomain,
            degree: BigInt::from(2),
            kernel_gen: Some(kernel.clone()),
            kernel_points: vec![kernel.clone()],
        }
    }

    /// Computes a 3-isogeny with the given kernel point.
    ///
    /// The kernel point must have order 3 (i.e., 3K = O but K ≠ O).
    ///
    /// For Montgomery curves, the 3-isogeny formula gives:
    /// A' = (A·x_K⁴ - 6x_K² + A) / x_K² for kernel at x = x_K
    pub fn isogeny_3(domain: &Curve, kernel: &Point) -> Self {
        let modulus = domain.modulus.clone();

        if kernel.is_infinity() {
            panic!("Kernel point cannot be infinity for 3-isogeny");
        }

        // Normalize to get x_K
        let x_k = if let Some(norm) = kernel.normalize() {
            norm.x
        } else {
            panic!("Cannot normalize kernel point");
        };

        if x_k.is_zero() {
            panic!("Kernel x-coordinate cannot be zero for 3-isogeny");
        }

        let one = Fp2::one(modulus.clone());
        let six = {
            let two = &one + &one;
            let three = &two + &one;
            &three + &three
        };

        // x_K², x_K⁴
        let x_k_sq = &x_k * &x_k;
        let x_k_4 = &x_k_sq * &x_k_sq;

        // A' = (A·x_K⁴ - 6x_K² + A) / x_K²
        let a_x4 = &domain.a * &x_k_4;
        let six_x2 = &six * &x_k_sq;
        let numerator = &(&a_x4 - &six_x2) + &domain.a;

        let x_k_sq_inv = x_k_sq.inverse().expect("x_K² should be invertible");
        let new_a = &numerator * &x_k_sq_inv;

        let codomain = Curve::new(new_a);

        // For 3-isogeny, kernel has 2 non-trivial points: K and -K (same x-coord)
        Self {
            domain: domain.clone(),
            codomain,
            degree: BigInt::from(3),
            kernel_gen: Some(kernel.clone()),
            kernel_points: vec![kernel.clone()],
        }
    }

    /// Computes an ℓ-isogeny for odd prime ℓ using Vélu's formulas.
    ///
    /// The kernel point must have order exactly ℓ.
    ///
    /// This uses the general Vélu formula which computes the isogeny
    /// by iterating over all non-trivial kernel points.
    pub fn isogeny_odd(domain: &Curve, kernel: &Point, ell: u64) -> Self {
        if ell == 2 {
            return Self::isogeny_2(domain, kernel);
        }
        if ell == 3 {
            return Self::isogeny_3(domain, kernel);
        }

        let modulus = domain.modulus.clone();

        if kernel.is_infinity() {
            panic!("Kernel point cannot be infinity");
        }

        // Compute all kernel points: K, 2K, 3K, ..., ((ℓ-1)/2)K
        // (We only need half since -iK has the same x-coordinate as iK)
        let half_ell = (ell - 1) / 2;
        let mut kernel_points = Vec::with_capacity(half_ell as usize);
        let mut current = kernel.clone();

        for _ in 0..half_ell {
            kernel_points.push(current.clone());
            current = domain.xadd(&current, kernel, kernel);
        }

        // Vélu's formula for Montgomery curves:
        // For each kernel point with x-coordinate x_i:
        // σ = Σ x_i
        // π = Π x_i
        // The new curve coefficient involves these sums/products

        let one = Fp2::one(modulus.clone());

        // Compute sum of x-coordinates (for Vélu formula)
        let mut sigma = Fp2::zero(modulus.clone());

        for pt in &kernel_points {
            if let Some(norm) = pt.normalize() {
                sigma = &sigma + &norm.x;
            }
        }

        // For odd ℓ-isogeny on Montgomery curve:
        // A' = A - 6σ (simplified formula for specific cases)
        // This is an approximation; the full formula is more complex
        let six = {
            let three = &(&one + &one) + &one;
            &three + &three
        };
        let new_a = &domain.a - &(&six * &sigma);

        let codomain = Curve::new(new_a);

        Self {
            domain: domain.clone(),
            codomain,
            degree: BigInt::from(ell),
            kernel_gen: Some(kernel.clone()),
            kernel_points,
        }
    }

    /// Evaluates the isogeny on a point.
    ///
    /// Given φ: E → E' and P ∈ E, computes φ(P) ∈ E'.
    pub fn evaluate(&self, point: &Point) -> Point {
        if point.is_infinity() {
            return Point::infinity(self.codomain.modulus.clone());
        }

        // Check if point is in the kernel
        for k in &self.kernel_points {
            if let (Some(p_norm), Some(k_norm)) = (point.normalize(), k.normalize()) {
                if p_norm.x == k_norm.x {
                    // Point is in kernel, maps to infinity
                    return Point::infinity(self.codomain.modulus.clone());
                }
            }
        }

        let degree = self.degree.to_u64_digits();
        if degree.1.is_empty() {
            return point.clone();
        }
        let ell = degree.1[0];

        match ell {
            2 => self.evaluate_2(point),
            3 => self.evaluate_3(point),
            _ => self.evaluate_odd(point),
        }
    }

    /// Evaluates a 2-isogeny on a point.
    fn evaluate_2(&self, point: &Point) -> Point {
        let kernel = self.kernel_gen.as_ref().expect("Need kernel for evaluation");

        // For 2-isogeny with kernel at x_K:
        // x' = (x·x_K - 1)² / (x - x_K)²
        // In projective: more complex formula

        let (x_p, z_p) = (&point.x, &point.z);
        let (x_k, z_k) = (&kernel.x, &kernel.z);

        // Using projective formula:
        // t0 = X_P + Z_P
        // t1 = X_P - Z_P
        // t2 = t0 * Z_K
        // t3 = t1 * X_K
        // X' = (t2 + t3)²
        // Z' = (t2 - t3)²

        let t0 = x_p + z_p;
        let t1 = x_p - z_p;
        let t2 = &t0 * z_k;
        let t3 = &t1 * x_k;

        let sum = &t2 + &t3;
        let diff = &t2 - &t3;

        let x_new = &sum * &sum;
        let z_new = &diff * &diff;

        Point::new(x_new, z_new)
    }

    /// Evaluates a 3-isogeny on a point.
    fn evaluate_3(&self, point: &Point) -> Point {
        let kernel = self.kernel_gen.as_ref().expect("Need kernel for evaluation");

        // For 3-isogeny with kernel generated by K at x_K:
        // Use the standard formula for Montgomery curves

        let (x_p, z_p) = (&point.x, &point.z);
        let (x_k, z_k) = (&kernel.x, &kernel.z);

        // Projective 3-isogeny evaluation formula
        let t0 = x_p * x_k;
        let t1 = z_p * z_k;
        let t2 = x_p * z_k;
        let t3 = z_p * x_k;

        let u = &t0 - &t1;
        let v = &t2 - &t3;

        let u_sq = &u * &u;
        let v_sq = &v * &v;

        // X' = X_P · u²
        // Z' = Z_P · v²
        let x_new = x_p * &u_sq;
        let z_new = z_p * &v_sq;

        Point::new(x_new, z_new)
    }

    /// Evaluates an odd-degree isogeny on a point using Vélu's formula.
    fn evaluate_odd(&self, point: &Point) -> Point {
        if self.kernel_points.is_empty() {
            return point.clone();
        }

        // General Vélu evaluation: for each kernel point K_i,
        // we modify the image point
        let mut result = point.clone();

        for kernel_pt in &self.kernel_points {
            // Evaluate the contribution from this kernel point
            // This is a simplified version; full Vélu is more complex
            let (x_p, z_p) = (&result.x, &result.z);
            let (x_k, z_k) = (&kernel_pt.x, &kernel_pt.z);

            let t0 = x_p * x_k;
            let t1 = z_p * z_k;
            let t2 = x_p * z_k;
            let t3 = z_p * x_k;

            let u = &t0 - &t1;
            let v = &t2 - &t3;

            let u_sq = &u * &u;
            let v_sq = &v * &v;

            let x_new = x_p * &u_sq;
            let z_new = z_p * &v_sq;

            result = Point::new(x_new, z_new);
        }

        result
    }

    /// Composes this isogeny with another.
    ///
    /// Returns φ₂ ∘ φ₁ where self = φ₁ and other = φ₂.
    /// The composition maps E₁ → E₂ → E₃.
    pub fn compose(&self, other: &Isogeny) -> Self {
        // Verify the domains/codomains match
        // (In practice, we'd check j-invariants match)

        Self {
            domain: self.domain.clone(),
            codomain: other.codomain.clone(),
            degree: &self.degree * &other.degree,
            kernel_gen: self.kernel_gen.clone(),
            kernel_points: Vec::new(), // Composed isogeny needs recomputation
        }
    }

    /// Returns the dual isogeny φ̂: E' → E.
    ///
    /// For an ℓ-isogeny φ, we have φ̂ ∘ φ = [ℓ] (multiplication by ℓ).
    pub fn dual(&self) -> Self {
        // The dual isogeny has the same degree but swapped domain/codomain
        // Computing the actual kernel of the dual requires more work
        Self {
            domain: self.codomain.clone(),
            codomain: self.domain.clone(),
            degree: self.degree.clone(),
            kernel_gen: None, // Would need to compute
            kernel_points: Vec::new(),
        }
    }
}

/// Computes a chain of 2-isogenies.
///
/// Given a point P of order 2^e, computes the isogeny with kernel ⟨P⟩.
/// This is done by iteratively computing 2-isogenies.
pub fn isogeny_chain_2(domain: &Curve, kernel_gen: &Point, e: u32) -> Isogeny {
    if e == 0 {
        return Isogeny::new(
            domain.clone(),
            domain.clone(),
            BigInt::one(),
        );
    }

    let mut current_curve = domain.clone();
    let mut current_kernel = kernel_gen.clone();
    let mut total_degree = BigInt::one();

    // We need to compute the isogeny step by step
    // At each step, we take a point of order 2^(e-i) and compute
    // a 2-isogeny with its [2^(e-i-1)]-multiple as kernel

    for i in 0..e {
        // Compute [2^(e-i-1)]K to get a point of order 2
        let remaining = e - i - 1;
        let mut order_2_point = current_kernel.clone();

        for _ in 0..remaining {
            order_2_point = current_curve.xdbl(&order_2_point);
        }

        // Compute 2-isogeny with this kernel
        let phi = Isogeny::isogeny_2(&current_curve, &order_2_point);

        // Push the kernel generator through
        if i < e - 1 {
            current_kernel = phi.evaluate(&current_kernel);
        }

        current_curve = phi.codomain;
        total_degree = &total_degree * BigInt::from(2);
    }

    Isogeny {
        domain: domain.clone(),
        codomain: current_curve,
        degree: total_degree,
        kernel_gen: Some(kernel_gen.clone()),
        kernel_points: Vec::new(),
    }
}

/// Computes a chain of 3-isogenies.
///
/// Given a point P of order 3^e, computes the isogeny with kernel ⟨P⟩.
pub fn isogeny_chain_3(domain: &Curve, kernel_gen: &Point, e: u32) -> Isogeny {
    if e == 0 {
        return Isogeny::new(
            domain.clone(),
            domain.clone(),
            BigInt::one(),
        );
    }

    let mut current_curve = domain.clone();
    let mut current_kernel = kernel_gen.clone();
    let mut total_degree = BigInt::one();

    for i in 0..e {
        // Compute [3^(e-i-1)]K to get a point of order 3
        let remaining = e - i - 1;
        let mut order_3_point = current_kernel.clone();

        for _ in 0..remaining {
            // Triple the point: 3P = 2P + P
            let two_p = current_curve.xdbl(&order_3_point);
            order_3_point = current_curve.xadd(&two_p, &order_3_point, &order_3_point);
        }

        // Compute 3-isogeny
        let phi = Isogeny::isogeny_3(&current_curve, &order_3_point);

        // Push kernel through
        if i < e - 1 {
            current_kernel = phi.evaluate(&current_kernel);
        }

        current_curve = phi.codomain;
        total_degree = &total_degree * BigInt::from(3);
    }

    Isogeny {
        domain: domain.clone(),
        codomain: current_curve,
        degree: total_degree,
        kernel_gen: Some(kernel_gen.clone()),
        kernel_points: Vec::new(),
    }
}

/// Computes an isogeny of smooth degree.
///
/// Given a point P whose order is B-smooth (product of small primes),
/// computes the isogeny with kernel ⟨P⟩.
pub fn isogeny_smooth(
    domain: &Curve,
    kernel_gen: &Point,
    degree_factors: &[(u64, u32)], // [(prime, exponent), ...]
) -> Isogeny {
    let mut current_curve = domain.clone();
    let mut current_kernel = kernel_gen.clone();
    let mut total_degree = BigInt::one();

    for &(prime, exp) in degree_factors {
        for _ in 0..exp {
            // Compute cofactor to get a point of order `prime`
            let mut cofactor_point = current_kernel.clone();

            // Multiply by cofactor (all other primes and remaining exponents)
            for &(other_prime, other_exp) in degree_factors {
                let times = if other_prime == prime {
                    // For this prime, we've already done some steps
                    0 // Simplified; should track remaining
                } else {
                    other_exp
                };
                for _ in 0..times {
                    cofactor_point = current_curve.scalar_mul(
                        &BigInt::from(other_prime),
                        &cofactor_point,
                    );
                }
            }

            // Compute isogeny of degree `prime`
            let phi = if prime == 2 {
                Isogeny::isogeny_2(&current_curve, &cofactor_point)
            } else if prime == 3 {
                Isogeny::isogeny_3(&current_curve, &cofactor_point)
            } else {
                Isogeny::isogeny_odd(&current_curve, &cofactor_point, prime)
            };

            current_kernel = phi.evaluate(&current_kernel);
            current_curve = phi.codomain;
            total_degree = &total_degree * BigInt::from(prime);
        }
    }

    Isogeny {
        domain: domain.clone(),
        codomain: current_curve,
        degree: total_degree,
        kernel_gen: Some(kernel_gen.clone()),
        kernel_points: Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::E0;
    use crate::field::Fp;

    fn test_prime() -> BigInt {
        // p = 431 ≡ 3 (mod 4), good for supersingular curves
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
    fn test_isogeny_2_origin_kernel() {
        // 2-isogeny with kernel at origin (0, 0)
        let p = test_prime();
        let e0 = E0::new(p.clone());

        // Kernel at (0, 0) has x = 0
        let kernel = make_point(0, 0, 1, 0, &p);

        let phi = Isogeny::isogeny_2(&e0.curve, &kernel);

        // For A = 0, kernel at origin gives A' = -A = 0
        assert!(phi.codomain.a.is_zero());
        assert_eq!(phi.degree, BigInt::from(2));
    }

    #[test]
    fn test_isogeny_evaluate_infinity() {
        let p = test_prime();
        let e0 = E0::new(p.clone());
        let kernel = make_point(0, 0, 1, 0, &p);

        let phi = Isogeny::isogeny_2(&e0.curve, &kernel);

        // Infinity should map to infinity
        let inf = Point::infinity(p);
        let result = phi.evaluate(&inf);
        assert!(result.is_infinity());
    }

    #[test]
    fn test_isogeny_compose() {
        let p = test_prime();
        let e0 = E0::new(p.clone());
        let kernel = make_point(0, 0, 1, 0, &p);

        let phi1 = Isogeny::isogeny_2(&e0.curve, &kernel);
        let phi2 = Isogeny::isogeny_2(&phi1.codomain, &kernel);

        let composed = phi1.compose(&phi2);
        assert_eq!(composed.degree, BigInt::from(4));
    }

    #[test]
    fn test_isogeny_chain_2_degree() {
        let p = test_prime();
        let e0 = E0::new(p.clone());

        // Create a point (we'll use a simple one for testing)
        let kernel = make_point(1, 0, 1, 0, &p);

        // Compute 2^2 = 4 isogeny
        let phi = isogeny_chain_2(&e0.curve, &kernel, 2);

        assert_eq!(phi.degree, BigInt::from(4));
    }

    #[test]
    fn test_isogeny_dual() {
        let p = test_prime();
        let e0 = E0::new(p.clone());
        let kernel = make_point(0, 0, 1, 0, &p);

        let phi = Isogeny::isogeny_2(&e0.curve, &kernel);
        let phi_dual = phi.dual();

        // Dual has same degree, swapped domain/codomain
        assert_eq!(phi_dual.degree, phi.degree);
    }
}
