//! Complex floating-point FFT for FALCON.
//!
//! This module implements the Fast Fourier Transform over complex numbers,
//! specifically for the negacyclic polynomial ring Z[X]/(X^n + 1).
//!
//! FALCON uses floating-point FFT (not integer NTT) for:
//! - Efficient polynomial multiplication in key generation
//! - The Fast Fourier Sampling algorithm during signing
//!
//! ## Ordering Convention
//!
//! This implementation uses the Falcon tree-ordered (interleaved) FFT
//! matching the reference implementation (tprest/falcon.py, `fft.py`).
//! The root ordering satisfies two key properties:
//!
//! 1. `w[2i+1] = -w[2i]` (adjacent pairs are negations)
//! 2. `w[i+n/2] = conj(w[i])` (second half conjugates the first half)
//!
//! Property 2 ensures that `split_fft` of a Hermitian-symmetric vector
//! produces Hermitian-symmetric children, preserving real-valuedness
//! down to the n=1 leaves. This is essential for ffSampling correctness.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use zeroize::Zeroize;

/// A complex number with f64 components.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Zeroize for Complex {
    fn zeroize(&mut self) {
        self.re.zeroize();
        self.im.zeroize();
    }
}

impl Complex {
    /// Zero complex number.
    pub const ZERO: Complex = Complex { re: 0.0, im: 0.0 };

    /// One (real unit).
    pub const ONE: Complex = Complex { re: 1.0, im: 0.0 };

    /// Imaginary unit i.
    pub const I: Complex = Complex { re: 0.0, im: 1.0 };

    /// Creates a new complex number.
    #[inline]
    pub const fn new(re: f64, im: f64) -> Self {
        Complex { re, im }
    }

    /// Creates a complex number from a real value.
    #[inline]
    pub const fn from_real(re: f64) -> Self {
        Complex { re, im: 0.0 }
    }

    /// Computes e^(i * theta) = cos(theta) + i*sin(theta).
    #[inline]
    pub fn exp_i(theta: f64) -> Self {
        Complex {
            re: theta.cos(),
            im: theta.sin(),
        }
    }

    /// Returns the complex conjugate.
    #[inline]
    pub fn conj(self) -> Self {
        Complex {
            re: self.re,
            im: -self.im,
        }
    }

    /// Returns the squared magnitude |z|^2 = re^2 + im^2.
    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Returns the magnitude |z|.
    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// Scales by a real number.
    #[inline]
    pub fn scale(self, s: f64) -> Self {
        Complex {
            re: self.re * s,
            im: self.im * s,
        }
    }

    /// Computes the multiplicative inverse 1/z.
    #[inline]
    pub fn inverse(self) -> Self {
        let norm_sq = self.norm_sq();
        Complex {
            re: self.re / norm_sq,
            im: -self.im / norm_sq,
        }
    }

    /// Divides by another complex number.
    #[inline]
    pub fn div(self, rhs: Complex) -> Self {
        self * rhs.inverse()
    }
}

impl Add for Complex {
    type Output = Complex;

    #[inline]
    fn add(self, rhs: Complex) -> Complex {
        Complex {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl AddAssign for Complex {
    #[inline]
    fn add_assign(&mut self, rhs: Complex) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

impl Sub for Complex {
    type Output = Complex;

    #[inline]
    fn sub(self, rhs: Complex) -> Complex {
        Complex {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl SubAssign for Complex {
    #[inline]
    fn sub_assign(&mut self, rhs: Complex) {
        self.re -= rhs.re;
        self.im -= rhs.im;
    }
}

impl Mul for Complex {
    type Output = Complex;

    #[inline]
    fn mul(self, rhs: Complex) -> Complex {
        Complex {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl MulAssign for Complex {
    #[inline]
    fn mul_assign(&mut self, rhs: Complex) {
        let re = self.re * rhs.re - self.im * rhs.im;
        let im = self.re * rhs.im + self.im * rhs.re;
        self.re = re;
        self.im = im;
    }
}

impl Neg for Complex {
    type Output = Complex;

    #[inline]
    fn neg(self) -> Complex {
        Complex {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl From<f64> for Complex {
    #[inline]
    fn from(re: f64) -> Self {
        Complex::from_real(re)
    }
}

// ============================================================================
// Roots of x^n + 1 in Falcon tree ordering
// ============================================================================

/// Computes the roots of x^n + 1 in the Falcon recursive tree ordering.
///
/// The ordering satisfies two key properties (matching tprest/falcon.py
/// `fft_constants.py`):
///
/// 1. `w[2i+1] = -w[2i]`  (adjacent pairs are negations)
/// 2. `w[i+n/2] = conj(w[i])`  (second half conjugates the first half)
///
/// Property 2 ensures that split_fft of a Hermitian-symmetric vector
/// produces Hermitian-symmetric children, propagating all the way down
/// to the leaves where values are purely real.
///
/// Uses iterative bottom-up construction to avoid deep recursion.
pub fn compute_roots(n: usize) -> Vec<Complex> {
    debug_assert!(n.is_power_of_two() && n >= 2);

    // Build iteratively from n=2 upward
    let mut current = vec![Complex::I, Complex::new(0.0, -1.0)]; // n=2 base
    let mut current_n = 2;

    while current_n < n {
        let next_n = current_n * 2;
        let hn = current_n; // = next_n / 2
        let mut roots = vec![Complex::ZERO; next_n];

        // First half: sqrt of parent roots (halving the argument)
        for i in 0..hn / 2 {
            let angle = current[i].im.atan2(current[i].re) / 2.0;
            roots[2 * i] = Complex::exp_i(angle);
            roots[2 * i + 1] = -roots[2 * i];
        }
        // Second half: conjugates of first half
        for i in 0..hn / 2 {
            roots[hn + 2 * i] = roots[2 * i].conj();
            roots[hn + 2 * i + 1] = -roots[hn + 2 * i];
        }

        current = roots;
        current_n = next_n;
    }

    current
}

// ============================================================================
// FFT Operations (Falcon recursive form)
// ============================================================================

/// Computes the forward negacyclic FFT.
///
/// Implements the Falcon reference recursive FFT for the ring Z[X]/(X^n + 1).
/// The output ordering matches the Falcon tree structure, ensuring that
/// adjacent entries (2i, 2i+1) are always conjugate pairs for real inputs.
/// This property is essential for ffSampling to produce correct signatures.
pub fn fft(a: &mut [Complex]) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    debug_assert!(n.is_power_of_two(), "FFT size must be power of 2");
    let roots = compute_roots(n);
    let result = fft_recursive(a, &roots, n);
    a.copy_from_slice(&result);
}

/// Recursive FFT implementation.
///
/// Splits coefficients into even/odd, recursively transforms each half,
/// then merges using the Falcon tree roots. Base case at n=2 directly
/// combines the two coefficients into a conjugate pair.
fn fft_recursive(coeffs: &[Complex], roots: &[Complex], n: usize) -> Vec<Complex> {
    if n == 2 {
        // Base case: [a, b] → [a + i*b, a - i*b]  (conjugate pair)
        return vec![
            Complex::new(coeffs[0].re - coeffs[1].im, coeffs[0].im + coeffs[1].re),
            Complex::new(coeffs[0].re + coeffs[1].im, coeffs[0].im - coeffs[1].re),
        ];
    }

    let hn = n / 2;
    // Split into even/odd coefficients
    let mut even = vec![Complex::ZERO; hn];
    let mut odd = vec![Complex::ZERO; hn];
    for i in 0..hn {
        even[i] = coeffs[2 * i];
        odd[i] = coeffs[2 * i + 1];
    }

    // Build child root tables
    let mut child_roots = vec![Complex::ZERO; hn];
    for i in 0..hn {
        // Child roots are the squares of the parent roots (at even indices)
        child_roots[i] = roots[2 * i] * roots[2 * i];
    }

    let f0 = fft_recursive(&even, &child_roots, hn);
    let f1 = fft_recursive(&odd, &child_roots, hn);

    // Merge: f_fft[2i] = f0[i] + w[2i]*f1[i], f_fft[2i+1] = f0[i] - w[2i]*f1[i]
    let mut result = vec![Complex::ZERO; n];
    for i in 0..hn {
        let t = roots[2 * i] * f1[i];
        result[2 * i] = f0[i] + t;
        result[2 * i + 1] = f0[i] - t;
    }
    result
}

/// Computes the inverse negacyclic FFT.
///
/// Reverses the Falcon recursive FFT, recovering real polynomial coefficients
/// from the FFT representation. Uses split_fft (FFT domain) and coefficient-
/// domain merge operations.
pub fn ifft(a: &mut [Complex]) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    debug_assert!(n.is_power_of_two(), "IFFT size must be power of 2");
    let roots = compute_roots(n);
    let result = ifft_recursive(a, &roots, n);
    a.copy_from_slice(&result);
}

/// Recursive inverse FFT implementation.
fn ifft_recursive(fft_vals: &[Complex], roots: &[Complex], n: usize) -> Vec<Complex> {
    if n == 2 {
        // Base case: [a+bi, a-bi] → [a, b]  (extract real and imaginary)
        return vec![
            Complex::from_real(fft_vals[0].re),
            Complex::from_real(fft_vals[0].im),
        ];
    }

    let hn = n / 2;

    // Split FFT representation into even/odd halves (interleaved)
    let mut f0 = vec![Complex::ZERO; hn];
    let mut f1 = vec![Complex::ZERO; hn];
    for i in 0..hn {
        f0[i] = (fft_vals[2 * i] + fft_vals[2 * i + 1]).scale(0.5);
        f1[i] = ((fft_vals[2 * i] - fft_vals[2 * i + 1]) * roots[2 * i].conj()).scale(0.5);
    }

    // Build child root tables
    let mut child_roots = vec![Complex::ZERO; hn];
    for i in 0..hn {
        child_roots[i] = roots[2 * i] * roots[2 * i];
    }

    let even = ifft_recursive(&f0, &child_roots, hn);
    let odd = ifft_recursive(&f1, &child_roots, hn);

    // Merge coefficients: interleave even and odd
    let mut result = vec![Complex::ZERO; n];
    for i in 0..hn {
        result[2 * i] = even[i];
        result[2 * i + 1] = odd[i];
    }
    result
}

// ============================================================================
// Split/Merge Operations (for FFT tree)
// ============================================================================

/// Splits an FFT representation into two halves (interleaved).
///
/// Uses the Falcon tree ordering where adjacent pairs (2i, 2i+1) are always
/// conjugate pairs for real polynomials. This ensures that children of a
/// Hermitian-symmetric parent are also Hermitian-symmetric, and at the
/// n=1 leaves, values are purely real.
///
/// The twiddle factor conj(w[2i]) is the conjugate of the tree root at
/// position 2i. For the tree root convention where w[2i+1] = -w[2i]:
///   f0[i] = (f[2i] + f[2i+1]) / 2
///   f1[i] = (f[2i] - f[2i+1]) * conj(w[2i]) / 2
pub fn split_fft(f: &[Complex]) -> (Vec<Complex>, Vec<Complex>) {
    let n = f.len();
    let hn = n / 2;

    let roots = compute_roots(n);
    let mut f0 = vec![Complex::ZERO; hn];
    let mut f1 = vec![Complex::ZERO; hn];

    for i in 0..hn {
        f0[i] = (f[2 * i] + f[2 * i + 1]).scale(0.5);
        f1[i] = ((f[2 * i] - f[2 * i + 1]) * roots[2 * i].conj()).scale(0.5);
    }

    (f0, f1)
}

/// Merges two FFT halves back into a full FFT representation (interleaved).
///
/// This is the inverse of split_fft.
pub fn merge_fft(f0: &[Complex], f1: &[Complex]) -> Vec<Complex> {
    let hn = f0.len();
    let n = hn * 2;

    let roots = compute_roots(n);
    let mut f = vec![Complex::ZERO; n];

    for i in 0..hn {
        let t = roots[2 * i] * f1[i];
        f[2 * i] = f0[i] + t;
        f[2 * i + 1] = f0[i] - t;
    }

    f
}

// ============================================================================
// Polynomial FFT Operations
// ============================================================================

/// Converts a polynomial from coefficient form to FFT form.
pub fn poly_to_fft(coeffs: &[i16], n: usize) -> Vec<Complex> {
    let mut f: Vec<Complex> = coeffs.iter().map(|&c| Complex::from_real(c as f64)).collect();
    f.resize(n, Complex::ZERO);
    fft(&mut f);
    f
}

/// Converts from FFT form back to (approximate) integer coefficients.
pub fn fft_to_poly(fft_coeffs: &[Complex]) -> Vec<i16> {
    let mut f = fft_coeffs.to_vec();
    ifft(&mut f);
    f.iter().map(|c| c.re.round() as i16).collect()
}

/// Multiplies two polynomials in FFT form (pointwise).
pub fn fft_mul(a: &[Complex], b: &[Complex]) -> Vec<Complex> {
    debug_assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).collect()
}

/// Adds two polynomials in FFT form.
pub fn fft_add(a: &[Complex], b: &[Complex]) -> Vec<Complex> {
    debug_assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect()
}

/// Subtracts two polynomials in FFT form.
pub fn fft_sub(a: &[Complex], b: &[Complex]) -> Vec<Complex> {
    debug_assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(&x, &y)| x - y).collect()
}

/// Negates a polynomial in FFT form.
pub fn fft_neg(a: &[Complex]) -> Vec<Complex> {
    a.iter().map(|&x| -x).collect()
}

/// Computes the adjoint (conjugate) of a polynomial in FFT form.
///
/// In the Falcon tree ordering, the adjoint permutation maps
/// index i to the index evaluating at the conjugate root.
/// Since w[i+n/2] = conj(w[i]), swapping the two halves with
/// conjugation gives the adjoint.
pub fn fft_adj(a: &[Complex]) -> Vec<Complex> {
    let n = a.len();
    let hn = n / 2;
    let mut result = vec![Complex::ZERO; n];
    for i in 0..hn {
        result[i] = a[i + hn].conj();
        result[i + hn] = a[i].conj();
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-9;

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < EPSILON
    }

    fn complex_approx_eq(a: Complex, b: Complex) -> bool {
        approx_eq(a.re, b.re) && approx_eq(a.im, b.im)
    }

    #[test]
    fn test_complex_add() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = a + b;
        assert!(approx_eq(c.re, 4.0));
        assert!(approx_eq(c.im, 6.0));
    }

    #[test]
    fn test_complex_mul() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        // (1+2i)(3+4i) = 3 + 4i + 6i + 8i^2 = 3 + 10i - 8 = -5 + 10i
        let c = a * b;
        assert!(approx_eq(c.re, -5.0));
        assert!(approx_eq(c.im, 10.0));
    }

    #[test]
    fn test_complex_conj() {
        let a = Complex::new(3.0, 4.0);
        let c = a.conj();
        assert!(approx_eq(c.re, 3.0));
        assert!(approx_eq(c.im, -4.0));
    }

    #[test]
    fn test_complex_norm() {
        let a = Complex::new(3.0, 4.0);
        assert!(approx_eq(a.norm_sq(), 25.0));
        assert!(approx_eq(a.norm(), 5.0));
    }

    #[test]
    fn test_complex_exp_i() {
        use std::f64::consts::PI;
        // e^(i*0) = 1
        let a = Complex::exp_i(0.0);
        assert!(complex_approx_eq(a, Complex::ONE));

        // e^(i*pi/2) = i
        let b = Complex::exp_i(PI / 2.0);
        assert!(complex_approx_eq(b, Complex::I));

        // e^(i*pi) = -1
        let c = Complex::exp_i(PI);
        assert!(complex_approx_eq(c, Complex::new(-1.0, 0.0)));
    }

    #[test]
    fn test_complex_inverse() {
        let a = Complex::new(3.0, 4.0);
        let a_inv = a.inverse();
        let product = a * a_inv;
        assert!(complex_approx_eq(product, Complex::ONE));
    }

    #[test]
    fn test_fft_ifft_roundtrip() {
        let n = 8;
        let original: Vec<Complex> = (0..n)
            .map(|i| Complex::from_real(i as f64))
            .collect();

        let mut transformed = original.clone();
        fft(&mut transformed);
        ifft(&mut transformed);

        for i in 0..n {
            assert!(
                (original[i].re - transformed[i].re).abs() < 1e-9 &&
                (original[i].im - transformed[i].im).abs() < 1e-9,
                "Mismatch at index {}: {:?} vs {:?}",
                i,
                original[i],
                transformed[i]
            );
        }
    }

    #[test]
    fn test_fft_ifft_roundtrip_512() {
        let n = 512;
        let original: Vec<Complex> = (0..n)
            .map(|i| Complex::from_real((i % 100) as f64 - 50.0))
            .collect();

        let mut transformed = original.clone();
        fft(&mut transformed);
        ifft(&mut transformed);

        for i in 0..n {
            assert!(
                (original[i].re - transformed[i].re).abs() < 1e-6,
                "Mismatch at index {}: {:?} vs {:?}",
                i,
                original[i],
                transformed[i]
            );
        }
    }

    #[test]
    fn test_split_merge_roundtrip() {
        let n = 8;
        let original: Vec<Complex> = (0..n)
            .map(|i| Complex::new(i as f64, (i * 2) as f64))
            .collect();

        let (f0, f1) = split_fft(&original);
        let merged = merge_fft(&f0, &f1);

        for i in 0..n {
            assert!(
                (original[i].re - merged[i].re).abs() < 1e-9 &&
                (original[i].im - merged[i].im).abs() < 1e-9,
                "Split/merge mismatch at {}: {:?} vs {:?}",
                i,
                original[i],
                merged[i]
            );
        }
    }

    #[test]
    fn test_poly_fft_roundtrip() {
        let coeffs: Vec<i16> = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let n = 8;

        let fft_form = poly_to_fft(&coeffs, n);
        let recovered = fft_to_poly(&fft_form);

        for i in 0..n {
            assert_eq!(
                coeffs[i], recovered[i],
                "Poly FFT roundtrip failed at {}",
                i
            );
        }
    }

    #[test]
    fn test_fft_mul_is_convolution() {
        // Test that FFT multiplication corresponds to polynomial multiplication
        // in the negacyclic ring Z[X]/(X^n + 1)
        let n = 4;

        // a(x) = 1 + 2x + 3x^2 + 4x^3
        let a_coeffs: Vec<i16> = vec![1, 2, 3, 4];
        // b(x) = 1 + x
        let b_coeffs: Vec<i16> = vec![1, 1, 0, 0];

        // Expected: a(x) * b(x) mod (x^4 + 1)
        // = (1 + 2x + 3x^2 + 4x^3)(1 + x) mod (x^4 + 1)
        // = 1 + x + 2x + 2x^2 + 3x^2 + 3x^3 + 4x^3 + 4x^4 mod (x^4 + 1)
        // = 1 + 3x + 5x^2 + 7x^3 + 4x^4 mod (x^4 + 1)
        // = 1 + 3x + 5x^2 + 7x^3 - 4 (since x^4 = -1)
        // = -3 + 3x + 5x^2 + 7x^3
        let expected: Vec<i16> = vec![-3, 3, 5, 7];

        let a_fft = poly_to_fft(&a_coeffs, n);
        let b_fft = poly_to_fft(&b_coeffs, n);
        let c_fft = fft_mul(&a_fft, &b_fft);
        let c_coeffs = fft_to_poly(&c_fft);

        for i in 0..n {
            assert_eq!(
                expected[i], c_coeffs[i],
                "Negacyclic multiplication failed at index {}",
                i
            );
        }
    }
}
