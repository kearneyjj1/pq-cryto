//! Complex floating-point FFT for FALCON.
//!
//! This module implements the Fast Fourier Transform over complex numbers,
//! specifically for the negacyclic polynomial ring Z[X]/(X^n + 1).
//!
//! FALCON uses floating-point FFT (not integer NTT) for:
//! - Efficient polynomial multiplication in key generation
//! - The Fast Fourier Sampling algorithm during signing

use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A complex number with f64 components.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
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
// FFT Operations
// ============================================================================

/// Computes the forward negacyclic FFT in-place.
///
/// For polynomials in Z[X]/(X^n + 1), we evaluate at the n-th roots of -1:
/// omega_k = exp(i * pi * (2k+1) / n) for k = 0..n-1.
///
/// Uses the twist method: pre-multiply by psi^j where psi = exp(i*pi/n),
/// then apply standard FFT.
pub fn fft(a: &mut [Complex]) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    debug_assert!(n.is_power_of_two(), "FFT size must be power of 2");

    // Pre-multiply by powers of psi = exp(i * pi / n) for negacyclic
    // psi is a primitive 2n-th root of unity, so psi^n = -1
    let psi = Complex::exp_i(PI / (n as f64));
    let mut psi_power = Complex::ONE;
    for i in 0..n {
        a[i] = a[i] * psi_power;
        psi_power = psi_power * psi;
    }

    // Standard Cooley-Tukey FFT
    fft_core(a);
}

/// Computes the inverse negacyclic FFT in-place.
///
/// Applies standard IFFT then post-multiplies by psi^(-j) to undo the twist.
pub fn ifft(a: &mut [Complex]) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    debug_assert!(n.is_power_of_two(), "IFFT size must be power of 2");

    // Standard inverse FFT
    ifft_core(a);

    // Post-multiply by inverse powers of psi = exp(-i * pi / n)
    let psi_inv = Complex::exp_i(-PI / (n as f64));
    let mut psi_power = Complex::ONE;
    for i in 0..n {
        a[i] = a[i] * psi_power;
        psi_power = psi_power * psi_inv;
    }
}

/// Standard forward FFT (Cooley-Tukey, in-place, radix-2).
fn fft_core(a: &mut [Complex]) {
    let n = a.len();
    let log_n = n.trailing_zeros() as usize;

    // Bit-reversal permutation
    for i in 0..n {
        let j = bit_reverse(i, log_n);
        if i < j {
            a.swap(i, j);
        }
    }

    // Cooley-Tukey butterflies
    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let angle = -2.0 * PI / (len as f64);
        let w_len = Complex::exp_i(angle);

        for start in (0..n).step_by(len) {
            let mut w = Complex::ONE;
            for j in 0..half {
                let u = a[start + j];
                let t = w * a[start + j + half];
                a[start + j] = u + t;
                a[start + j + half] = u - t;
                w = w * w_len;
            }
        }
        len *= 2;
    }
}

/// Standard inverse FFT (in-place).
fn ifft_core(a: &mut [Complex]) {
    let n = a.len();
    let log_n = n.trailing_zeros() as usize;

    // Bit-reversal permutation
    for i in 0..n {
        let j = bit_reverse(i, log_n);
        if i < j {
            a.swap(i, j);
        }
    }

    // Inverse butterflies (conjugate twiddle factors)
    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let angle = 2.0 * PI / (len as f64);  // Positive angle for inverse
        let w_len = Complex::exp_i(angle);

        for start in (0..n).step_by(len) {
            let mut w = Complex::ONE;
            for j in 0..half {
                let u = a[start + j];
                let t = w * a[start + j + half];
                a[start + j] = u + t;
                a[start + j + half] = u - t;
                w = w * w_len;
            }
        }
        len *= 2;
    }

    // Scale by 1/n
    let scale = 1.0 / (n as f64);
    for x in a.iter_mut() {
        *x = x.scale(scale);
    }
}

/// Bit-reversal of an index.
#[inline]
fn bit_reverse(mut x: usize, bits: usize) -> usize {
    let mut result = 0;
    for _ in 0..bits {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    result
}

// ============================================================================
// Split/Merge Operations (for FFT tree)
// ============================================================================

/// Splits an FFT representation into two halves.
///
/// Given f in FFT form, computes (f0, f1) where f(x) = f0(x^2) + x * f1(x^2).
/// This is the key operation for tree-based sampling.
pub fn split_fft(f: &[Complex]) -> (Vec<Complex>, Vec<Complex>) {
    let n = f.len();
    let hn = n / 2;

    let mut f0 = vec![Complex::ZERO; hn];
    let mut f1 = vec![Complex::ZERO; hn];

    for i in 0..hn {
        let a_plus = f[i];
        let a_minus = f[i + hn];

        // Twiddle factor for the split
        let angle = -PI * (i as f64) / (n as f64);
        let zeta = Complex::exp_i(angle);

        f0[i] = (a_plus + a_minus).scale(0.5);
        f1[i] = ((a_plus - a_minus) * zeta).scale(0.5);
    }

    (f0, f1)
}

/// Merges two FFT halves back into a full FFT representation.
///
/// This is the inverse of split_fft.
pub fn merge_fft(f0: &[Complex], f1: &[Complex]) -> Vec<Complex> {
    let hn = f0.len();
    let n = hn * 2;

    let mut f = vec![Complex::ZERO; n];

    for i in 0..hn {
        // Inverse twiddle factor
        let angle = PI * (i as f64) / (n as f64);
        let zeta = Complex::exp_i(angle);

        let t = f1[i] * zeta;
        f[i] = f0[i] + t;
        f[i + hn] = f0[i] - t;
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
pub fn fft_adj(a: &[Complex]) -> Vec<Complex> {
    let n = a.len();
    let mut result = vec![Complex::ZERO; n];
    result[0] = a[0].conj();
    for i in 1..n {
        result[i] = a[n - i].conj();
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
    fn test_bit_reverse() {
        assert_eq!(bit_reverse(0, 3), 0);
        assert_eq!(bit_reverse(1, 3), 4);
        assert_eq!(bit_reverse(2, 3), 2);
        assert_eq!(bit_reverse(3, 3), 6);
        assert_eq!(bit_reverse(4, 3), 1);
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
