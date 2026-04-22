//! Polynomial arithmetic over Z_q for FALCON.
//!
//! This module provides polynomial operations in the ring Z_q[X]/(X^n + 1)
//! where q = 12289 and n is 512 or 1024.

use crate::fft::{fft, ifft, Complex};
use crate::field::Zq;
use crate::params::Q;
use zeroize::Zeroize;

/// A polynomial in Z_q[X]/(X^n + 1).
///
/// Coefficients are stored in standard order: coeffs[i] is the coefficient of x^i.
#[derive(Clone, Debug, PartialEq)]
pub struct Poly {
    /// Polynomial coefficients in [0, q-1].
    pub coeffs: Vec<Zq>,
}

impl Drop for Poly {
    fn drop(&mut self) {
        self.coeffs.zeroize();
    }
}

impl Poly {
    /// Creates a new zero polynomial of degree n-1.
    pub fn zero(n: usize) -> Self {
        Poly {
            coeffs: vec![Zq::ZERO; n],
        }
    }

    /// Creates a polynomial from i16 coefficients.
    pub fn from_i16(coeffs: &[i16]) -> Self {
        Poly {
            coeffs: coeffs.iter().map(|&c| Zq::new(c as i32)).collect(),
        }
    }

    /// Creates a polynomial from Zq coefficients.
    pub fn from_zq(coeffs: Vec<Zq>) -> Self {
        Poly { coeffs }
    }

    /// Returns the number of coefficients.
    #[inline]
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Returns true if all coefficients are zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| c.is_zero())
    }

    /// Gets a coefficient.
    #[inline]
    pub fn get(&self, i: usize) -> Zq {
        self.coeffs[i]
    }

    /// Sets a coefficient.
    #[inline]
    pub fn set(&mut self, i: usize, val: Zq) {
        self.coeffs[i] = val;
    }

    /// Returns the coefficients as i16 values.
    pub fn to_i16(&self) -> Vec<i16> {
        self.coeffs.iter().map(|c| c.value()).collect()
    }

    /// Returns the centered coefficients in [-(q-1)/2, (q-1)/2].
    pub fn to_centered(&self) -> Vec<i16> {
        self.coeffs.iter().map(|c| c.centered()).collect()
    }

    /// Computes the squared L2 norm using centered representation.
    ///
    /// Uses saturating arithmetic to prevent overflow for very large polynomials.
    pub fn norm_sq(&self) -> i64 {
        self.coeffs
            .iter()
            .map(|c| {
                let v = c.centered() as i64;
                v * v
            })
            .fold(0i64, |acc, x| acc.saturating_add(x))
    }

    /// Computes the L2 norm.
    pub fn norm(&self) -> f64 {
        (self.norm_sq() as f64).sqrt()
    }

    /// Computes the infinity norm (max absolute coefficient).
    pub fn norm_inf(&self) -> i16 {
        self.coeffs
            .iter()
            .map(|c| c.centered().abs())
            .max()
            .unwrap_or(0)
    }

    /// Adds two polynomials coefficient-wise.
    pub fn add(&self, other: &Poly) -> Poly {
        debug_assert_eq!(self.len(), other.len());
        Poly {
            coeffs: self
                .coeffs
                .iter()
                .zip(other.coeffs.iter())
                .map(|(&a, &b)| a + b)
                .collect(),
        }
    }

    /// Subtracts two polynomials coefficient-wise.
    pub fn sub(&self, other: &Poly) -> Poly {
        debug_assert_eq!(self.len(), other.len());
        Poly {
            coeffs: self
                .coeffs
                .iter()
                .zip(other.coeffs.iter())
                .map(|(&a, &b)| a - b)
                .collect(),
        }
    }

    /// Negates the polynomial.
    pub fn neg(&self) -> Poly {
        Poly {
            coeffs: self.coeffs.iter().map(|&c| -c).collect(),
        }
    }

    /// Multiplies by a scalar.
    pub fn scale(&self, s: Zq) -> Poly {
        Poly {
            coeffs: self.coeffs.iter().map(|&c| c * s).collect(),
        }
    }

    /// Multiplies two polynomials in Z_q[X]/(X^n + 1) using FFT.
    ///
    /// Uses floating-point FFT — fast but subject to rounding error.
    /// For exact arithmetic, use [`mul_ntt`](Self::mul_ntt) instead.
    pub fn mul(&self, other: &Poly) -> Poly {
        debug_assert_eq!(self.len(), other.len());

        // Convert to FFT form
        let mut a_fft: Vec<Complex> = self
            .coeffs
            .iter()
            .map(|c| Complex::from_real(c.value() as f64))
            .collect();
        let mut b_fft: Vec<Complex> = other
            .coeffs
            .iter()
            .map(|c| Complex::from_real(c.value() as f64))
            .collect();

        fft(&mut a_fft);
        fft(&mut b_fft);

        // Pointwise multiply
        let mut c_fft: Vec<Complex> = a_fft
            .iter()
            .zip(b_fft.iter())
            .map(|(&a, &b)| a * b)
            .collect();

        // Convert back
        ifft(&mut c_fft);

        // Round and reduce with overflow protection
        Poly {
            coeffs: c_fft
                .iter()
                .map(|c| {
                    let rounded = c.re.round();
                    // Clamp to i32 range to prevent overflow on cast
                    let clamped = rounded.clamp(i32::MIN as f64, i32::MAX as f64) as i32;
                    Zq::new(clamped)
                })
                .collect(),
        }
    }

    /// Multiplies two polynomials in Z_q[X]/(X^n + 1) using NTT.
    ///
    /// Uses exact modular arithmetic via the Number Theoretic Transform.
    /// No floating-point rounding — results are exact in Z_q. This is the
    /// correct method for verification where exactness matters.
    ///
    /// Complexity: O(n log n) via Cooley-Tukey / Gentleman-Sande butterflies.
    pub fn mul_ntt(&self, other: &Poly) -> Poly {
        debug_assert_eq!(self.len(), other.len());
        let n = self.len();

        let a_ntt = ntt(&self.coeffs, n);
        let b_ntt = ntt(&other.coeffs, n);
        let c_ntt: Vec<Zq> = (0..n).map(|i| a_ntt[i] * b_ntt[i]).collect();

        Poly {
            coeffs: intt(&c_ntt, n),
        }
    }
}

// ============================================================================
// Polynomial FFT Representation (test-only)
// ============================================================================

/// A polynomial in FFT form (evaluated at the n-th roots of -1).
///
/// Operations in FFT form are pointwise and thus O(n) instead of O(n^2).
/// Only used in tests — production code works directly with `Vec<Complex>`.
#[cfg(test)]
#[derive(Clone, Debug)]
pub struct PolyFft {
    pub coeffs: Vec<Complex>,
}

#[cfg(test)]
impl Drop for PolyFft {
    fn drop(&mut self) {
        self.coeffs.zeroize();
    }
}

#[cfg(test)]
impl PolyFft {
    pub fn from_poly(p: &Poly) -> Self {
        let mut coeffs: Vec<Complex> = p
            .coeffs
            .iter()
            .map(|c| Complex::from_real(c.value() as f64))
            .collect();
        fft(&mut coeffs);
        PolyFft { coeffs }
    }

    pub fn to_poly(&self) -> Poly {
        let mut coeffs = self.coeffs.clone();
        ifft(&mut coeffs);
        Poly {
            coeffs: coeffs
                .iter()
                .map(|c| {
                    let rounded = c.re.round();
                    let clamped = rounded.clamp(i32::MIN as f64, i32::MAX as f64) as i32;
                    Zq::new(clamped)
                })
                .collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    pub fn mul(&self, other: &PolyFft) -> PolyFft {
        PolyFft {
            coeffs: self.coeffs.iter().zip(other.coeffs.iter()).map(|(&x, &y)| x * y).collect(),
        }
    }

    pub fn adj(&self) -> PolyFft {
        PolyFft {
            coeffs: self.coeffs.iter().map(|c| c.conj()).collect(),
        }
    }
}

// ============================================================================
// NTT Operations (for verification in Z_q)
// ============================================================================

/// Primitive n-th root of unity in Z_q for NTT.
/// For q = 12289 and n = 512, 1024, 2048, or 4096 we need roots of unity.
fn find_primitive_root(n: usize) -> Zq {
    // For q = 12289:
    // - Order of multiplicative group: q-1 = 12288 = 2^12 * 3
    // - For n = 1024: need 1024th root of unity (for negacyclic 512)
    // - For n = 2048: need 2048th root of unity (for negacyclic 1024)
    //
    // We need to find a generator of Z_q*.
    // g = 11 is a primitive root mod 12289.
    // omega_n = g^((q-1)/n) mod q

    assert!(
        (Q - 1) as usize % n == 0,
        "n={} must divide q-1={}",
        n,
        Q - 1
    );

    // 11 is a known primitive root mod 12289 (order = q-1 = 12288).
    // Verified: 11^((q-1)/p) != 1 for each prime factor p of 12288 = 2^12 * 3.
    // Uses pow (not pow_ct) because the exponent (q-1)/n is public.
    let g = Zq::new(11);
    let exp = (Q - 1) as u32 / (n as u32);
    g.pow(exp)
}

/// Reverses the bottom `log_n` bits of `x`.
#[inline]
fn bit_reverse(x: usize, log_n: u32) -> usize {
    x.reverse_bits() >> (usize::BITS - log_n)
}

/// Computes the negacyclic NTT in O(n log n) via Cooley-Tukey (DIT) butterflies.
///
/// Evaluates `a(X)` at `ψ, ψ³, ψ⁵, …, ψ^(2n-1)` where `ψ` is a primitive
/// 2n-th root of unity in Z_q. This is equivalent to:
/// 1. Twist: `a'[j] = a[j] * ψ^j`
/// 2. Standard DIT NTT on `a'` with `ω = ψ²`
pub fn ntt(a: &[Zq], n: usize) -> Vec<Zq> {
    let log_n = n.trailing_zeros();
    debug_assert_eq!(1 << log_n, n, "n must be a power of 2");

    let psi = find_primitive_root(2 * n); // 2n-th root of unity
    let omega = psi.square();             // n-th root of unity

    // Step 1: Pre-multiply by ψ^j (negacyclic twist)
    let mut tw = Vec::with_capacity(n);
    let mut psi_pow = Zq::ONE;
    for j in 0..n {
        tw.push(a[j] * psi_pow);
        psi_pow = psi_pow * psi;
    }

    // Step 2: Bit-reversal permutation
    let mut result = vec![Zq::ZERO; n];
    for i in 0..n {
        result[bit_reverse(i, log_n)] = tw[i];
    }

    // Step 3: Cooley-Tukey DIT butterflies
    let mut len = 1;
    while len < n {
        let step = omega.pow((n / (2 * len)) as u32);
        let mut k = 0;
        while k < n {
            let mut w = Zq::ONE;
            for j in 0..len {
                let u = result[k + j];
                let t = w * result[k + j + len];
                result[k + j] = u + t;
                result[k + j + len] = u - t;
                w = w * step;
            }
            k += 2 * len;
        }
        len <<= 1;
    }

    result
}

/// Computes the inverse negacyclic NTT in O(n log n) via Gentleman-Sande (DIF) butterflies.
///
/// Recovers coefficients from NTT-domain values:
/// 1. Gentleman-Sande DIF butterflies with `ω⁻¹ = ψ⁻²`
/// 2. Bit-reversal permutation
/// 3. Un-twist: `a[j] *= ψ⁻ʲ · n⁻¹`
pub fn intt(a: &[Zq], n: usize) -> Vec<Zq> {
    let log_n = n.trailing_zeros();
    debug_assert_eq!(1 << log_n, n, "n must be a power of 2");

    let psi = find_primitive_root(2 * n);
    let omega_inv = psi.square().inverse(); // ω⁻¹
    let psi_inv = psi.inverse();
    let n_inv = Zq::new(n as i32).inverse();

    let mut result = a.to_vec();

    // Step 1: Gentleman-Sande DIF butterflies
    let mut len = n / 2;
    while len >= 1 {
        let step = omega_inv.pow((n / (2 * len)) as u32);
        let mut k = 0;
        while k < n {
            let mut w = Zq::ONE;
            for j in 0..len {
                let u = result[k + j];
                let t = result[k + j + len];
                result[k + j] = u + t;
                result[k + j + len] = (u - t) * w;
                w = w * step;
            }
            k += 2 * len;
        }
        len >>= 1;
    }

    // Step 2: Bit-reversal permutation
    let mut permuted = vec![Zq::ZERO; n];
    for i in 0..n {
        permuted[bit_reverse(i, log_n)] = result[i];
    }

    // Step 3: Un-twist and divide by n: a[j] *= ψ⁻ʲ · n⁻¹
    let mut psi_inv_pow = Zq::ONE;
    for j in 0..n {
        permuted[j] = permuted[j] * psi_inv_pow * n_inv;
        psi_inv_pow = psi_inv_pow * psi_inv;
    }

    permuted
}

/// O(n²) naive negacyclic NTT, retained for comparison testing.
#[cfg(test)]
fn ntt_naive(a: &[Zq], n: usize) -> Vec<Zq> {
    let omega = find_primitive_root(2 * n);
    let mut result = vec![Zq::ZERO; n];
    for i in 0..n {
        let mut sum = Zq::ZERO;
        let omega_i = omega.pow((2 * i + 1) as u32);
        let mut omega_ij = Zq::ONE;
        for j in 0..n {
            sum += a[j] * omega_ij;
            omega_ij = omega_ij * omega_i;
        }
        result[i] = sum;
    }
    result
}

/// O(n²) naive inverse negacyclic NTT, retained for comparison testing.
#[cfg(test)]
fn intt_naive(a: &[Zq], n: usize) -> Vec<Zq> {
    let omega = find_primitive_root(2 * n);
    let omega_inv = omega.inverse();
    let n_inv = Zq::new(n as i32).inverse();
    let mut result = vec![Zq::ZERO; n];
    for i in 0..n {
        let mut sum = Zq::ZERO;
        for k in 0..n {
            let exp = ((2 * k + 1) * i) as u32;
            sum += a[k] * omega_inv.pow(exp);
        }
        result[i] = sum * n_inv;
    }
    result
}

/// Multiplies two polynomials in Z_q[X]/(X^n + 1) using schoolbook method.
/// Used for testing and small polynomials.
pub fn poly_mul_schoolbook(a: &[Zq], b: &[Zq]) -> Vec<Zq> {
    let n = a.len();
    debug_assert_eq!(n, b.len());

    let mut result = vec![Zq::ZERO; n];

    for i in 0..n {
        for j in 0..n {
            let idx = i + j;
            if idx < n {
                result[idx] += a[i] * b[j];
            } else {
                // x^n = -1, so x^(n+k) = -x^k
                result[idx - n] -= a[i] * b[j];
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_add() {
        let a = Poly::from_i16(&[1, 2, 3, 4]);
        let b = Poly::from_i16(&[5, 6, 7, 8]);
        let c = a.add(&b);
        assert_eq!(c.to_i16(), vec![6, 8, 10, 12]);
    }

    #[test]
    fn test_poly_sub() {
        let a = Poly::from_i16(&[10, 20, 30, 40]);
        let b = Poly::from_i16(&[1, 2, 3, 4]);
        let c = a.sub(&b);
        assert_eq!(c.to_i16(), vec![9, 18, 27, 36]);
    }

    #[test]
    fn test_poly_neg() {
        let a = Poly::from_i16(&[1, 2, 3, 4]);
        let neg_a = a.neg();
        let sum = a.add(&neg_a);
        assert!(sum.is_zero());
    }

    #[test]
    fn test_poly_scale() {
        let a = Poly::from_i16(&[1, 2, 3, 4]);
        let b = a.scale(Zq::new(3));
        assert_eq!(b.to_i16(), vec![3, 6, 9, 12]);
    }

    #[test]
    fn test_poly_mul_simple() {
        // (1 + x) * (1 - x) = 1 - x^2 in normal polynomial ring
        // But in Z[X]/(X^n + 1), we need to consider wrapping
        let _n = 4;
        let a = Poly::from_i16(&[1, 1, 0, 0]);
        let b = Poly::from_i16(&[1, -1 + Q as i16, 0, 0]); // 1 - x, using positive representation
        let c = a.mul(&b);

        // (1 + x)(1 - x) = 1 - x^2
        let expected = Poly::from_i16(&[1, 0, Q as i16 - 1, 0]); // 1 - x^2
        assert_eq!(c.coeffs, expected.coeffs);
    }

    #[test]
    fn test_poly_mul_matches_schoolbook() {
        // Test that FFT multiplication matches schoolbook
        let a_coeffs: Vec<Zq> = vec![1, 2, 3, 4, 5, 6, 7, 8]
            .into_iter()
            .map(Zq::new)
            .collect();
        let b_coeffs: Vec<Zq> = vec![8, 7, 6, 5, 4, 3, 2, 1]
            .into_iter()
            .map(Zq::new)
            .collect();

        let a = Poly::from_zq(a_coeffs.clone());
        let b = Poly::from_zq(b_coeffs.clone());

        let fft_result = a.mul(&b);
        let schoolbook_result = poly_mul_schoolbook(&a_coeffs, &b_coeffs);

        for i in 0..8 {
            assert_eq!(
                fft_result.coeffs[i], schoolbook_result[i],
                "Mismatch at index {}: fft={:?}, schoolbook={:?}",
                i, fft_result.coeffs[i], schoolbook_result[i]
            );
        }
    }

    #[test]
    fn test_poly_mul_ntt_matches_schoolbook() {
        // Test that NTT multiplication matches schoolbook (exact)
        let a_coeffs: Vec<Zq> = vec![1, 2, 3, 4, 5, 6, 7, 8]
            .into_iter()
            .map(Zq::new)
            .collect();
        let b_coeffs: Vec<Zq> = vec![8, 7, 6, 5, 4, 3, 2, 1]
            .into_iter()
            .map(Zq::new)
            .collect();

        let a = Poly::from_zq(a_coeffs.clone());
        let b = Poly::from_zq(b_coeffs.clone());

        let ntt_result = a.mul_ntt(&b);
        let schoolbook_result = poly_mul_schoolbook(&a_coeffs, &b_coeffs);

        for i in 0..8 {
            assert_eq!(
                ntt_result.coeffs[i], schoolbook_result[i],
                "NTT mismatch at index {}: ntt={:?}, schoolbook={:?}",
                i, ntt_result.coeffs[i], schoolbook_result[i]
            );
        }
    }

    #[test]
    fn test_poly_fft_roundtrip() {
        let p = Poly::from_i16(&[1, 2, 3, 4, 5, 6, 7, 8]);
        let fft = PolyFft::from_poly(&p);
        let p2 = fft.to_poly();
        assert_eq!(p.coeffs, p2.coeffs);
    }

    #[test]
    fn test_poly_fft_mul() {
        let a = Poly::from_i16(&[1, 2, 3, 4, 5, 6, 7, 8]);
        let b = Poly::from_i16(&[8, 7, 6, 5, 4, 3, 2, 1]);

        // Multiply in coefficient form
        let c_coeff = a.mul(&b);

        // Multiply in FFT form
        let a_fft = PolyFft::from_poly(&a);
        let b_fft = PolyFft::from_poly(&b);
        let c_fft = a_fft.mul(&b_fft);
        let c_from_fft = c_fft.to_poly();

        assert_eq!(c_coeff.coeffs, c_from_fft.coeffs);
    }

    #[test]
    fn test_poly_norm() {
        let p = Poly::from_i16(&[1, -1 + Q as i16, 2, -2 + Q as i16]); // [1, -1, 2, -2]
        assert_eq!(p.norm_sq(), 1 + 1 + 4 + 4); // 10
        assert_eq!(p.norm_inf(), 2);
    }

    #[test]
    fn test_schoolbook_negacyclic() {
        // Test schoolbook multiplication handles x^n = -1 correctly
        let n = 4;

        // x^3 * x = x^4 = -1
        let a: Vec<Zq> = vec![0, 0, 0, 1].into_iter().map(Zq::new).collect();
        let b: Vec<Zq> = vec![0, 1, 0, 0].into_iter().map(Zq::new).collect();
        let c = poly_mul_schoolbook(&a, &b);

        // x^3 * x = x^4 = -1 in Z[X]/(X^4 + 1)
        assert_eq!(c[0], Zq::new(Q - 1)); // -1
        for i in 1..n {
            assert!(c[i].is_zero());
        }
    }

    #[test]
    fn test_poly_fft_adj() {
        let p = Poly::from_i16(&[1, 2, 3, 4, 5, 6, 7, 8]);
        let p_fft = PolyFft::from_poly(&p);
        let adj = p_fft.adj();

        // The adjoint in FFT form is just the complex conjugate
        for i in 0..p_fft.len() {
            assert_eq!(adj.coeffs[i], p_fft.coeffs[i].conj());
        }
    }

    #[test]
    fn test_primitive_root() {
        // For negacyclic NTT of size n, we need 2n-th roots of unity
        // Check that omega^(2n) = 1 and omega^n != 1
        for &n in &[512, 1024] {
            let omega = find_primitive_root(2 * n);

            // omega^(2n) should be 1
            assert_eq!(
                omega.pow((2 * n) as u32),
                Zq::ONE,
                "omega^(2n) should be 1 for n={}",
                n
            );

            // omega^n should NOT be 1 (primitive 2n-th root)
            assert_ne!(
                omega.pow(n as u32),
                Zq::ONE,
                "omega should be primitive 2n-th root for n={}",
                n
            );

            // omega^n should be -1 (for negacyclic: x^n = -1)
            assert_eq!(
                omega.pow(n as u32),
                Zq::new(Q - 1),
                "omega^n should be -1 for negacyclic n={}",
                n
            );
        }
    }

    #[test]
    fn test_ntt_fast_vs_naive() {
        // Compare O(n log n) Cooley-Tukey NTT against O(n²) naive at several sizes
        for &n in &[8, 16, 64, 512] {
            let a: Vec<Zq> = (0..n).map(|i| Zq::new((i * 37 + 13) as i32 % Q)).collect();
            let fast = ntt(&a, n);
            let naive = ntt_naive(&a, n);
            assert_eq!(fast, naive, "ntt mismatch at n={}", n);
        }
    }

    #[test]
    fn test_intt_fast_vs_naive() {
        // Compare O(n log n) Gentleman-Sande INTT against O(n²) naive at several sizes
        for &n in &[8, 16, 64, 512] {
            let a: Vec<Zq> = (0..n).map(|i| Zq::new((i * 53 + 7) as i32 % Q)).collect();
            let fast = intt(&a, n);
            let naive = intt_naive(&a, n);
            assert_eq!(fast, naive, "intt mismatch at n={}", n);
        }
    }

    #[test]
    fn test_ntt_intt_roundtrip() {
        // Verify ntt then intt recovers original polynomial
        for &n in &[8, 16, 64, 512, 1024] {
            let a: Vec<Zq> = (0..n).map(|i| Zq::new((i * 41 + 3) as i32 % Q)).collect();
            let transformed = ntt(&a, n);
            let recovered = intt(&transformed, n);
            assert_eq!(a, recovered, "roundtrip failed at n={}", n);
        }
    }

    #[test]
    fn test_ntt_mul_matches_schoolbook() {
        // Verify fast NTT-based multiplication matches schoolbook at n=512
        let n = 512;
        let a: Vec<Zq> = (0..n).map(|i| Zq::new((i * 17 + 5) as i32 % Q)).collect();
        let b: Vec<Zq> = (0..n).map(|i| Zq::new((i * 23 + 11) as i32 % Q)).collect();

        let a_ntt = ntt(&a, n);
        let b_ntt = ntt(&b, n);
        let c_ntt: Vec<Zq> = (0..n).map(|i| a_ntt[i] * b_ntt[i]).collect();
        let ntt_result = intt(&c_ntt, n);

        let schoolbook = poly_mul_schoolbook(&a, &b);
        assert_eq!(ntt_result, schoolbook, "NTT mul != schoolbook at n=512");
    }
}
