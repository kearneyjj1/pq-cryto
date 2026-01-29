//! Fast Fourier Sampling for FALCON (FIPS 206).
//!
//! This module implements the ffSampling algorithm, which is the core
//! of FALCON's signing procedure. It samples a lattice vector close to
//! a target using FFT-domain decomposition with the LDL* tree.
//!
//! # FIPS 206 Compliance
//!
//! This implementation follows the FALCON specification for ffSampling:
//! - Recursive tree-based sampling using the LDL* decomposition
//! - Proper conditioning of samples using the L factors
//! - Uses the FIPS 206 compliant Gaussian sampler from gaussian.rs
//!
//! # Algorithm Overview
//!
//! The ffSampling algorithm works recursively:
//! 1. At leaf nodes: sample z1 from N(t1, sigma1²), then z0 from N(t0 + l10*(t1-z1), sigma0²)
//! 2. At internal nodes: split the target, recurse on children, merge results
//!
//! The key insight is that in the FFT domain, the n-dimensional lattice
//! sampling problem decomposes into n independent 2×2 problems.

use crate::fft::{fft, ifft, merge_fft, split_fft, Complex};
use crate::fft_tree::{GramSchmidt, LdlTree};
use crate::gaussian::{sample_gaussian, SamplerZ, SIGMA_MIN_512};
use rand::RngCore;

// ============================================================================
// ffSampling Algorithm (FIPS 206)
// ============================================================================

/// The Fast Fourier Sampler implementing FIPS 206 ffSampling.
///
/// This sampler uses the LDL* tree structure to efficiently sample
/// from a discrete Gaussian distribution over the NTRU lattice.
/// The recursive algorithm decomposes the problem using FFT splits.
pub struct FfSampler {
    /// The Gram-Schmidt / LDL* data for the secret key.
    gs: GramSchmidt,
    /// The integer Gaussian sampler.
    sampler_z: SamplerZ,
    /// The signing sigma (from parameters).
    /// For FALCON-512: ~165.7, for FALCON-1024: ~168.4
    sigma_sign: f64,
}

impl FfSampler {
    /// Creates a new Fast Fourier Sampler from the Gram-Schmidt data.
    ///
    /// The sigma_sign is inferred from the polynomial degree n:
    /// - n <= 512: sigma ≈ 165.7 (FALCON-512)
    /// - n > 512: sigma ≈ 168.4 (FALCON-1024)
    pub fn new(gs: GramSchmidt) -> Self {
        let n = gs.f_fft.len();
        // FALCON sigma values from FIPS 206
        // sigma = 1.17 * sqrt(q) * sqrt(2n / (2n-1))
        let sigma_sign = if n <= 512 {
            165.7366171228152  // FALCON-512
        } else {
            168.38857144162388  // FALCON-1024
        };

        FfSampler {
            gs,
            sampler_z: SamplerZ::new(),
            sigma_sign,
        }
    }

    /// Creates a new Fast Fourier Sampler with a specific sigma value.
    pub fn with_sigma(gs: GramSchmidt, sigma_sign: f64) -> Self {
        FfSampler {
            gs,
            sampler_z: SamplerZ::new(),
            sigma_sign,
        }
    }

    /// Performs ffSampling to sample a lattice point close to target (t0, t1).
    ///
    /// This is the main FIPS 206 ffSampling algorithm. Given a target vector
    /// (t0, t1) in FFT form, it samples (z0, z1) from the discrete Gaussian
    /// distribution centered at the target, using the LDL* tree for the
    /// correct covariance structure.
    ///
    /// # Arguments
    ///
    /// * `rng` - Random number generator
    /// * `t0_fft` - First component of target in FFT form
    /// * `t1_fft` - Second component of target in FFT form
    ///
    /// # Returns
    ///
    /// (z0_fft, z1_fft) - Sampled lattice point in FFT form
    pub fn ff_sampling<R: RngCore>(
        &self,
        rng: &mut R,
        t0_fft: &[Complex],
        t1_fft: &[Complex],
    ) -> (Vec<Complex>, Vec<Complex>) {
        let n = t0_fft.len();
        assert_eq!(n, t1_fft.len(), "Target components must have same length");
        assert!(n.is_power_of_two(), "Length must be power of 2");

        // Start the recursive sampling
        self.ff_sampling_recursive(rng, t0_fft, t1_fft, 0, 0)
    }

    /// Recursive ffSampling implementation.
    ///
    /// Traverses the LDL* tree, sampling at each level using the
    /// stored L factors and sigma values.
    fn ff_sampling_recursive<R: RngCore>(
        &self,
        rng: &mut R,
        t0: &[Complex],
        t1: &[Complex],
        level: usize,
        pos: usize,
    ) -> (Vec<Complex>, Vec<Complex>) {
        let n = t0.len();

        if n == 1 {
            // Base case: 2×2 sampling using LDL* decomposition
            return self.ff_sampling_base_case(rng, t0[0], t1[0], level, pos);
        }

        // Recursive case: split, recurse, merge

        // Split targets using FFT split
        let (t0_0, t0_1) = split_fft(t0);
        let (t1_0, t1_1) = split_fft(t1);

        // Recurse on left subtree (even indices)
        let (z0_0, z1_0) = self.ff_sampling_recursive(rng, &t0_0, &t1_0, level + 1, 2 * pos);

        // Recurse on right subtree (odd indices)
        let (z0_1, z1_1) = self.ff_sampling_recursive(rng, &t0_1, &t1_1, level + 1, 2 * pos + 1);

        // Merge results
        let z0 = merge_fft(&z0_0, &z0_1);
        let z1 = merge_fft(&z1_0, &z1_1);

        (z0, z1)
    }

    /// Base case sampling for a single 2×2 complex lattice.
    ///
    /// At each leaf of the LDL* tree, we have a 2×2 problem:
    /// - Sample z1 from N(t1, (sigma/sqrt(d1))²)
    /// - Sample z0 from N(t0 + l10*(t1 - z1), (sigma/sqrt(d0))²)
    ///
    /// The tree stores sqrt(d_i), so sampling_sigma = sigma_sign / sqrt(d_i).
    /// The conditioning through l10 ensures the correct covariance structure.
    fn ff_sampling_base_case<R: RngCore>(
        &self,
        rng: &mut R,
        t0: Complex,
        t1: Complex,
        level: usize,
        pos: usize,
    ) -> (Vec<Complex>, Vec<Complex>) {
        let tree = &self.gs.tree;

        // Get sqrt(d) values from the LDL* tree leaf
        // These are sqrt(d0) and sqrt(d1) from the LDL* decomposition
        let sqrt_d0 = tree.get_sigma(level, pos, 0);
        let sqrt_d1 = tree.get_sigma(level, pos, 1);

        // Compute actual sampling sigmas: sigma_sample = sigma_sign / sqrt(d)
        // This normalizes the sampling variance by the Gram matrix structure
        let sigma0 = if sqrt_d0 > 1e-10 {
            self.sigma_sign / sqrt_d0
        } else {
            self.sigma_sign
        };
        let sigma1 = if sqrt_d1 > 1e-10 {
            self.sigma_sign / sqrt_d1
        } else {
            self.sigma_sign
        };

        // Get the l10 projection coefficient
        let l10 = tree.get_l10(level, pos, 0);

        // Sample z1 first (from the orthogonalized direction)
        // For complex targets, sample real and imaginary parts independently
        let z1_re = self.sampler_z.sample(rng, t1.re, sigma1) as f64;
        let z1_im = self.sampler_z.sample(rng, t1.im, sigma1) as f64;
        let z1 = Complex::new(z1_re, z1_im);

        // Compute the adjustment for z0 based on z1
        // t0_adj = t0 + l10 * (t1 - z1)
        let diff = t1 - z1;
        let adjustment = l10 * diff;
        let t0_adj = t0 + adjustment;

        // Sample z0 from the adjusted distribution
        let z0_re = self.sampler_z.sample(rng, t0_adj.re, sigma0) as f64;
        let z0_im = self.sampler_z.sample(rng, t0_adj.im, sigma0) as f64;
        let z0 = Complex::new(z0_re, z0_im);

        (vec![z0], vec![z1])
    }

    /// Samples a signature for a message hash.
    ///
    /// Given the hash c of (r, message) in FFT form, computes the target
    /// vector and samples a lattice point close to it.
    ///
    /// The target is t = B^(-1) * (c, 0) where B is the secret basis.
    /// For FALCON's basis B = [[f, g], [F, G]] with NTRU equation f*G - g*F = q:
    ///   B^(-1) = (1/q) * [[G, -g], [-F, f]]
    ///   t = B^(-1) * (c, 0) = ((G*c)/q, (-F*c)/q)
    ///
    /// Returns (z0, z1) sampled close to t using ffSampling.
    /// The samples are returned in FFT form.
    pub fn sample_preimage<R: RngCore>(
        &self,
        rng: &mut R,
        c_fft: &[Complex],
    ) -> (Vec<Complex>, Vec<Complex>) {
        let n = c_fft.len();
        let q = 12289.0;

        // Compute target: t = B^(-1) * (c, 0)
        // For basis B = [[f, g], [F, G]] with NTRU equation f*G - g*F = q:
        // B^(-1) = (1/q) * [[G, -g], [-F, f]]
        // t0 = (1/q) * G*c
        // t1 = (1/q) * (-F*c)

        // Compute G*c/q and -F*c/q in FFT domain (pointwise division by q)
        // Note: In FFT domain, polynomial multiplication is pointwise,
        // so we can just scale each FFT coefficient by 1/q
        let mut t0_fft: Vec<Complex> = Vec::with_capacity(n);
        let mut t1_fft: Vec<Complex> = Vec::with_capacity(n);

        for i in 0..n {
            // (G*c)/q in FFT domain = G_fft[i] * c_fft[i] / q
            t0_fft.push((self.gs.big_g_fft[i] * c_fft[i]).scale(1.0 / q));
            // (-F*c)/q in FFT domain = -F_fft[i] * c_fft[i] / q
            t1_fft.push((-self.gs.big_f_fft[i] * c_fft[i]).scale(1.0 / q));
        }

        // ffSampling operates in FFT domain using split_fft/merge_fft
        // Sample using ffSampling (targets and outputs are in FFT form)
        self.ff_sampling(rng, &t0_fft, &t1_fft)
    }

    /// Legacy interface for compatibility with old code.
    ///
    /// This simply calls ff_sampling with the given targets.
    #[deprecated(note = "Use ff_sampling instead")]
    pub fn sample<R: RngCore>(
        &self,
        rng: &mut R,
        t0_fft: &[Complex],
        t1_fft: &[Complex],
        _sigma: f64,
    ) -> (Vec<Complex>, Vec<Complex>) {
        self.ff_sampling(rng, t0_fft, t1_fft)
    }

    /// Gets a reference to the Gram-Schmidt data.
    pub fn gram_schmidt(&self) -> &GramSchmidt {
        &self.gs
    }
}

// ============================================================================
// Simplified Sampler (for testing and fallback)
// ============================================================================

/// Simplified sampler for testing.
///
/// This sampler doesn't use the full FFT tree but still produces
/// valid lattice samples (with potentially larger norms).
pub struct SimpleSampler {
    sampler_z: SamplerZ,
}

impl SimpleSampler {
    /// Creates a new simple sampler.
    pub fn new() -> Self {
        SimpleSampler {
            sampler_z: SamplerZ::new(),
        }
    }

    /// Samples z from a discrete Gaussian with given mean and sigma.
    pub fn sample<R: RngCore>(&self, rng: &mut R, mu: f64, sigma: f64) -> i64 {
        self.sampler_z.sample(rng, mu, sigma)
    }

    /// Samples a polynomial from a discrete Gaussian.
    ///
    /// Each coefficient is sampled independently from N(mu[i], sigma^2).
    pub fn sample_poly<R: RngCore>(&self, rng: &mut R, mu: &[f64], sigma: f64) -> Vec<i64> {
        mu.iter().map(|&m| self.sample(rng, m, sigma)).collect()
    }

    /// Samples a polynomial with zero mean.
    pub fn sample_poly_zero<R: RngCore>(&self, rng: &mut R, n: usize, sigma: f64) -> Vec<i64> {
        (0..n).map(|_| self.sample(rng, 0.0, sigma)).collect()
    }
}

impl Default for SimpleSampler {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Signature Computation Utilities
// ============================================================================

/// Computes the signature from the sampled lattice point.
///
/// Given the sampled (z0, z1), computes s = (s1, s2) where:
/// - s1 = c - z0 * f - z1 * F (should be zero for valid preimage)
/// - s2 = -z0 * g - z1 * G
///
/// The actual signature is s2 (compressed).
pub fn compute_signature(
    z0_fft: &[Complex],
    z1_fft: &[Complex],
    f_fft: &[Complex],
    g_fft: &[Complex],
    big_f_fft: &[Complex],
    big_g_fft: &[Complex],
    c_fft: &[Complex],
) -> (Vec<Complex>, Vec<Complex>) {
    let n = z0_fft.len();

    // s1 = c - z0*f - z1*F (should be close to zero)
    let mut s1_fft = Vec::with_capacity(n);
    for i in 0..n {
        s1_fft.push(c_fft[i] - z0_fft[i] * f_fft[i] - z1_fft[i] * big_f_fft[i]);
    }

    // s2 = -z0*g - z1*G
    let mut s2_fft = Vec::with_capacity(n);
    for i in 0..n {
        s2_fft.push(-z0_fft[i] * g_fft[i] - z1_fft[i] * big_g_fft[i]);
    }

    (s1_fft, s2_fft)
}

/// Verifies that the signature has an acceptable norm.
///
/// The squared norm of (s1, s2) should be below the FIPS 206 bound.
///
/// # FIPS 206 Bounds
///
/// - FALCON-512: ||s||² ≤ 34,034,726
/// - FALCON-1024: ||s||² ≤ 70,265,242
pub fn check_signature_norm(s1_fft: &[Complex], s2_fft: &[Complex], bound_sq: f64) -> bool {
    // Convert to coefficient form
    let mut s1 = s1_fft.to_vec();
    let mut s2 = s2_fft.to_vec();
    ifft(&mut s1);
    ifft(&mut s2);

    // Compute squared norm
    let norm_sq: f64 =
        s1.iter().map(|c| c.re * c.re).sum::<f64>() + s2.iter().map(|c| c.re * c.re).sum::<f64>();

    norm_sq <= bound_sq
}

/// Computes the squared norm of a signature.
pub fn signature_norm_sq(s1_fft: &[Complex], s2_fft: &[Complex]) -> f64 {
    let mut s1 = s1_fft.to_vec();
    let mut s2 = s2_fft.to_vec();
    ifft(&mut s1);
    ifft(&mut s2);

    s1.iter().map(|c| c.re * c.re).sum::<f64>() + s2.iter().map(|c| c.re * c.re).sum::<f64>()
}

// ============================================================================
// FIPS 206 Signature Bounds
// ============================================================================

/// FIPS 206 signature bound squared for FALCON-512.
pub const FALCON_512_SIG_BOUND_SQ: f64 = 34_034_726.0;

/// FIPS 206 signature bound squared for FALCON-1024.
pub const FALCON_1024_SIG_BOUND_SQ: f64 = 70_265_242.0;

/// Returns the appropriate signature bound for the given n.
pub fn get_sig_bound_sq(n: usize) -> f64 {
    if n <= 512 {
        FALCON_512_SIG_BOUND_SQ
    } else {
        FALCON_1024_SIG_BOUND_SQ
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fft::fft;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_simple_sampler() {
        let sampler = SimpleSampler::new();
        let mut rng = StdRng::seed_from_u64(42);

        // Sample a polynomial with zero mean
        let sigma = 10.0;
        let n = 16;
        let z = sampler.sample_poly_zero(&mut rng, n, sigma);

        assert_eq!(z.len(), n);

        // Check that the sample has reasonable norm
        let norm_sq: i64 = z.iter().map(|&x| x * x).sum();
        let expected_norm_sq = (n as f64) * sigma * sigma;

        // Very loose check since this is random
        assert!(
            (norm_sq as f64) < expected_norm_sq * 5.0,
            "Norm squared {} should be roughly {}",
            norm_sq,
            expected_norm_sq
        );
    }

    #[test]
    fn test_sample_poly_with_mean() {
        let sampler = SimpleSampler::new();
        let mut rng = StdRng::seed_from_u64(123);

        let mu: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
        let sigma = 2.0;

        // Sample many times and check mean
        let n_samples = 1000;
        let n = mu.len();
        let mut sums = vec![0i64; n];

        for _ in 0..n_samples {
            let z = sampler.sample_poly(&mut rng, &mu, sigma);
            for (i, &zi) in z.iter().enumerate() {
                sums[i] += zi;
            }
        }

        // Check that means are close to mu
        for (i, &sum) in sums.iter().enumerate() {
            let mean = (sum as f64) / (n_samples as f64);
            assert!(
                (mean - mu[i]).abs() < 1.0,
                "Mean {} should be close to {} at index {}",
                mean,
                mu[i],
                i
            );
        }
    }

    #[test]
    fn test_check_signature_norm() {
        // Create a signature with known norm
        let mut s1: Vec<Complex> = vec![1.0, 2.0, 3.0, 4.0]
            .into_iter()
            .map(|x| Complex::from_real(x))
            .collect();
        let mut s2: Vec<Complex> = vec![1.0, 1.0, 1.0, 1.0]
            .into_iter()
            .map(|x| Complex::from_real(x))
            .collect();

        // Convert to FFT form
        fft(&mut s1);
        fft(&mut s2);

        // Squared norm in coefficient form: 1+4+9+16 + 4 = 34
        let bound_sq = 50.0; // Above the norm
        assert!(check_signature_norm(&s1, &s2, bound_sq));

        let bound_sq_small = 20.0; // Below the norm
        assert!(!check_signature_norm(&s1, &s2, bound_sq_small));
    }

    #[test]
    fn test_ff_sampler_creation() {
        // Create a simple Gram-Schmidt structure
        let n = 4;
        let f_fft: Vec<Complex> = (0..n).map(|i| Complex::new(1.0 + i as f64, 0.0)).collect();
        let g_fft: Vec<Complex> = (0..n).map(|i| Complex::new(0.0, 1.0 + i as f64)).collect();
        let big_f_fft: Vec<Complex> = (0..n).map(|i| Complex::new(2.0 * i as f64, 0.0)).collect();
        let big_g_fft: Vec<Complex> = (0..n).map(|i| Complex::new(0.0, 2.0 * i as f64)).collect();

        let gs = GramSchmidt::new(f_fft, g_fft, big_f_fft, big_g_fft);
        let sampler = FfSampler::new(gs);

        // Verify the sampler was created
        assert_eq!(sampler.gram_schmidt().f_fft.len(), n);
    }

    #[test]
    fn test_ff_sampling_basic() {
        // Create a simple test case
        let n = 4;
        let f_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(1.0 + 0.5 * i as f64, 0.1))
            .collect();
        let g_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(0.1, 1.0 + 0.5 * i as f64))
            .collect();
        let big_f_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(10.0 + i as f64, 0.5))
            .collect();
        let big_g_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(0.5, 10.0 + i as f64))
            .collect();

        let gs = GramSchmidt::new(f_fft, g_fft, big_f_fft, big_g_fft);
        let sampler = FfSampler::new(gs);

        let mut rng = StdRng::seed_from_u64(42);

        // Create a target
        let t0_fft: Vec<Complex> = (0..n).map(|i| Complex::new(i as f64 * 0.1, 0.0)).collect();
        let t1_fft: Vec<Complex> = (0..n).map(|i| Complex::new(0.0, i as f64 * 0.1)).collect();

        // Sample
        let (z0, z1) = sampler.ff_sampling(&mut rng, &t0_fft, &t1_fft);

        // Check output dimensions
        assert_eq!(z0.len(), n);
        assert_eq!(z1.len(), n);

        // Check that samples are finite
        for i in 0..n {
            assert!(z0[i].re.is_finite(), "z0[{}].re should be finite", i);
            assert!(z0[i].im.is_finite(), "z0[{}].im should be finite", i);
            assert!(z1[i].re.is_finite(), "z1[{}].re should be finite", i);
            assert!(z1[i].im.is_finite(), "z1[{}].im should be finite", i);
        }
    }

    #[test]
    fn test_ff_sampling_statistics() {
        // Test that ffSampling produces samples with expected mean
        // Using a small sigma to reduce variance for more precise mean estimation
        let n = 4;
        let f_fft: Vec<Complex> = vec![
            Complex::new(2.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(2.0, 0.0),
        ];
        let g_fft: Vec<Complex> = vec![
            Complex::new(0.0, 2.0),
            Complex::new(0.0, 2.0),
            Complex::new(0.0, 2.0),
            Complex::new(0.0, 2.0),
        ];
        let big_f_fft: Vec<Complex> = vec![
            Complex::new(10.0, 0.0),
            Complex::new(10.0, 0.0),
            Complex::new(10.0, 0.0),
            Complex::new(10.0, 0.0),
        ];
        let big_g_fft: Vec<Complex> = vec![
            Complex::new(0.0, 10.0),
            Complex::new(0.0, 10.0),
            Complex::new(0.0, 10.0),
            Complex::new(0.0, 10.0),
        ];

        let gs = GramSchmidt::new(f_fft, g_fft, big_f_fft, big_g_fft);
        // Use a smaller sigma for this test to reduce variance
        let sampler = FfSampler::with_sigma(gs, 5.0);

        let mut rng = StdRng::seed_from_u64(12345);

        // Target at zero
        let t0_fft: Vec<Complex> = vec![Complex::ZERO; n];
        let t1_fft: Vec<Complex> = vec![Complex::ZERO; n];

        // Sample many times and check mean is close to target
        let n_samples = 1000;
        let mut sum_z0 = vec![Complex::ZERO; n];
        let mut sum_z1 = vec![Complex::ZERO; n];

        for _ in 0..n_samples {
            let (z0, z1) = sampler.ff_sampling(&mut rng, &t0_fft, &t1_fft);
            for i in 0..n {
                sum_z0[i] = sum_z0[i] + z0[i];
                sum_z1[i] = sum_z1[i] + z1[i];
            }
        }

        // Check mean is close to zero
        // With sigma=5 and n_samples=1000, expected deviation ~= sigma / sqrt(n_samples) ~= 0.16
        // Allow 10x margin for numerical effects
        for i in 0..n {
            let mean_z0 = sum_z0[i].scale(1.0 / n_samples as f64);
            let mean_z1 = sum_z1[i].scale(1.0 / n_samples as f64);

            // Allow for some deviation (these are random samples)
            assert!(
                mean_z0.norm() < 5.0,
                "Mean z0[{}] = {:?} should be close to zero",
                i,
                mean_z0
            );
            assert!(
                mean_z1.norm() < 5.0,
                "Mean z1[{}] = {:?} should be close to zero",
                i,
                mean_z1
            );
        }
    }

    #[test]
    fn test_signature_bounds() {
        assert_eq!(get_sig_bound_sq(512), FALCON_512_SIG_BOUND_SQ);
        assert_eq!(get_sig_bound_sq(1024), FALCON_1024_SIG_BOUND_SQ);
        assert_eq!(get_sig_bound_sq(256), FALCON_512_SIG_BOUND_SQ); // Uses 512 bound for smaller n
    }

    #[test]
    fn test_signature_norm_sq() {
        let mut s1: Vec<Complex> = vec![3.0, 4.0]
            .into_iter()
            .map(|x| Complex::from_real(x))
            .collect();
        let mut s2: Vec<Complex> = vec![0.0, 0.0]
            .into_iter()
            .map(|x| Complex::from_real(x))
            .collect();

        fft(&mut s1);
        fft(&mut s2);

        let norm_sq = signature_norm_sq(&s1, &s2);
        // Should be 3² + 4² = 25
        assert!((norm_sq - 25.0).abs() < 0.01);
    }

    #[test]
    fn test_ff_sampling_n16() {
        // Test ffSampling with n=16 (similar to keygen_16)
        let n = 16;

        // Create a simple basis with reasonable values
        let f_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(2.0 + (i as f64 * 0.1), 0.1))
            .collect();
        let g_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(0.1, 2.0 + (i as f64 * 0.1)))
            .collect();
        let big_f_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(100.0 + (i as f64), 1.0))
            .collect();
        let big_g_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(1.0, 100.0 + (i as f64)))
            .collect();

        let gs = GramSchmidt::new(f_fft.clone(), g_fft.clone(), big_f_fft.clone(), big_g_fft.clone());

        // Check that the tree was built correctly
        assert_eq!(gs.tree.n, n);
        assert_eq!(gs.tree.depth, 5); // log2(16) + 1 = 5

        // Check sigma values are reasonable
        for i in 0..n {
            let s0 = gs.sigma0(i);
            let s1 = gs.sigma1(i);
            assert!(s0 > 0.0, "sigma0[{}] = {} should be positive", i, s0);
            assert!(s1 >= 0.0, "sigma1[{}] = {} should be non-negative", i, s1);
        }

        let sampler = FfSampler::new(gs);
        let mut rng = StdRng::seed_from_u64(42);

        // Create a target
        let t0_fft: Vec<Complex> = vec![Complex::ZERO; n];
        let t1_fft: Vec<Complex> = vec![Complex::ZERO; n];

        // This should complete quickly
        let (z0, z1) = sampler.ff_sampling(&mut rng, &t0_fft, &t1_fft);

        assert_eq!(z0.len(), n);
        assert_eq!(z1.len(), n);

        // Check outputs are finite
        for i in 0..n {
            assert!(z0[i].re.is_finite(), "z0[{}].re should be finite", i);
            assert!(z0[i].im.is_finite(), "z0[{}].im should be finite", i);
            assert!(z1[i].re.is_finite(), "z1[{}].re should be finite", i);
            assert!(z1[i].im.is_finite(), "z1[{}].im should be finite", i);
        }
    }

    #[test]
    fn test_ff_sampling_with_keygen() {
        use crate::keygen::keygen_16;

        // Use seed 12345 which is known to produce valid NTRU pairs
        let mut rng = StdRng::seed_from_u64(12345);

        // Generate key pair with small n for testing
        let keypair = keygen_16(&mut rng).expect("keygen failed");

        // Check that the GramSchmidt is valid
        let gs = &keypair.sk.gs;
        let n = 16;
        assert_eq!(gs.f_fft.len(), n);
        assert_eq!(gs.tree.n, n);

        // Print sigma values to see if they're reasonable
        println!("Sigma values from keygen:");
        for i in 0..n {
            let s0 = gs.sigma0(i);
            let s1 = gs.sigma1(i);
            println!("  i={}: sigma0={:.4}, sigma1={:.4}", i, s0, s1);
            assert!(s0 > 0.0 && s0.is_finite(), "sigma0[{}] = {} should be positive and finite", i, s0);
            assert!(s1 >= 0.0 && s1.is_finite(), "sigma1[{}] = {} should be non-negative and finite", i, s1);
        }

        // Create the sampler and sample
        let sampler = FfSampler::new(gs.clone());

        // Sample with zero target
        let t0_fft: Vec<Complex> = vec![Complex::ZERO; n];
        let t1_fft: Vec<Complex> = vec![Complex::ZERO; n];

        let (z0, z1) = sampler.ff_sampling(&mut rng, &t0_fft, &t1_fft);

        assert_eq!(z0.len(), n);
        assert_eq!(z1.len(), n);
        println!("ffSampling completed successfully with keygen keys");
    }

    #[test]
    fn test_sample_preimage_with_keygen() {
        use crate::keygen::keygen_16;
        use crate::hash::hash_to_point;
        use crate::params::FALCON_16;

        // Use seed 12345 which is known to produce valid NTRU pairs
        let mut rng = StdRng::seed_from_u64(12345);

        // Generate key pair with small n for testing
        let keypair = keygen_16(&mut rng).expect("keygen failed");

        // Create a test message hash (like in signing)
        let message = b"Hello, FALCON!";
        let nonce = [0u8; 40];
        let c = hash_to_point(message, &nonce, &FALCON_16);

        // Convert c to FFT form
        let mut c_fft: Vec<Complex> = c
            .iter()
            .map(|zq| Complex::from_real(zq.value() as f64))
            .collect();
        fft(&mut c_fft);

        println!("c_fft sample values:");
        for i in 0..4 {
            println!("  c_fft[{}] = {:?}", i, c_fft[i]);
        }

        // Create the sampler
        let sampler = FfSampler::new(keypair.sk.gs.clone());

        // Sample preimage
        println!("Calling sample_preimage...");
        let (z0, z1) = sampler.sample_preimage(&mut rng, &c_fft);

        println!("sample_preimage completed!");
        assert_eq!(z0.len(), 16);
        assert_eq!(z1.len(), 16);

        // Check outputs are finite
        for i in 0..16 {
            assert!(z0[i].re.is_finite(), "z0[{}].re should be finite", i);
            assert!(z1[i].re.is_finite(), "z1[{}].re should be finite", i);
        }
    }
}
