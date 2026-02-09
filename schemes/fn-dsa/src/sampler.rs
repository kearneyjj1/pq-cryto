//! Fast Fourier Sampling for FALCON.
//!
//! This module implements the ffSampling algorithm, which is the core
//! of FALCON's signing procedure. It samples a lattice vector close to
//! a target using FFT-domain decomposition.
//!
//! The key insight is that in the FFT domain, the n-dimensional lattice
//! sampling problem decomposes into n independent 2×2 problems.
//!
//! ## Basis Convention
//!
//! We use the NTRU secret basis with row-vector convention:
//!
//! ```text
//!   B = [[g, -f], [G, -F]]       where fG - gF = q
//! ```
//!
//! The target for ffSampling is `t = (c, 0)·B⁻¹`:
//!
//! ```text
//!   B⁻¹ = (1/q) · [[-F, f], [-G, g]]
//!   t = (c, 0) · B⁻¹ = (-F·c/q,  f·c/q)
//! ```
//!
//! After sampling `(z0, z1) ≈ t`, the signature is `(s1, s2) = (c, 0) - (z0, z1)·B`:
//!
//! ```text
//!   s2 = z0·f + z1·F
//!   s1 = c - z0·g - z1·G  =  c - s2·h   (since h = g·f⁻¹ mod q)
//! ```

use crate::fft::{merge_fft, split_fft, Complex};
#[cfg(test)]
use crate::fft::ifft;
use crate::fft_tree::GramSchmidt;
use crate::gaussian::SamplerZ;
use rand::RngCore;
use zeroize::Zeroize;

/// The Fast Fourier Sampler (FIPS 206 ffSampling).
///
/// Implements the recursive FFT-domain lattice sampling algorithm.
/// Uses the LDL* tree to decompose the n-dimensional Gaussian sampling
/// into independent 1D problems at the leaves, with conditioning at each level.
pub struct FfSampler<'a> {
    /// Borrowed reference to the Gram-Schmidt / LDL* tree data.
    /// Avoids cloning secret key material on every signing call.
    gs: &'a GramSchmidt,
    /// The FIPS 206 integer Gaussian sampler.
    sampler_z: SamplerZ,
    /// The signing sigma parameter (sigma_sign).
    sigma_sign: f64,
    /// The per-parameter sigma_min.
    sigma_min: f64,
}

impl Drop for FfSampler<'_> {
    fn drop(&mut self) {
        // Clear sigma values as defense-in-depth.
        // GramSchmidt is borrowed, not owned — its Drop runs with SecretKey.
        self.sigma_sign = 0.0;
        self.sigma_min = 0.0;
    }
}

impl<'a> FfSampler<'a> {
    /// Creates a new Fast Fourier Sampler borrowing the Gram-Schmidt data.
    pub fn new(gs: &'a GramSchmidt, sigma_sign: f64, sigma_min: f64) -> Self {
        FfSampler {
            gs,
            sampler_z: SamplerZ::with_sigma_min(sigma_min),
            sigma_sign,
            sigma_min,
        }
    }

    /// Samples a lattice point close to the target `(c, 0)` using ffSampling.
    ///
    /// Computes the target `t = (c, 0)·B⁻¹` where `B = [[g, -f], [G, -F]]`:
    ///
    /// ```text
    ///   t0 = -F·c/q,   t1 = f·c/q     (pointwise in FFT domain)
    /// ```
    ///
    /// Then runs the recursive ffSampling algorithm to sample `(z0, z1) ≈ t`.
    /// See module-level documentation for the basis convention derivation.
    ///
    /// Returns `(z0_fft, z1_fft)` — the sampled lattice point in FFT form.
    pub fn sample_signature<R: RngCore>(
        &self,
        rng: &mut R,
        c_fft: &[Complex],
    ) -> (Vec<Complex>, Vec<Complex>) {
        let n = c_fft.len();
        let q = 12289.0;

        // Target t = (c, 0)·B⁻¹ = (c, 0)·(1/q)·[[-F, f], [-G, g]]
        // First column: t0 = c·(-F)/q = -F·c/q
        // Second column: t1 = c·f/q = f·c/q
        let t0_fft: Vec<Complex> = (0..n)
            .map(|i| (-self.gs.big_f_fft[i] * c_fft[i]).scale(1.0 / q))
            .collect();
        let t1_fft: Vec<Complex> = (0..n)
            .map(|i| (self.gs.f_fft[i] * c_fft[i]).scale(1.0 / q))
            .collect();

        // Run recursive ffSampling
        self.ff_sampling_recursive(rng, &t0_fft, &t1_fft, 0, 0)
    }

    /// Recursive ffSampling per FIPS 206.
    ///
    /// At each level, retrieves the LDL* factor l10 from the tree,
    /// samples z1 via the right subtree, conditions t0, then samples z0
    /// via the left subtree.
    ///
    /// At the leaves (n=1), samples both real and imaginary parts
    /// from the discrete Gaussian using SamplerZ.
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
            // Base case: 2D sampling at a leaf node
            let node = self.gs.tree.get_node(level, pos);

            // sigma_eff = sigma_sign / tree_sigma (the Gram-Schmidt norm)
            // For standard FALCON-512/1024, sigma_eff ≈ sigma_min at all leaves.
            // For degenerate bases (toy n=16), tree_sigma can be near 0, making
            // sigma_eff enormous and ccs = sigma_min/sigma_eff ≈ 0 (very slow).
            // Cap sigma_eff at 4*sigma_min to ensure ccs >= 0.25 (~4 trials avg).
            // This cap never triggers for properly conditioned FALCON keys.
            let max_sigma_eff = 4.0 * self.sigma_min;
            // Branchless sigma computation: .max(1e-300) prevents div-by-zero
            // for degenerate tree nodes, and .clamp() bounds the result.
            // For well-conditioned FALCON keys, node.sigma[i] >> 1e-300 always.
            let sigma_1 = (self.sigma_sign / node.sigma[1].max(1e-300))
                .clamp(self.sigma_min, max_sigma_eff);
            let sigma_0 = (self.sigma_sign / node.sigma[0].max(1e-300))
                .clamp(self.sigma_min, max_sigma_eff);
            let l10_val = node.l10[0];

            // Sample z1 (real part only).
            //
            // With the Falcon tree-ordered FFT, the interleaved split_fft
            // decomposes conjugate pairs at each level.  At the n=1 leaves
            // the target values are purely real, so only the real part is
            // sampled.  This yields 2n total samples (matching the reference
            // implementation) and the correct expected signature norm.
            let z1_re = self.sampler_z.sample(rng, t1[0].re, sigma_1) as f64;
            let z1 = Complex::from_real(z1_re);

            // Condition: t0' = t0 + l10 * (t1 - z1)
            let t0_cond = t0[0] + l10_val * (t1[0] - z1);

            // Sample z0 (real part only, same reasoning as z1)
            let z0_re = self.sampler_z.sample(rng, t0_cond.re, sigma_0) as f64;
            let z0 = Complex::from_real(z0_re);

            return (vec![z0], vec![z1]);
        }

        // Recursive case: split, condition, and recurse
        let node = self.gs.tree.get_node(level, pos);
        // Clone l10 since we need it after the mutable borrow in recursion
        let mut l10: Vec<Complex> = node.l10.clone();

        // 1. Split t1 and recurse on RIGHT subtree
        let (mut t1_even, mut t1_odd) = split_fft(t1);
        let (mut z1_even, mut z1_odd) = self.ff_sampling_recursive(
            rng, &t1_even, &t1_odd, level + 1, 2 * pos + 1,
        );
        t1_even.zeroize();
        t1_odd.zeroize();
        let z1 = merge_fft(&z1_even, &z1_odd);
        z1_even.zeroize();
        z1_odd.zeroize();

        // 2. Condition: t0' = t0 + l10 * (t1 - z1) pointwise
        let mut t0_cond: Vec<Complex> = (0..n)
            .map(|i| t0[i] + l10[i] * (t1[i] - z1[i]))
            .collect();
        l10.zeroize();

        // 3. Split conditioned t0' and recurse on LEFT subtree
        let (mut t0p_even, mut t0p_odd) = split_fft(&t0_cond);
        t0_cond.zeroize();
        let (mut z0_even, mut z0_odd) = self.ff_sampling_recursive(
            rng, &t0p_even, &t0p_odd, level + 1, 2 * pos,
        );
        t0p_even.zeroize();
        t0p_odd.zeroize();
        let z0 = merge_fft(&z0_even, &z0_odd);
        z0_even.zeroize();
        z0_odd.zeroize();

        (z0, z1)
    }
}

/// Simplified sampler for testing.
///
/// This sampler doesn't use the full FFT tree but still produces
/// valid lattice samples (with potentially larger norms).
#[cfg(test)]
pub struct SimpleSampler {
    sampler_z: SamplerZ,
}

#[cfg(test)]
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

#[cfg(test)]
impl Default for SimpleSampler {
    fn default() -> Self {
        Self::new()
    }
}

/// Verifies that the signature has an acceptable norm.
///
/// The norm of (s1, s2) should be below the bound for the given parameters.
#[cfg(test)]
pub fn check_signature_norm(
    s1_fft: &[Complex],
    s2_fft: &[Complex],
    bound_sq: f64,
) -> bool {
    // Convert to coefficient form
    let mut s1 = s1_fft.to_vec();
    let mut s2 = s2_fft.to_vec();
    ifft(&mut s1);
    ifft(&mut s2);

    // Compute squared norm
    let norm_sq: f64 = s1.iter().map(|c| c.re * c.re).sum::<f64>()
        + s2.iter().map(|c| c.re * c.re).sum::<f64>();

    norm_sq <= bound_sq
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fft::fft;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

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
        let _n = 4;
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
}
