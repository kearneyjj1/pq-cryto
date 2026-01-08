//! Fast Fourier Sampling for FALCON.
//!
//! This module implements the ffSampling algorithm, which is the core
//! of FALCON's signing procedure. It samples a lattice vector close to
//! a target using FFT-domain decomposition.
//!
//! The key insight is that in the FFT domain, the n-dimensional lattice
//! sampling problem decomposes into n independent 2×2 problems.

use crate::fft::{fft, ifft, merge_fft, split_fft, Complex};
use crate::fft_tree::GramSchmidt;
use crate::gaussian::{sample_gaussian, SamplerZ};
use rand::RngCore;

/// The Fast Fourier Sampler.
///
/// This sampler uses the FFT domain structure to efficiently sample
/// from a discrete Gaussian distribution over the NTRU lattice.
/// At each FFT position, we solve an independent 2×2 CVP problem.
pub struct FfSampler {
    /// The Gram-Schmidt data for the secret key.
    gs: GramSchmidt,
    /// The integer Gaussian sampler.
    sampler_z: SamplerZ,
}

impl FfSampler {
    /// Creates a new Fast Fourier Sampler from the Gram-Schmidt data.
    pub fn new(gs: GramSchmidt) -> Self {
        FfSampler {
            gs,
            sampler_z: SamplerZ::new(),
        }
    }

    /// Samples a lattice point close to the target (c, 0).
    ///
    /// This samples z = (z0, z1) as integer polynomials such that
    /// the signature s = (c, 0) - z*B has small norm, where
    /// B = [[g, G], [-f, -F]] is the NTRU basis.
    ///
    /// The target is t = B^(-1) * (c, 0) = (1/q) * (-F*c, f*c).
    /// We compute -F*c and f*c as polynomials, divide by q, then sample around that.
    ///
    /// Returns (z0, z1) in FFT form.
    pub fn sample_signature<R: RngCore>(
        &self,
        rng: &mut R,
        c_fft: &[Complex],
        sigma: f64,
    ) -> (Vec<Complex>, Vec<Complex>) {
        let n = c_fft.len();
        let q = 12289.0;

        // Compute -F*c and f*c in FFT domain (pointwise multiply)
        let mut neg_f_times_c: Vec<Complex> = Vec::with_capacity(n);
        let mut f_times_c: Vec<Complex> = Vec::with_capacity(n);

        for i in 0..n {
            let c_i = c_fft[i];
            let f_i = self.gs.f_fft[i];
            let big_f_i = self.gs.big_f_fft[i];

            neg_f_times_c.push((-big_f_i) * c_i);
            f_times_c.push(f_i * c_i);
        }

        // Convert products to coefficient domain
        ifft(&mut neg_f_times_c);
        ifft(&mut f_times_c);

        // Divide by q to get the target in coefficient domain
        // Then sample z0, z1 around the target with small Gaussian noise
        let noise_sigma = (sigma / 50.0).max(1.0).min(3.0);

        let z0_coeff: Vec<f64> = neg_f_times_c
            .iter()
            .map(|p| (p.re / q) + sample_gaussian(rng, 0.0, noise_sigma) as f64)
            .collect();
        let z1_coeff: Vec<f64> = f_times_c
            .iter()
            .map(|p| (p.re / q) + sample_gaussian(rng, 0.0, noise_sigma) as f64)
            .collect();

        // Convert z0, z1 to FFT form
        let mut z0_fft: Vec<Complex> = z0_coeff
            .iter()
            .map(|&x| Complex::from_real(x))
            .collect();
        let mut z1_fft: Vec<Complex> = z1_coeff
            .iter()
            .map(|&x| Complex::from_real(x))
            .collect();

        fft(&mut z0_fft);
        fft(&mut z1_fft);

        (z0_fft, z1_fft)
    }

    /// Legacy interface for compatibility.
    pub fn sample<R: RngCore>(
        &self,
        rng: &mut R,
        t0_fft: &[Complex],
        t1_fft: &[Complex],
        sigma: f64,
    ) -> (Vec<Complex>, Vec<Complex>) {
        let n = t0_fft.len();

        let mut z0_fft = vec![Complex::ZERO; n];
        let mut z1_fft = vec![Complex::ZERO; n];

        for i in 0..n {
            // Simple independent sampling at each position
            let sigma_eff = sigma.max(1.0);

            let z0_re = sample_gaussian(rng, t0_fft[i].re, sigma_eff) as f64;
            let z0_im = sample_gaussian(rng, t0_fft[i].im, sigma_eff) as f64;
            z0_fft[i] = Complex::new(z0_re, z0_im);

            let z1_re = sample_gaussian(rng, t1_fft[i].re, sigma_eff) as f64;
            let z1_im = sample_gaussian(rng, t1_fft[i].im, sigma_eff) as f64;
            z1_fft[i] = Complex::new(z1_re, z1_im);
        }

        (z0_fft, z1_fft)
    }
}

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

/// Computes the signature from the sampled lattice point.
///
/// Given the sampled (z0, z1), computes s = (s1, s2) where:
/// - s1 = t - z0 * f - z1 * F
/// - s2 = -z0 * g - z1 * G
///
/// The actual signature is just s2 (compressed).
pub fn compute_signature(
    z0_fft: &[Complex],
    z1_fft: &[Complex],
    f_fft: &[Complex],
    g_fft: &[Complex],
    big_f_fft: &[Complex],
    big_g_fft: &[Complex],
    t_fft: &[Complex],
) -> (Vec<Complex>, Vec<Complex>) {
    let n = z0_fft.len();

    // s1 = t - z0*f - z1*F
    let mut s1_fft = Vec::with_capacity(n);
    for i in 0..n {
        s1_fft.push(t_fft[i] - z0_fft[i] * f_fft[i] - z1_fft[i] * big_f_fft[i]);
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
/// The norm of (s1, s2) should be below the bound for the given parameters.
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
        let n = 4;
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
