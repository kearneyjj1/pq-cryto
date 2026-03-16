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
use crate::gaussian::SamplerZ;
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
    /// We then use recursive ffSampling via the FFT tree to sample around
    /// that target.
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

        // Compute target: t = B^(-1) * (c, 0) = (1/q) * (-F*c, f*c)
        let mut t0: Vec<Complex> = Vec::with_capacity(n);
        let mut t1: Vec<Complex> = Vec::with_capacity(n);

        for i in 0..n {
            let c_i = c_fft[i];
            t0.push((-self.gs.big_f_fft[i] * c_i).scale(1.0 / q));
            t1.push((self.gs.f_fft[i] * c_i).scale(1.0 / q));
        }

        // Recursive ffSampling using the FFT tree
        self.ff_sample_recursive(rng, &t0, &t1, 0, 0)
    }

    /// Recursive ffSampling using the FFT tree (FALCON spec Section 3.8.2).
    ///
    /// At each level, the problem is split into two half-size subproblems
    /// using the FFT split/merge operations. At the leaves (n=1), a 2x2
    /// Gram-Schmidt sampling is performed.
    ///
    /// The tree stores the Gram-Schmidt data at each node:
    /// - Leaves: sigma[0] = ||b1||, sigma[1] = ||b2*|| (these ARE the sampling widths)
    /// - Internal: sigma[i] = ||b1[i]||, gs_coeffs[i] = <(F,G),(f,g)*>[i]
    fn ff_sample_recursive<R: RngCore>(
        &self,
        rng: &mut R,
        t0: &[Complex],
        t1: &[Complex],
        level: usize,
        pos: usize,
    ) -> (Vec<Complex>, Vec<Complex>) {
        let n = t0.len();
        let node = self.gs.tree.get_node(level, pos);

        if n == 1 {
            // TEMPORARY DEBUG: pure rounding, no GS adjustment
            let z1_val = t1[0].re.round();
            let z0_val = t0[0].re.round();

            return (
                vec![Complex::from_real(z0_val)],
                vec![Complex::from_real(z1_val)],
            );
        }

        // Recursive case

        // Split t1 and recurse on RIGHT child
        let (t1_lo, t1_hi) = split_fft(t1);
        let (z1_lo, z1_hi) = self.ff_sample_recursive(
            rng, &t1_lo, &t1_hi, level + 1, 2 * pos + 1
        );
        let z1_full = merge_fft(&z1_lo, &z1_hi);

        // Adjust t0 using this node's Gram-Schmidt coefficients
        // l[i] = gs_coeffs[i] / sigma[i]^2
        // t0_adj[i] = t0[i] + l[i] * (t1[i] - z1_full[i])
        let mut t0_adj = Vec::with_capacity(n);
        for i in 0..n {
            let sigma_sq = node.sigma[i] * node.sigma[i];
            let l_i = if sigma_sq > 1e-30 {
                node.gs_coeffs[i].scale(1.0 / sigma_sq)
            } else {
                Complex::ZERO
            };
            let diff = t1[i] - z1_full[i];
            t0_adj.push(t0[i] + l_i * diff);
        }

        // Split adjusted t0 and recurse on LEFT child
        let (t0_lo, t0_hi) = split_fft(&t0_adj);
        let (z0_lo, z0_hi) = self.ff_sample_recursive(
            rng, &t0_lo, &t0_hi, level + 1, 2 * pos
        );
        let z0_full = merge_fft(&z0_lo, &z0_hi);

        (z0_full, z1_full)
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

#[cfg(test)]
mod ff_sampling_tests {
    use crate::keygen::keygen_16;
    use crate::params::FALCON_16;
    use crate::sign::sign;
    use crate::verify::verify;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_split_fft_produces_real_at_leaves() {
        use crate::fft::{fft, split_fft, merge_fft, Complex};

        // Start with a real polynomial
        let n = 16;
        let mut poly: Vec<Complex> = (0..n).map(|i| Complex::from_real((i as f64) - 8.0)).collect();
        fft(&mut poly);

        // Check imaginary parts at each split level
        fn check_splits(v: &[Complex], depth: usize, label: &str) {
            if v.len() == 1 {
                eprintln!("  depth={} {}: re={:.6}, im={:.6}", depth, label, v[0].re, v[0].im);
                return;
            }
            let (lo, hi) = split_fft(v);
            check_splits(&lo, depth + 1, &format!("{}/lo", label));
            check_splits(&hi, depth + 1, &format!("{}/hi", label));
        }

        eprintln!("Splitting FFT of real polynomial:");
        check_splits(&poly, 0, "root");
    }

    #[test]
    fn test_ff_sampling_debug_norms() {
        use crate::fft::{fft, ifft, Complex};
        use crate::hash::{generate_nonce, hash_to_point};
        use crate::sampler::FfSampler;
        use crate::params::Q;

        let mut rng = StdRng::seed_from_u64(42);
        let keypair = keygen_16(&mut rng).expect("keygen_16 failed");

        // Print tree sigma values at leaves
        let tree = &keypair.sk.gs.tree;
        let n = keypair.sk.params.n;
        let depth = tree.depth;
        let leaf_level = depth - 1;
        eprintln!("n={}, depth={}, leaf_level={}", n, depth, leaf_level);
        for pos in 0..n {
            let node = tree.get_node(leaf_level, pos);
            eprintln!("  leaf[{}]: sigma={:?}", pos, node.sigma);
        }
        eprintln!("sig_bound_sq = {}", keypair.sk.params.sig_bound_sq);
        eprintln!("sigma = {}", keypair.sk.params.sigma);

        // Now try to sign manually and print the norm
        let sk = &keypair.sk;
        let ff_sampler = FfSampler::new(sk.gs.clone());

        let mut f_fft: Vec<Complex> = sk.f.iter().map(|&x| Complex::from_real(x as f64)).collect();
        let mut g_fft: Vec<Complex> = sk.g.iter().map(|&x| Complex::from_real(x as f64)).collect();
        let mut big_f_fft: Vec<Complex> = sk.big_f.iter().map(|&x| Complex::from_real(x as f64)).collect();
        let mut big_g_fft: Vec<Complex> = sk.big_g.iter().map(|&x| Complex::from_real(x as f64)).collect();
        fft(&mut f_fft);
        fft(&mut g_fft);
        fft(&mut big_f_fft);
        fft(&mut big_g_fft);

        let nonce = generate_nonce(&mut rng);
        let c = hash_to_point(b"test", &nonce, &sk.params);
        let mut c_fft: Vec<Complex> = c.iter().map(|zq| Complex::from_real(zq.value() as f64)).collect();
        fft(&mut c_fft);

        // Print the target values
        let q = 12289.0;
        let mut t0_dbg: Vec<Complex> = Vec::with_capacity(n);
        let mut t1_dbg: Vec<Complex> = Vec::with_capacity(n);
        for i in 0..n {
            t0_dbg.push((-sk.gs.big_f_fft[i] * c_fft[i]).scale(1.0 / q));
            t1_dbg.push((sk.gs.f_fft[i] * c_fft[i]).scale(1.0 / q));
        }
        let mut t0_coeff = t0_dbg.clone();
        let mut t1_coeff = t1_dbg.clone();
        ifft(&mut t0_coeff);
        ifft(&mut t1_coeff);
        eprintln!("t0 (coeff domain): {:?}", t0_coeff.iter().map(|c| c.re).collect::<Vec<_>>());
        eprintln!("t1 (coeff domain): {:?}", t1_coeff.iter().map(|c| c.re).collect::<Vec<_>>());

        let (z0_fft, z1_fft) = ff_sampler.sample_signature(&mut rng, &c_fft, sk.params.sigma);

        // Compute s2 = z0*f + z1*F
        let mut s2_fft: Vec<Complex> = Vec::with_capacity(n);
        for i in 0..n {
            s2_fft.push(z0_fft[i] * f_fft[i] + z1_fft[i] * big_f_fft[i]);
        }
        ifft(&mut s2_fft);
        let s2: Vec<i16> = s2_fft.iter()
            .map(|c| {
                let val = c.re.round() as i32;
                let reduced = ((val % (Q as i32)) + (Q as i32)) % (Q as i32);
                if reduced > (Q as i32) / 2 { (reduced - (Q as i32)) as i16 } else { reduced as i16 }
            })
            .collect();

        // Convert z0, z1 to coefficient domain and print
        let mut z0_coeff = z0_fft.clone();
        let mut z1_coeff = z1_fft.clone();
        ifft(&mut z0_coeff);
        ifft(&mut z1_coeff);
        eprintln!("z0 coefficients: {:?}", z0_coeff.iter().map(|c| c.re).collect::<Vec<_>>());
        eprintln!("z1 coefficients: {:?}", z1_coeff.iter().map(|c| c.re).collect::<Vec<_>>());

        let s2_norm_sq: i64 = s2.iter().map(|&x| (x as i64) * (x as i64)).sum();
        eprintln!("s2 coefficients: {:?}", s2);
        eprintln!("s2_norm_sq = {}", s2_norm_sq);

        // Also compute s1 for full norm
        use crate::poly::Poly;
        let c_poly = Poly::from_zq(c.clone());
        let s2_poly = Poly::from_i16(&s2);
        let h_poly = Poly::from_i16(&sk.h);
        let s2h = s2_poly.mul(&h_poly);
        let s1_poly = c_poly.sub(&s2h);
        let s1_norm_sq = s1_poly.norm_sq();
        let total = s1_norm_sq + s2_norm_sq;
        eprintln!("s1_norm_sq = {}, total = {}, bound = {}", s1_norm_sq, total, sk.params.sig_bound_sq);
    }

    #[test]
    fn test_ff_sampling_falcon_16() {
        let mut rng = StdRng::seed_from_u64(42);

        // Generate key and sign
        let keypair = keygen_16(&mut rng).expect("keygen_16 failed");
        let msg = b"test message";
        let sig = sign(&mut rng, &keypair.sk, msg).expect("signing failed");

        // Verify
        assert!(verify(&keypair.pk, msg, &sig).is_ok(), "Signature should verify");
    }

    #[test]
    fn test_ff_sampling_falcon_16_multiple() {
        let mut rng = StdRng::seed_from_u64(123);

        let keypair = keygen_16(&mut rng).expect("keygen_16 failed");

        // Sign multiple messages and check all verify
        for i in 0..5 {
            let msg = format!("message {}", i);
            let sig = sign(&mut rng, &keypair.sk, msg.as_bytes()).expect("signing failed");
            assert!(
                verify(&keypair.pk, msg.as_bytes(), &sig).is_ok(),
                "Signature {} should verify", i
            );
        }
    }
}
