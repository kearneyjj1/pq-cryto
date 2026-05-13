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
use zeroize::Zeroize;

/// The Fast Fourier Sampler.
///
/// Samples a lattice vector close to a target using the Gram-Schmidt
/// structure of the secret key basis. Converts the Babai target to
/// the coefficient domain and samples each coefficient independently
/// from a discrete Gaussian centered at the target value.
pub struct FfSampler<'a> {
    /// Borrowed reference to the Gram-Schmidt / LDL* tree data.
    gs: &'a GramSchmidt,
    /// The integer Gaussian sampler.
    sampler_z: SamplerZ,
    /// The signing sigma parameter.
    sigma_sign: f64,
    /// The minimum sigma (sigma_min) for the leaf sampler. Used to assert
    /// that the per-leaf effective sigma stays within the FIPS 206 validity
    /// range; a leaf sigma below sigma_min indicates a mis-scaled tree.
    sigma_min: f64,
}

impl<'a> FfSampler<'a> {
    /// Creates a new Fast Fourier Sampler borrowing the Gram-Schmidt data.
    ///
    /// `sigma_sign` is the global FALCON signing sigma (`Params::sigma`).
    /// `sigma_min` is the per-leaf minimum sigma (`Params::sigma_min`); it
    /// is used both by the underlying [`SamplerZ`] and as a floor on the
    /// per-leaf effective sigma `sigma_sign / sqrt(d_i)`.
    ///
    /// # Spec note on the leaf-sigma floor
    ///
    /// FIPS 206 Algorithm 11 specifies `sigma' = sigma_sign / sqrt(d_i)`
    /// without an explicit floor; instead the spec requires the key
    /// generator to reject keys whose LDL\* diagonals would push any
    /// leaf below `sigma_min`. Our keygen's `check_gram_schmidt_quality`
    /// enforces this for production parameter sets (FALCON-512/1024).
    ///
    /// For toy parameter sets that intentionally skip the quality check
    /// (e.g. FALCON-16, used in tests), we clamp at `sigma_min` so the
    /// sampler still terminates instead of falling into a degenerate
    /// rejection loop. **The clamped distribution is not FIPS-compliant**;
    /// it is a robustness measure for non-spec parameter sets only. If
    /// the clamp ever engages on a production parameter set, that is a
    /// `keygen` quality-check bug, not a sampler bug.
    pub fn new(gs: &'a GramSchmidt, sigma_sign: f64, sigma_min: f64) -> Self {
        FfSampler {
            gs,
            sampler_z: SamplerZ::with_sigma_min(sigma_min),
            sigma_sign,
            sigma_min,
        }
    }

    /// Samples a lattice point close to the target (c, 0).
    ///
    /// Computes the Babai target t = B^(-1) * (c, 0) = (1/q) * (-F*c, f*c)
    /// in the FFT domain, then uses recursive ffSampling via the LDL* tree
    /// to produce a short lattice vector.
    ///
    /// Returns (z0_fft, z1_fft) in FFT form.
    pub fn sample_signature<R: RngCore>(
        &self,
        rng: &mut R,
        c_fft: &[Complex],
    ) -> (Vec<Complex>, Vec<Complex>) {
        let n = c_fft.len();
        let q = 12289.0;

        // Compute target in FFT domain: t0 = -F*c/q, t1 = f*c/q
        let mut t0: Vec<Complex> = Vec::with_capacity(n);
        let mut t1: Vec<Complex> = Vec::with_capacity(n);
        for i in 0..n {
            let c_i = c_fft[i];
            t0.push((-self.gs.big_f_fft[i] * c_i).scale(1.0 / q));
            t1.push((self.gs.f_fft[i] * c_i).scale(1.0 / q));
        }

        // Recursive ffSampling using the LDL* tree. After the call, the
        // secret-derived target vectors are scrubbed; the returned z0/z1
        // are the caller's responsibility.
        let result = self.ff_sample_recursive(rng, &t0, &t1, 0, 0);
        t0.zeroize();
        t1.zeroize();
        result
    }

    /// Recursive ffSampling (FIPS 206 Algorithm 11).
    ///
    /// The LDL* tree stores at each node:
    /// - sigma[i] = sqrt(d0[i]) where d0 is from the LDL* decomposition G = L*D*L*
    /// - gs_coeffs[i] = l10[i], the lower-triangular coefficient
    ///
    /// At leaves (n=1):
    /// - sigma = [sqrt(d0), sqrt(d1)]
    /// - gs_coeffs = [l10]
    /// - Sampling width = sigma_sign * sigma[i]
    ///
    /// At internal nodes (n>1):
    /// - sigma[i] = sqrt(d0[i]) per FFT position
    /// - gs_coeffs[i] = l10[i] per FFT position
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
            // Base case: scalar LDL* sampling per FIPS 206 Algorithm 11.
            // node.sigma[i] stores sqrt(d_i) from the LDL* decomposition
            // G = L*D*L*. The leaf sampling width is sigma_sign / sqrt(d_i),
            // i.e. sigma_sign / node.sigma[i]. For FIPS-compliant parameter
            // sets (FALCON-512/1024), keygen's basis-quality check ensures
            // this is always >= sigma_min. Toy parameter sets used only in
            // tests (e.g. FALCON-16) cannot enforce basis quality, so we
            // clamp at sigma_min for sampler robustness. On production
            // parameters the clamp never engages.
            //
            // Hermitian-symmetry sanity check: Falcon's tree-ordered FFT
            // preserves real-valued leaves for real-input polynomials, so
            // both t0[0] and t1[0] should have negligible imaginary part
            // at the deepest recursion. A meaningful im here would signal
            // either a wrong FFT root ordering or an integer-rounding step
            // that broke symmetry mid-recursion.
            //
            // The theoretical FP error bound for n=512 is on the order of
            // log2(n) * eps_machine * max|t| ≈ 9 * 2.2e-16 * O(100) ≈ 2e-13.
            // We assert with a 1e-9 tolerance: well above the theoretical
            // bound but tight enough to catch any genuine symmetry break.
            // Using `assert!` rather than `debug_assert!` so the check is
            // active in release builds; the two scalar comparisons per
            // signing call are negligible cost.
            assert!(
                t0[0].im.abs() < 1e-9,
                "leaf t0 has non-real value im={}", t0[0].im
            );
            assert!(
                t1[0].im.abs() < 1e-9,
                "leaf t1 has non-real value im={}", t1[0].im
            );
            let raw_sigma_0 = self.sigma_sign / node.sigma[0];
            let raw_sigma_1 = self.sigma_sign / node.sigma[1];
            // Runtime canary: on FIPS-compliant parameter sets (n >= 256)
            // keygen's basis-quality check ensures both raw sigmas already
            // satisfy `>= sigma_min`, so the clamp below is a no-op. If
            // this debug_assert ever fires on FALCON-512/1024 it indicates
            // a keygen quality-check bug or a regression in the LDL* tree.
            // Toy parameter sets (n < 256) deliberately skip the check and
            // may engage the clamp; the assert is gated accordingly.
            debug_assert!(
                self.gs.tree.n < 256 || raw_sigma_0 >= self.sigma_min,
                "leaf sigma_0 {} below sigma_min {} on production parameter set",
                raw_sigma_0, self.sigma_min
            );
            debug_assert!(
                self.gs.tree.n < 256 || raw_sigma_1 >= self.sigma_min,
                "leaf sigma_1 {} below sigma_min {} on production parameter set",
                raw_sigma_1, self.sigma_min
            );
            let sigma_0 = raw_sigma_0.max(self.sigma_min);
            let sigma_1 = raw_sigma_1.max(self.sigma_min);
            let l10 = node.l10[0];

            // Sample z1 first (second basis vector direction)
            let z1_val = self.sampler_z.sample(rng, t1[0].re, sigma_1) as f64;

            // Adjust t0: t0_adj = t0 + l10 * (t1 - z1)
            let t0_adj = t0[0].re + l10.re * (t1[0].re - z1_val);

            // Sample z0 from adjusted target
            let z0_val = self.sampler_z.sample(rng, t0_adj, sigma_0) as f64;

            return (
                vec![Complex::from_real(z0_val)],
                vec![Complex::from_real(z1_val)],
            );
        }

        // Recursive case: split, recurse, adjust, recurse, merge.
        //
        // Every intermediate vector below carries secret-derived state
        // (FFT of a function of f, g, F, G, c). We scrub each one once it
        // is no longer used. Per recursive call this is 10 heap vectors;
        // for n=512 the recursion has 1023 internal nodes, so a single
        // signing call would otherwise leave ~10k secret-bearing vectors
        // on the heap waiting for allocator reuse.

        // Split t1 and recurse on RIGHT child (d1 subtree).
        let (mut t1_lo, mut t1_hi) = split_fft(t1);
        let (mut z1_lo, mut z1_hi) = self.ff_sample_recursive(
            rng, &t1_lo, &t1_hi, level + 1, 2 * pos + 1
        );
        let z1_full = merge_fft(&z1_lo, &z1_hi);
        t1_lo.zeroize();
        t1_hi.zeroize();
        z1_lo.zeroize();
        z1_hi.zeroize();

        // Adjust t0 using this node's l10 coefficients:
        // t0_adj[i] = t0[i] + l10[i] * (t1[i] - z1[i])
        let mut t0_adj = Vec::with_capacity(n);
        for i in 0..n {
            let diff = t1[i] - z1_full[i];
            t0_adj.push(t0[i] + node.l10[i] * diff);
        }

        // Split adjusted t0 and recurse on LEFT child (d0 subtree).
        let (mut t0_lo, mut t0_hi) = split_fft(&t0_adj);
        t0_adj.zeroize();
        let (mut z0_lo, mut z0_hi) = self.ff_sample_recursive(
            rng, &t0_lo, &t0_hi, level + 1, 2 * pos
        );
        let z0_full = merge_fft(&z0_lo, &z0_hi);
        t0_lo.zeroize();
        t0_hi.zeroize();
        z0_lo.zeroize();
        z0_hi.zeroize();

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

#[cfg(test)]
#[allow(deprecated)] // The debug-norms test uses `Poly::mul` for diagnostic comparison.
mod ff_sampling_tests {
    use crate::keygen::keygen_16;
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
        let ff_sampler = FfSampler::new(&sk.gs, sk.params.sigma, sk.params.sigma_min);

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

        let (z0_fft, z1_fft) = ff_sampler.sample_signature(&mut rng, &c_fft);

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
