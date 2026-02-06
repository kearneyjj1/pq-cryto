//! Discrete Gaussian sampling for FALCON.
//!
//! This module implements the discrete Gaussian sampler used in FALCON's
//! signing algorithm. The sampler produces integers from a discrete Gaussian
//! distribution centered at a given mean with a given standard deviation.
//!
//! # Security Notes
//!
//! This implementation uses CDT-based sampling with constant-time table lookups
//! for the base sampler. However, some operations (especially for large sigma)
//! may still have timing variations. For maximum security, hardware-specific
//! constant-time implementations should be used.

use rand::{Rng, RngCore};

/// Base standard deviation for the half-Gaussian CDT sampler.
/// This is sigma_0 = 1.8205 from the FALCON / FIPS 206 specification.
pub const SIGMA_0: f64 = 1.8205;

/// Minimum sigma for leaf nodes in the LDL* tree (FALCON-512).
/// sigma_min = 1.277833697 per FIPS 206.
pub const SIGMA_MIN: f64 = 1.277833697;

/// Precomputed 1 / (2 * sigma_0^2) for the BerExp rejection test.
const INV_2SIGMA0_SQ: f64 = 0.150_868_809_571_805_55;

/// Maximum deviation from the mean (in units of sigma) that we consider.
/// Values beyond this are essentially zero probability.
const MAX_SIGMA_MULT: f64 = 10.0;

/// CDT table precision: number of bits of precision.
const CDT_BITS: usize = 72;

// ============================================================================
// Base Gaussian Sampler
// ============================================================================

/// Samples from a discrete Gaussian distribution N(0, sigma^2) over Z.
///
/// Uses the cumulative distribution table (CDT) method for efficiency.
/// Returns an integer z sampled according to exp(-z^2 / (2*sigma^2)).
pub fn sample_z_gaussian<R: RngCore>(rng: &mut R, sigma: f64) -> i64 {
    if sigma <= 0.0 {
        return 0;
    }

    // For small sigma, use rejection sampling
    // For larger sigma, use the half-Gaussian CDT method

    // Sample sign: +1 or -1
    let sign: i64 = if rng.gen_bool(0.5) { 1 } else { -1 };

    // Sample |z| from the half-Gaussian
    let abs_z = sample_half_gaussian(rng, sigma);

    // Handle z=0 specially: it should appear with probability 1/(1 + 2*sum)
    // For now, use the sign convention
    if abs_z == 0 {
        0
    } else {
        sign * abs_z
    }
}

/// Samples |z| from a half-Gaussian (z >= 0).
fn sample_half_gaussian<R: RngCore>(rng: &mut R, sigma: f64) -> i64 {
    // Use rejection sampling for the half-Gaussian
    let max_z = (sigma * MAX_SIGMA_MULT).ceil() as i64;

    loop {
        // Sample uniformly from [0, max_z]
        let z = rng.gen_range(0..=max_z);

        // Accept with probability proportional to exp(-z^2 / (2*sigma^2))
        let prob = (-((z * z) as f64) / (2.0 * sigma * sigma)).exp();

        // Use a uniform random number to decide
        let u: f64 = rng.gen();
        if u < prob {
            return z;
        }
    }
}

/// Samples from a discrete Gaussian N(mu, sigma^2) centered at mu.
///
/// This samples z such that the probability is proportional to exp(-(z-mu)^2 / (2*sigma^2)).
pub fn sample_gaussian<R: RngCore>(rng: &mut R, mu: f64, sigma: f64) -> i64 {
    if sigma <= 0.0 {
        return mu.round() as i64;
    }

    // Sample z from N(0, sigma^2), then shift by mu and round
    // This is an approximation; proper FALCON uses a more sophisticated approach

    // For FALCON, we use rejection sampling around mu
    let max_offset = (sigma * MAX_SIGMA_MULT).ceil() as i64;
    let mu_floor = mu.floor() as i64;

    loop {
        // Sample uniformly around mu
        let z = mu_floor + rng.gen_range(-max_offset..=max_offset);

        // Accept with probability proportional to exp(-(z-mu)^2 / (2*sigma^2))
        let diff = z as f64 - mu;
        let prob = (-(diff * diff) / (2.0 * sigma * sigma)).exp();

        let u: f64 = rng.gen();
        if u < prob {
            return z;
        }
    }
}

// ============================================================================
// Bernoullí Sampler
// ============================================================================

/// Samples a Bernoulli variable: returns true with probability exp(-x) for x >= 0.
pub fn sample_bernoulli_exp<R: RngCore>(rng: &mut R, x: f64) -> bool {
    if x <= 0.0 {
        return true;
    }
    if x >= 1.0 {
        // For x >= 1, recurse: P(accept) = exp(-x) = exp(-1) * exp(-(x-1))
        // So accept if both exp(-1) and exp(-(x-1)) accept
        if !sample_bernoulli_exp_minus_one(rng) {
            return false;
        }
        return sample_bernoulli_exp(rng, x - 1.0);
    }

    // For 0 < x < 1, use direct comparison
    // This is simpler and more reliable for the use case
    let prob = (-x).exp();
    rng.gen::<f64>() < prob
}

/// Samples Bernoulli with probability exp(-1) ≈ 0.368.
fn sample_bernoulli_exp_minus_one<R: RngCore>(rng: &mut R) -> bool {
    // exp(-1) ≈ 0.36787944117144233
    // Use a precomputed table or direct comparison
    const EXP_MINUS_1: f64 = 0.36787944117144233;
    rng.gen::<f64>() < EXP_MINUS_1
}

/// Samples Bernoulli with probability cosh(x)^(-1) for computing rejection.
fn sample_bernoulli_cosh<R: RngCore>(rng: &mut R, x: f64) -> bool {
    // P(accept) = 1/cosh(x) = 2/(e^x + e^(-x))
    let prob = 1.0 / x.cosh();
    rng.gen::<f64>() < prob
}

// ============================================================================
// BaseSampler for FALCON
// ============================================================================

/// The base sampler for FALCON's signature scheme.
///
/// This samples from a discrete Gaussian with standard deviation sigma_min = 1.277833,
/// which is the minimum Gaussian width used in the FALCON sampler.
pub struct BaseSampler {
    /// Precomputed CDT for the base Gaussian.
    cdt: Vec<u64>,
    /// Sigma for the base sampler.
    sigma: f64,
}

impl BaseSampler {
    /// Creates a new base sampler at sigma_0 = 1.8205 (FIPS 206).
    pub fn new() -> Self {
        let sigma = SIGMA_0;
        let cdt = Self::compute_cdt(sigma);
        BaseSampler { cdt, sigma }
    }

    /// Creates a base sampler with a custom sigma.
    pub fn with_sigma(sigma: f64) -> Self {
        let cdt = Self::compute_cdt(sigma);
        BaseSampler { cdt, sigma }
    }

    /// Computes the CDT for a given sigma.
    fn compute_cdt(sigma: f64) -> Vec<u64> {
        // CDT[i] = floor(2^63 * sum_{j=0}^{i} P(z=j))
        // where P(z=j) = exp(-j^2 / (2*sigma^2)) / Z
        // and Z = sum_{j=-inf}^{inf} exp(-j^2 / (2*sigma^2))

        // First compute the normalizing constant (for half-Gaussian)
        let max_z = (sigma * MAX_SIGMA_MULT).ceil() as usize;
        let mut probs = Vec::with_capacity(max_z + 1);

        for z in 0..=max_z {
            let zf = z as f64;
            probs.push((-(zf * zf) / (2.0 * sigma * sigma)).exp());
        }

        // Normalize (for half-Gaussian, z=0 has weight 1, z>0 has weight 2)
        let z: f64 = probs[0] + 2.0 * probs[1..].iter().sum::<f64>();

        // Build CDT
        let mut cdt = Vec::with_capacity(max_z + 1);
        let mut cumsum = 0.0;
        let scale = (1u64 << 63) as f64;

        for (i, &p) in probs.iter().enumerate() {
            let weight = if i == 0 { 1.0 } else { 2.0 };
            cumsum += weight * p / z;
            cdt.push((cumsum * scale) as u64);
        }

        cdt
    }

    /// Samples from the base discrete Gaussian.
    ///
    /// # Security
    ///
    /// This uses a constant-time linear scan through the CDT table,
    /// avoiding data-dependent branches that could leak information.
    /// All comparisons use arithmetic masking instead of conditional branches.
    pub fn sample<R: RngCore>(&self, rng: &mut R) -> i64 {
        // Sample from the half-Gaussian using CDT
        let u: u64 = rng.gen();
        let half_u = u >> 1; // Use top 63 bits

        // Constant-time linear scan through CDT
        // Count how many thresholds half_u exceeds
        let mut abs_z: i64 = 0;
        for &threshold in &self.cdt {
            // Constant-time comparison: (half_u >= threshold) without branching
            // If half_u >= threshold, then (half_u.wrapping_sub(threshold)) has high bit 0
            // If half_u < threshold, then wrapping_sub produces a large number with high bit 1
            let diff = half_u.wrapping_sub(threshold);
            let exceeded = ((!(diff >> 63)) & 1) as i64; // 1 if half_u >= threshold, 0 otherwise
            abs_z += exceeded;
        }

        // Sample sign using the LSB of u (constant-time)
        let sign_bit = (u & 1) as i64;
        let sign = 1 - 2 * sign_bit; // 1 if bit=0, -1 if bit=1

        // Apply sign, but return 0 if abs_z == 0
        // Use constant-time selection: result = sign * abs_z (which is 0 when abs_z is 0)
        sign * abs_z
    }

    /// Samples z0 >= 0 from the half-Gaussian at sigma_0.
    ///
    /// This is used as the proposal distribution in the FIPS 206 SamplerZ.
    pub fn sample_half<R: RngCore>(&self, rng: &mut R) -> i64 {
        let u: u64 = rng.gen();
        let half_u = u >> 1; // 63-bit precision for CDT

        // Constant-time linear scan through CDT
        let mut z0: i64 = 0;
        for &threshold in &self.cdt {
            let diff = half_u.wrapping_sub(threshold);
            let exceeded = ((!(diff >> 63)) & 1) as i64;
            z0 += exceeded;
        }

        z0
    }
}

impl Default for BaseSampler {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// BerExp: Bernoulli exponential acceptance test
// ============================================================================

/// Accepts with probability ccs * exp(-x).
///
/// For x >= 0: tests Bernoulli(exp(-x)) AND Bernoulli(ccs) independently.
/// For x < 0: exp(-x) > 1, so computes the product directly and caps at 1.
fn ber_exp<R: RngCore>(rng: &mut R, x: f64, ccs: f64) -> bool {
    if x < 0.0 {
        // exp(-x) > 1; compute P = min(1, ccs * exp(-x)) directly
        let p = (ccs * (-x).exp()).min(1.0);
        return rng.gen::<f64>() < p;
    }

    // x >= 0: P = ccs * exp(-x). Test each factor independently.
    if !sample_bernoulli_exp(rng, x) {
        return false;
    }
    rng.gen::<f64>() < ccs
}

// ============================================================================
// SamplerZ: FIPS 206 Algorithm 12 — discrete Gaussian sampler
// ============================================================================

/// Samples an integer from a discrete Gaussian distribution N(mu, sigma^2).
///
/// Implements FIPS 206 Algorithm 12 (SamplerZ) using the half-Gaussian
/// base sampler at sigma_0, sign randomization, and BerExp rejection.
pub struct SamplerZ {
    base: BaseSampler,
    /// The sigma_min parameter for this parameter set.
    sigma_min: f64,
}

impl SamplerZ {
    /// Creates a new SamplerZ with the FIPS 206 base sampler at sigma_0.
    /// Uses the default FALCON-512 sigma_min.
    pub fn new() -> Self {
        SamplerZ {
            base: BaseSampler::new(),
            sigma_min: SIGMA_MIN,
        }
    }

    /// Creates a new SamplerZ with a specific sigma_min value.
    /// Use this for FALCON-1024 which has a different sigma_min.
    pub fn with_sigma_min(sigma_min: f64) -> Self {
        SamplerZ {
            base: BaseSampler::new(),
            sigma_min,
        }
    }

    /// Samples z from the discrete Gaussian N(mu, sigma^2) per FIPS 206 Algorithm 12.
    ///
    /// For sigma >= sigma_min, uses the FIPS 206 algorithm. For sigma < sigma_min
    /// (which can occur with toy test parameters), falls back to simple rejection
    /// sampling.
    pub fn sample<R: RngCore>(&self, rng: &mut R, mu: f64, sigma: f64) -> i64 {
        // Fallback for sigma below sigma_min (toy parameters with poorly conditioned basis)
        if sigma < self.sigma_min {
            return sample_gaussian(rng, mu, sigma.max(0.01));
        }

        let s = mu.floor() as i64;
        let r = mu - (s as f64);
        let dss = 1.0 / (2.0 * sigma * sigma);
        let ccs = self.sigma_min / sigma;

        loop {
            // 1. Sample z0 >= 0 from half-Gaussian at sigma_0
            let z0 = self.base.sample_half(rng);

            // 2. Random sign bit
            let b: i64 = if rng.gen_bool(0.5) { 1 } else { 0 };

            // 3. z = b + (2b - 1) * z0
            //    b=1 → z = 1 + z0 (positive side)
            //    b=0 → z = -z0    (negative side)
            let z = b + (2 * b - 1) * z0;

            // 4. Reject (z=0, b=0) to avoid double-counting zero
            if z == 0 && b == 0 {
                continue;
            }

            // 5. Compute exponent for rejection test
            //    x = (z - r)^2 / (2*sigma^2) - z0^2 / (2*sigma_0^2)
            let zr = (z as f64) - r;
            let x = zr * zr * dss - (z0 as f64) * (z0 as f64) * INV_2SIGMA0_SQ;

            // 6. Accept with probability ccs * exp(-x)
            if ber_exp(rng, x, ccs) {
                return s + z;
            }
        }
    }
}

impl Default for SamplerZ {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_sample_z_gaussian_basic() {
        let mut rng = StdRng::seed_from_u64(42);
        let sigma = 2.0;

        // Sample many values and check basic properties
        let mut sum: i64 = 0;
        let mut sum_sq: i64 = 0;
        let n = 10000;

        for _ in 0..n {
            let z = sample_z_gaussian(&mut rng, sigma);
            sum += z;
            sum_sq += z * z;
        }

        // Mean should be close to 0
        let mean = (sum as f64) / (n as f64);
        assert!(mean.abs() < 0.5, "Mean {} should be close to 0", mean);

        // Variance should be close to sigma^2
        let variance = (sum_sq as f64) / (n as f64) - mean * mean;
        let expected_var = sigma * sigma;
        assert!(
            (variance - expected_var).abs() < 1.0,
            "Variance {} should be close to {}",
            variance,
            expected_var
        );
    }

    #[test]
    fn test_sample_gaussian_mean() {
        let mut rng = StdRng::seed_from_u64(123);
        let mu = 5.0;
        let sigma = 2.0;

        let mut sum: i64 = 0;
        let n = 5000;

        for _ in 0..n {
            sum += sample_gaussian(&mut rng, mu, sigma);
        }

        let mean = (sum as f64) / (n as f64);
        assert!(
            (mean - mu).abs() < 0.5,
            "Mean {} should be close to {}",
            mean,
            mu
        );
    }

    #[test]
    fn test_base_sampler() {
        let sampler = BaseSampler::new();
        let mut rng = StdRng::seed_from_u64(456);

        // Sample and check distribution
        let mut counts = [0i32; 21]; // -10 to +10
        let n = 10000;

        for _ in 0..n {
            let z = sampler.sample(&mut rng);
            let idx = (z + 10) as usize;
            if idx < counts.len() {
                counts[idx] += 1;
            }
        }

        // z=0 should be most common
        let z0_count = counts[10];
        let z1_count = counts[11];
        let zm1_count = counts[9];

        assert!(
            z0_count > z1_count && z0_count > zm1_count,
            "z=0 should be most common"
        );
    }

    #[test]
    fn test_sampler_z() {
        let sampler = SamplerZ::new();
        let mut rng = StdRng::seed_from_u64(789);

        // Test with sigma near sigma_min (the regime SamplerZ is used in)
        // In the ffSampling tree, sigma_eff ≈ sigma_sign / tree_sigma ≈ 1.3
        let mu = 3.5;
        let sigma = 1.5; // Slightly above sigma_min

        let mut sum: i64 = 0;
        let mut sum_sq: i64 = 0;
        let n = 5000;

        for _ in 0..n {
            let z = sampler.sample(&mut rng, mu, sigma);
            sum += z;
            sum_sq += z * z;
        }

        let mean = (sum as f64) / (n as f64);
        let variance = (sum_sq as f64) / (n as f64) - mean * mean;

        // Mean should be close to mu
        assert!(
            (mean - mu).abs() < 0.5,
            "Mean {} should be close to {}",
            mean,
            mu
        );

        // Variance should be in the right ballpark
        let expected_var = sigma * sigma;
        assert!(
            variance > expected_var * 0.3 && variance < expected_var * 3.0,
            "Variance {} should be roughly {}",
            variance,
            expected_var
        );
    }

    #[test]
    fn test_sampler_z_small_sigma() {
        let sampler = SamplerZ::new();
        let mut rng = StdRng::seed_from_u64(999);

        // Test at sigma_min boundary
        let mu = 0.0;
        let sigma = SIGMA_MIN;

        let mut sum: i64 = 0;
        let n = 5000;

        for _ in 0..n {
            let z = sampler.sample(&mut rng, mu, sigma);
            sum += z;
            // Should produce small integers centered at 0
            assert!(
                z.abs() < 20,
                "z={} too large for sigma_min={}",
                z,
                sigma
            );
        }

        let mean = (sum as f64) / (n as f64);
        assert!(mean.abs() < 0.5, "Mean {} should be close to 0", mean);
    }

    #[test]
    fn test_bernoulli_exp() {
        let mut rng = StdRng::seed_from_u64(111);

        // For x=0, should always return true
        for _ in 0..100 {
            assert!(sample_bernoulli_exp(&mut rng, 0.0));
        }

        // For x=0.5, should return true about exp(-0.5) ≈ 0.607 of the time
        let mut count = 0;
        let n = 10000;
        for _ in 0..n {
            if sample_bernoulli_exp(&mut rng, 0.5) {
                count += 1;
            }
        }
        let rate = (count as f64) / (n as f64);
        let expected = (-0.5f64).exp();
        assert!(
            (rate - expected).abs() < 0.05,
            "Rate {} should be close to {}",
            rate,
            expected
        );
    }
}
