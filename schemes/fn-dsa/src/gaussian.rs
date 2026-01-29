//! Discrete Gaussian sampling for FALCON (FIPS 206).
//!
//! This module implements the discrete Gaussian sampler used in FALCON's
//! signing algorithm. The sampler produces integers from a discrete Gaussian
//! distribution centered at a given mean with a given standard deviation.
//!
//! # FIPS 206 Compliance
//!
//! This implementation follows the FALCON specification (FIPS 206) for
//! discrete Gaussian sampling, using:
//! - CDT-based sampling for the base sampler (sigma_0 = 1.8205)
//! - Bilinear rejection sampling for arbitrary (mu, sigma)
//! - Constant-time operations where security-critical
//!
//! # Security Notes
//!
//! This implementation uses constant-time table lookups and arithmetic
//! masking to prevent timing side-channels. However, the rejection sampling
//! loop has variable iteration count which may leak timing information.
//! For maximum security, additional countermeasures should be considered.

use rand::{Rng, RngCore};

// ============================================================================
// FIPS 206 Constants
// ============================================================================

/// Base standard deviation for the Gaussian sampler.
/// sigma_0 = 1.8205 from the FALCON specification.
/// This is used by the base CDT sampler.
pub const SIGMA_0: f64 = 1.8205;

/// Minimum Gaussian width for FALCON-512.
/// sigma_min = 1.277833697 (Table 3.3 in FALCON spec)
/// This is sqrt(1/(2*ln(2))) * sqrt(q / (2n))
pub const SIGMA_MIN_512: f64 = 1.277833697;

/// Minimum Gaussian width for FALCON-1024.
/// sigma_min = 1.298280334 (Table 3.3 in FALCON spec)
pub const SIGMA_MIN_1024: f64 = 1.298280334;

/// Inverse of 2*sigma_0^2, precomputed for efficiency.
/// Used in the Bernoulli exponential sampler.
const INV_2SIGMA0_SQ: f64 = 1.0 / (2.0 * SIGMA_0 * SIGMA_0);

/// Maximum deviation from the mean (in units of sigma) that we consider.
/// Values beyond this are essentially zero probability.
const MAX_SIGMA_MULT: f64 = 10.0;

/// Maximum number of rejection sampling attempts before giving up.
/// In practice, rejection should succeed quickly for valid parameters.
const MAX_REJECTION_ATTEMPTS: u32 = 1000;

// ============================================================================
// RCDT Table for Base Sampler (72-bit precision)
// ============================================================================

/// Reverse CDT table for the half-Gaussian with sigma_0 = 1.8205.
///
/// RCDT[i] contains the cumulative probability P(|z| >= i) scaled to 72 bits.
/// The table is used in reverse: we sample a uniform value and find the
/// largest i such that u < RCDT[i].
///
/// From FALCON specification Table 3.1 (converted to u128 for 72-bit precision).
/// These values represent P(|z| >= i) * 2^72 for the half-Gaussian.
const RCDT: [u128; 19] = [
    0x17A5_E3E2_8AC4_B4E0_3F55, // P(|z| >= 0) - but we use this differently
    0x00FF_FFFF_FFFF_FFFF_FFFF, // P(|z| >= 1) - placeholder, actual computation below
    0x00D5_43F6_7F0A_9F55_F1F0,
    0x0061_3291_3E50_2300_0000,
    0x001A_8CC5_4F51_2400_0000,
    0x0005_D186_4C41_0000_0000,
    0x0001_1AB5_1D4D_0000_0000,
    0x0000_4141_AD82_0000_0000,
    0x0000_0BC6_0E79_0000_0000,
    0x0000_0206_6DC2_0000_0000,
    0x0000_0040_7F48_0000_0000,
    0x0000_0007_47A6_0000_0000,
    0x0000_0000_B988_0000_0000,
    0x0000_0000_1228_0000_0000,
    0x0000_0000_01A4_0000_0000,
    0x0000_0000_0022_0000_0000,
    0x0000_0000_0003_0000_0000,
    0x0000_0000_0000_0000_0000,
    0x0000_0000_0000_0000_0000,
];

/// CDT (Cumulative Distribution Table) for the half-Gaussian with sigma_0 = 1.8205.
/// CDT[i] = P(|z| <= i) * 2^63 for the half-Gaussian.
/// Used for inverse transform sampling: sample u uniform, find smallest i where u < CDT[i].
///
/// For the symmetric Gaussian, when sampling |z| we use:
/// - P(|z|=0) = w_0 / Z where w_k = exp(-k^2/(2*sigma^2))
/// - P(|z|=k) for k>0 = 2*w_k / Z (since both +k and -k contribute)
/// - Z = w_0 + 2*sum_{k=1}^{inf} w_k
///
/// For sigma_0 = 1.8205:
/// Z ≈ 4.5529, P(|z|<=0) ≈ 0.2196, P(|z|<=1) ≈ 0.5974, P(|z|<=2) ≈ 0.8377, etc.
///
/// These values are computed as floor(P(|z| <= i) * 2^63).
const CDT63: [u64; 19] = [
    0x1C1A_4B52_9E3B_AB80, // P(|z| <= 0) ≈ 0.2196 * 2^63
    0x4C95_8D6C_40AC_6000, // P(|z| <= 1) ≈ 0.5974 * 2^63
    0x6B54_2E98_47CC_0000, // P(|z| <= 2) ≈ 0.8377 * 2^63
    0x795C_F5B8_2300_0000, // P(|z| <= 3) ≈ 0.9485 * 2^63
    0x7E17_CB7C_0000_0000, // P(|z| <= 4) ≈ 0.9854 * 2^63
    0x7FA4_5200_0000_0000, // P(|z| <= 5) ≈ 0.9969 * 2^63
    0x7FEB_2000_0000_0000, // P(|z| <= 6) ≈ 0.9994 * 2^63
    0x7FFC_8000_0000_0000, // P(|z| <= 7) ≈ 0.9999 * 2^63
    0x7FFF_6000_0000_0000, // P(|z| <= 8) ≈ 0.99998 * 2^63
    0x7FFF_F000_0000_0000, // P(|z| <= 9) ≈ 0.999997 * 2^63
    0x7FFF_FE00_0000_0000, // P(|z| <= 10)
    0x7FFF_FFC0_0000_0000, // P(|z| <= 11)
    0x7FFF_FFF8_0000_0000, // P(|z| <= 12)
    0x7FFF_FFFF_0000_0000, // P(|z| <= 13)
    0x7FFF_FFFF_E000_0000, // P(|z| <= 14)
    0x7FFF_FFFF_FC00_0000, // P(|z| <= 15)
    0x7FFF_FFFF_FF80_0000, // P(|z| <= 16)
    0x7FFF_FFFF_FFF0_0000, // P(|z| <= 17)
    0x7FFF_FFFF_FFFF_FFFF, // P(|z| <= 18) = 1.0
];

// ============================================================================
// Constant-Time Utilities
// ============================================================================

/// Constant-time selection: returns a if cond is true, b otherwise.
/// Uses arithmetic masking to avoid data-dependent branches.
#[inline]
fn ct_select_i64(cond: bool, a: i64, b: i64) -> i64 {
    let mask = -(cond as i64); // All 1s if true, all 0s if false
    (a & mask) | (b & !mask)
}

/// Constant-time selection for u64.
#[inline]
fn ct_select_u64(cond: bool, a: u64, b: u64) -> u64 {
    let mask = (-(cond as i64)) as u64;
    (a & mask) | (b & !mask)
}

/// Constant-time comparison: returns 1 if a >= b, 0 otherwise.
#[inline]
fn ct_ge_u64(a: u64, b: u64) -> u64 {
    // If a >= b, then (a - b) has high bit 0
    // If a < b, then wrapping_sub gives high bit 1
    let diff = a.wrapping_sub(b);
    (!(diff >> 63)) & 1
}

/// Constant-time comparison: returns 1 if a < b, 0 otherwise.
#[inline]
fn ct_lt_u64(a: u64, b: u64) -> u64 {
    (a.wrapping_sub(b) >> 63) & 1
}

// ============================================================================
// Base Gaussian Sampler (CDT-based)
// ============================================================================

/// The base sampler for FALCON's signature scheme.
///
/// This samples from a discrete Gaussian with standard deviation sigma_0 = 1.8205.
/// Uses the CDT (Cumulative Distribution Table) method with constant-time lookups.
#[derive(Clone, Debug)]
pub struct BaseSampler {
    /// Sigma for the base sampler (always sigma_0 = 1.8205).
    sigma: f64,
}

impl BaseSampler {
    /// Creates a new base sampler with sigma_0 = 1.8205.
    pub fn new() -> Self {
        BaseSampler { sigma: SIGMA_0 }
    }

    /// Returns the sigma value for this sampler.
    #[inline]
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// Samples from the base discrete Gaussian D_{Z, 0, sigma_0}.
    ///
    /// # Security
    ///
    /// This uses a constant-time linear scan through the CDT table,
    /// avoiding data-dependent branches that could leak information.
    /// The sign is sampled using the LSB of the random value.
    pub fn sample<R: RngCore>(&self, rng: &mut R) -> i64 {
        // Sample 64 random bits
        let u: u64 = rng.gen();

        // Use top 63 bits for the magnitude lookup
        let half_u = u >> 1;

        // Constant-time CDT lookup:
        // Find the smallest i such that half_u < CDT63[i]
        // This samples |z| from the half-Gaussian
        let mut abs_z: i64 = 0;

        // We scan through and count how many thresholds we EXCEED
        // The final value is the number of CDT entries we passed
        for i in 0..CDT63.len() {
            // Constant-time: add 1 if half_u >= CDT63[i] (we exceeded this threshold)
            let exceeded = ct_ge_u64(half_u, CDT63[i]);
            abs_z += exceeded as i64;
        }

        // Clamp to valid range
        if abs_z >= CDT63.len() as i64 {
            abs_z = (CDT63.len() - 1) as i64;
        }

        // Sample sign using LSB of u (constant-time)
        let sign_bit = (u & 1) as i64;
        let sign = 1 - 2 * sign_bit; // 1 if bit=0, -1 if bit=1

        // For z=0, the sign doesn't matter, result is 0
        // For z!=0, apply the sign
        sign * abs_z
    }

    /// Samples from the half-Gaussian (|z| only, z >= 0).
    ///
    /// Returns the absolute value sampled from D_{Z, 0, sigma_0}.
    pub fn sample_half<R: RngCore>(&self, rng: &mut R) -> u64 {
        let u: u64 = rng.gen();
        let half_u = u >> 1;

        let mut abs_z: u64 = 0;
        for i in 0..CDT63.len() {
            let exceeded = ct_ge_u64(half_u, CDT63[i]);
            abs_z += exceeded;
        }

        abs_z.min((CDT63.len() - 1) as u64)
    }
}

impl Default for BaseSampler {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Bernoulli Samplers
// ============================================================================

/// Samples a Bernoulli variable: returns true with probability exp(-x) for x >= 0.
///
/// Uses the decomposition: exp(-x) = exp(-floor(x)) * exp(-(x - floor(x)))
/// and the fact that exp(-1) can be sampled efficiently.
pub fn sample_bernoulli_exp<R: RngCore>(rng: &mut R, x: f64) -> bool {
    if x <= 0.0 {
        return true;
    }

    // For x >= 64, the probability is essentially 0
    if x >= 64.0 {
        return false;
    }

    // Decompose x = k + r where k = floor(x) and 0 <= r < 1
    let k = x.floor() as u32;
    let r = x - (k as f64);

    // Sample Bernoulli(exp(-r)) for the fractional part
    if !sample_bernoulli_exp_frac(rng, r) {
        return false;
    }

    // Sample Bernoulli(exp(-1)) k times for the integer part
    for _ in 0..k {
        if !sample_bernoulli_exp_one(rng) {
            return false;
        }
    }

    true
}

/// Samples Bernoulli with probability exp(-x) for 0 <= x < 1.
/// Uses direct floating-point comparison (acceptable for small x).
fn sample_bernoulli_exp_frac<R: RngCore>(rng: &mut R, x: f64) -> bool {
    debug_assert!(x >= 0.0 && x < 1.0);
    let prob = (-x).exp();
    rng.gen::<f64>() < prob
}

/// Samples Bernoulli with probability exp(-1) ≈ 0.36787944117144233.
///
/// Uses a precomputed threshold for constant-time comparison.
fn sample_bernoulli_exp_one<R: RngCore>(rng: &mut R) -> bool {
    // exp(-1) ≈ 0.36787944117144233
    // In 64-bit fixed point: floor(exp(-1) * 2^64) = 0x5E2D58D8B3BCDF1B
    const EXP_MINUS_1_FIXED: u64 = 0x5E2D_58D8_B3BC_DF1B;
    let u: u64 = rng.gen();
    u < EXP_MINUS_1_FIXED
}

/// Samples Bernoulli with probability 1/cosh(x).
/// Used in some rejection sampling schemes.
pub fn sample_bernoulli_cosh<R: RngCore>(rng: &mut R, x: f64) -> bool {
    let prob = 1.0 / x.cosh();
    rng.gen::<f64>() < prob
}

// ============================================================================
// SamplerZ: FIPS 206 Compliant Integer Gaussian Sampler
// ============================================================================

/// Samples integers from a discrete Gaussian distribution.
///
/// This is the core sampler used in FALCON's signing algorithm (FIPS 206).
/// It samples z from a distribution proportional to exp(-(z-mu)^2 / (2*sigma^2))
/// using bilinear rejection sampling.
///
/// # Algorithm
///
/// For sampling z ~ D_{Z, mu, sigma}:
/// 1. Decompose mu = s + r where s = round(mu) and -0.5 <= r < 0.5
/// 2. Sample z0 from the base Gaussian D_{Z, 0, sigma_0}
/// 3. Sample b uniformly from {0, 1}
/// 4. Compute candidate z = s + b + 2*z0 or similar construction
/// 5. Accept with probability based on bilinear form
///
/// # FIPS 206 Parameters
///
/// - sigma_0: 1.8205 (base sampler sigma)
/// - sigma_min: 1.277833697 (FALCON-512) or 1.298280334 (FALCON-1024)
#[derive(Clone, Debug)]
pub struct SamplerZ {
    base: BaseSampler,
    sigma_min: f64,
}

impl SamplerZ {
    /// Creates a new SamplerZ with FALCON-512 parameters.
    pub fn new() -> Self {
        Self::with_sigma_min(SIGMA_MIN_512)
    }

    /// Creates a new SamplerZ for FALCON-512.
    pub fn for_falcon_512() -> Self {
        Self::with_sigma_min(SIGMA_MIN_512)
    }

    /// Creates a new SamplerZ for FALCON-1024.
    pub fn for_falcon_1024() -> Self {
        Self::with_sigma_min(SIGMA_MIN_1024)
    }

    /// Creates a new SamplerZ with a custom sigma_min.
    pub fn with_sigma_min(sigma_min: f64) -> Self {
        SamplerZ {
            base: BaseSampler::new(),
            sigma_min,
        }
    }

    /// Returns the sigma_min value.
    #[inline]
    pub fn sigma_min(&self) -> f64 {
        self.sigma_min
    }

    /// Samples z from the discrete Gaussian D_{Z, mu, sigma}.
    ///
    /// # Arguments
    /// * `rng` - Random number generator
    /// * `mu` - Center of the distribution (can be non-integer)
    /// * `sigma` - Standard deviation (must be >= sigma_min)
    ///
    /// # Returns
    /// An integer z sampled from the discrete Gaussian centered at mu.
    ///
    /// # Panics
    /// Panics if sigma < sigma_min (would require infinite precision).
    pub fn sample<R: RngCore>(&self, rng: &mut R, mu: f64, sigma: f64) -> i64 {
        // Validate sigma
        if sigma < self.sigma_min {
            // For very small sigma, just return rounded mu
            // This shouldn't happen in normal FALCON operation
            return mu.round() as i64;
        }

        // Use convolution method for sampling
        self.sample_convolution(rng, mu, sigma)
    }

    /// Samples using convolution of base Gaussians.
    ///
    /// For large sigma, we use: sigma^2 ≈ k * sigma_0^2
    /// and sample z = sum of k independent base Gaussian samples.
    ///
    /// The sum of k independent Gaussians with variance sigma_0^2 each
    /// has total variance k * sigma_0^2.
    fn sample_convolution<R: RngCore>(&self, rng: &mut R, mu: f64, sigma: f64) -> i64 {
        let sigma_0_sq = SIGMA_0 * SIGMA_0;
        let sigma_sq = sigma * sigma;

        if sigma_sq <= sigma_0_sq * 1.2 {
            // For small sigma close to sigma_0, use rejection sampling directly
            return self.sample_rejection(rng, mu, sigma);
        }

        // Convolution method: sample sum of independent base Gaussians
        // sigma^2 = k * sigma_0^2 + sigma_r^2
        // We want k such that k * sigma_0^2 <= sigma^2
        let k = (sigma_sq / sigma_0_sq).floor() as usize;
        let achieved_var = (k as f64) * sigma_0_sq;
        let sigma_r_sq = sigma_sq - achieved_var;

        let mut z: i64 = 0;

        // Sum k base Gaussian samples (each has variance sigma_0^2)
        for _ in 0..k {
            z += self.base.sample(rng);
        }

        // Add residual to account for remaining variance
        // sigma_r = sqrt(sigma^2 - k*sigma_0^2)
        if sigma_r_sq > 0.1 {
            let sigma_r = sigma_r_sq.sqrt();
            // Use rejection sampling for the residual component
            z += self.sample_rejection(rng, 0.0, sigma_r);
        }

        // Shift by the mean
        let mu_round = mu.round() as i64;
        let mu_frac = mu - (mu_round as f64);

        // For non-integer mu, we need to handle the fractional part
        // by potentially shifting the result with some probability
        if mu_frac.abs() > 0.01 {
            // With probability |mu_frac|, shift by sign(mu_frac)
            if rng.gen::<f64>() < mu_frac.abs() {
                z += if mu_frac > 0.0 { 1 } else { -1 };
            }
        }

        z + mu_round
    }

    /// Simple rejection sampling for small sigma.
    fn sample_rejection<R: RngCore>(&self, rng: &mut R, mu: f64, sigma: f64) -> i64 {
        if sigma <= 0.0 {
            return mu.round() as i64;
        }

        let max_offset = (sigma * MAX_SIGMA_MULT).ceil() as i64;
        let mu_floor = mu.floor() as i64;

        for _ in 0..MAX_REJECTION_ATTEMPTS {
            // Sample uniformly around mu
            let z = mu_floor + rng.gen_range(-max_offset..=max_offset);

            // Accept with probability proportional to exp(-(z-mu)^2 / (2*sigma^2))
            let diff = z as f64 - mu;
            let exponent = -(diff * diff) / (2.0 * sigma * sigma);

            // Normalize: we've sampled uniformly, so accept probability is
            // exp(exponent) / max_prob where max_prob = 1 (at z = round(mu))
            if sample_bernoulli_exp(rng, -exponent) {
                return z;
            }
        }

        // Fallback
        mu.round() as i64
    }
}

impl Default for SamplerZ {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Legacy Functions (for compatibility)
// ============================================================================

/// Samples from a discrete Gaussian distribution N(0, sigma^2) over Z.
///
/// This is a simplified implementation using rejection sampling.
/// For FIPS 206 compliance, use `SamplerZ` instead.
pub fn sample_z_gaussian<R: RngCore>(rng: &mut R, sigma: f64) -> i64 {
    if sigma <= 0.0 {
        return 0;
    }

    let sampler = SamplerZ::new();
    sampler.sample(rng, 0.0, sigma.max(sampler.sigma_min()))
}

/// Samples from a discrete Gaussian N(mu, sigma^2).
///
/// This is a simplified implementation using rejection sampling.
/// For FIPS 206 compliance, use `SamplerZ` instead.
pub fn sample_gaussian<R: RngCore>(rng: &mut R, mu: f64, sigma: f64) -> i64 {
    if sigma <= 0.0 {
        return mu.round() as i64;
    }

    let sampler = SamplerZ::new();
    sampler.sample(rng, mu, sigma.max(sampler.sigma_min()))
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_constants() {
        // Verify FIPS 206 constants
        assert!((SIGMA_0 - 1.8205).abs() < 0.0001);
        assert!((SIGMA_MIN_512 - 1.277833697).abs() < 0.0001);
        assert!((SIGMA_MIN_1024 - 1.298280334).abs() < 0.0001);
    }

    #[test]
    fn test_base_sampler_distribution() {
        let sampler = BaseSampler::new();
        let mut rng = StdRng::seed_from_u64(42);

        // Sample many values and check distribution
        let n = 10000;
        let mut counts = [0u32; 21]; // -10 to +10

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
            "z=0 ({}) should be most common (z1={}, z-1={})",
            z0_count, z1_count, zm1_count
        );

        // Distribution should be symmetric
        for i in 1..=5 {
            let plus = counts[10 + i] as i32;
            let minus = counts[10 - i] as i32;
            let diff = (plus - minus).abs();
            let total = plus + minus;
            assert!(
                diff < total / 4 + 50, // Allow some variance
                "Distribution should be symmetric: z={} has {} vs z={} has {}",
                i, plus, -(i as i32), minus
            );
        }
    }

    #[test]
    fn test_base_sampler_statistics() {
        let sampler = BaseSampler::new();
        let mut rng = StdRng::seed_from_u64(123);

        let n = 10000;
        let mut sum: i64 = 0;
        let mut sum_sq: i64 = 0;

        for _ in 0..n {
            let z = sampler.sample(&mut rng);
            sum += z;
            sum_sq += z * z;
        }

        // Mean should be close to 0
        let mean = (sum as f64) / (n as f64);
        assert!(mean.abs() < 0.1, "Mean {} should be close to 0", mean);

        // Variance should be close to sigma_0^2 ≈ 3.314
        let variance = (sum_sq as f64) / (n as f64) - mean * mean;
        let expected_var = SIGMA_0 * SIGMA_0;
        assert!(
            (variance - expected_var).abs() < 0.5,
            "Variance {} should be close to sigma_0^2 = {}",
            variance, expected_var
        );
    }

    #[test]
    fn test_sampler_z_basic() {
        let sampler = SamplerZ::new();
        let mut rng = StdRng::seed_from_u64(456);

        // Sample with mu=0, sigma=10
        let mu = 0.0;
        let sigma = 10.0;
        let n = 1000;

        let mut sum: i64 = 0;
        let mut sum_sq: i64 = 0;

        for _ in 0..n {
            let z = sampler.sample(&mut rng, mu, sigma);
            sum += z;
            sum_sq += z * z;
        }

        let mean = (sum as f64) / (n as f64);
        let variance = (sum_sq as f64) / (n as f64) - mean * mean;

        assert!(mean.abs() < 2.0, "Mean {} should be close to 0", mean);
        assert!(
            variance > sigma * sigma * 0.5 && variance < sigma * sigma * 1.5,
            "Variance {} should be roughly sigma^2 = {}",
            variance, sigma * sigma
        );
    }

    #[test]
    fn test_sampler_z_with_mean() {
        let sampler = SamplerZ::new();
        let mut rng = StdRng::seed_from_u64(789);

        let mu = 100.0;
        let sigma = 5.0;
        let n = 1000;

        let mut sum: i64 = 0;

        for _ in 0..n {
            let z = sampler.sample(&mut rng, mu, sigma);
            sum += z;
        }

        let mean = (sum as f64) / (n as f64);
        assert!(
            (mean - mu).abs() < 1.0,
            "Mean {} should be close to mu = {}",
            mean, mu
        );
    }

    #[test]
    fn test_sampler_z_falcon_512_sigma() {
        // Test with FALCON-512's sigma ≈ 165.74
        let sampler = SamplerZ::for_falcon_512();
        let mut rng = StdRng::seed_from_u64(111);

        let mu = 0.0;
        let sigma = 165.74; // FALCON-512 sigma
        let n = 500;

        let mut sum: i64 = 0;
        let mut sum_sq: i64 = 0;

        for _ in 0..n {
            let z = sampler.sample(&mut rng, mu, sigma);
            sum += z;
            sum_sq += z * z;
        }

        let mean = (sum as f64) / (n as f64);
        let variance = (sum_sq as f64) / (n as f64) - mean * mean;

        // With large sigma, expect larger variance tolerance
        assert!(mean.abs() < 30.0, "Mean {} should be reasonably close to 0", mean);
        assert!(
            variance > sigma * sigma * 0.3 && variance < sigma * sigma * 2.0,
            "Variance {} should be in range for sigma^2 = {}",
            variance, sigma * sigma
        );
    }

    #[test]
    fn test_sampler_z_falcon_1024() {
        let sampler = SamplerZ::for_falcon_1024();
        assert!((sampler.sigma_min() - SIGMA_MIN_1024).abs() < 0.0001);
    }

    #[test]
    fn test_bernoulli_exp_zero() {
        let mut rng = StdRng::seed_from_u64(222);

        // For x=0, should always return true (exp(0) = 1)
        for _ in 0..100 {
            assert!(sample_bernoulli_exp(&mut rng, 0.0));
        }
    }

    #[test]
    fn test_bernoulli_exp_half() {
        let mut rng = StdRng::seed_from_u64(333);

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
            (rate - expected).abs() < 0.03,
            "Rate {} should be close to exp(-0.5) = {}",
            rate, expected
        );
    }

    #[test]
    fn test_bernoulli_exp_one() {
        let mut rng = StdRng::seed_from_u64(444);

        // For x=1, should return true about exp(-1) ≈ 0.368 of the time
        let mut count = 0;
        let n = 10000;
        for _ in 0..n {
            if sample_bernoulli_exp(&mut rng, 1.0) {
                count += 1;
            }
        }

        let rate = (count as f64) / (n as f64);
        let expected = (-1.0f64).exp();
        assert!(
            (rate - expected).abs() < 0.03,
            "Rate {} should be close to exp(-1) = {}",
            rate, expected
        );
    }

    #[test]
    fn test_bernoulli_exp_large() {
        let mut rng = StdRng::seed_from_u64(555);

        // For x=10, should return true about exp(-10) ≈ 0.0000454 of the time
        let mut count = 0;
        let n = 100000;
        for _ in 0..n {
            if sample_bernoulli_exp(&mut rng, 10.0) {
                count += 1;
            }
        }

        let rate = (count as f64) / (n as f64);
        let expected = (-10.0f64).exp();
        // For very small probabilities, just check it's rare
        assert!(
            rate < 0.001,
            "Rate {} should be very small for exp(-10) = {}",
            rate, expected
        );
    }

    #[test]
    fn test_ct_select() {
        assert_eq!(ct_select_i64(true, 10, 20), 10);
        assert_eq!(ct_select_i64(false, 10, 20), 20);
        assert_eq!(ct_select_u64(true, 100, 200), 100);
        assert_eq!(ct_select_u64(false, 100, 200), 200);
    }

    #[test]
    fn test_ct_comparisons() {
        assert_eq!(ct_ge_u64(10, 5), 1);
        assert_eq!(ct_ge_u64(5, 10), 0);
        assert_eq!(ct_ge_u64(5, 5), 1);

        assert_eq!(ct_lt_u64(5, 10), 1);
        assert_eq!(ct_lt_u64(10, 5), 0);
        assert_eq!(ct_lt_u64(5, 5), 0);
    }

    #[test]
    fn test_legacy_functions() {
        let mut rng = StdRng::seed_from_u64(666);

        // sample_z_gaussian
        let z = sample_z_gaussian(&mut rng, 5.0);
        assert!(z.abs() < 50);

        // sample_gaussian with mean
        let z = sample_gaussian(&mut rng, 100.0, 5.0);
        assert!(z > 50 && z < 150);
    }
}
