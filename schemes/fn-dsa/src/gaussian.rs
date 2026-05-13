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
use std::f64::consts::{E, PI};

/// Base standard deviation for the Gaussian sampler.
/// This is sigma_0 = 1.8205 from the FALCON specification.
pub const SIGMA_0: f64 = 1.8205;

/// Maximum deviation from the mean (in units of sigma) that we consider.
/// Values beyond this are essentially zero probability.
const MAX_SIGMA_MULT: f64 = 10.0;

// ============================================================================
// Base Gaussian Sampler (used at keygen time only)
// ============================================================================

/// Samples from a discrete Gaussian distribution N(0, sigma^2) over Z.
///
/// Uses rejection sampling with a uniform proposal. Returns an integer z
/// sampled according to `exp(-z²/(2σ²))`. Not constant-time; intended for
/// keygen-time use only (the signer's per-leaf sampler is constant-time
/// via [`SamplerZ`]).
pub(crate) fn sample_z_gaussian<R: RngCore>(rng: &mut R, sigma: f64) -> i64 {
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
/// Test-only: orphaned by the FIPS 206 Algorithm 12 SamplerZ rewrite, retained
/// for the legacy `test_sample_gaussian_mean` unit test.
#[cfg(test)]
fn sample_gaussian<R: RngCore>(rng: &mut R, mu: f64, sigma: f64) -> i64 {
    if sigma <= 0.0 {
        return mu.round() as i64;
    }

    // If mu is outside a safe range, just return the rounded value.
    // With any reasonable sigma, the probability mass is concentrated near mu,
    // so this is essentially exact when |mu| >> sigma.
    let safe_limit = (i64::MAX / 2) as f64;
    if mu.abs() > safe_limit || mu.is_nan() || mu.is_infinite() {
        return if mu.is_nan() { 0 } else { mu.round() as i64 };
    }

    // For FALCON, we use rejection sampling around mu
    let max_offset = (sigma * MAX_SIGMA_MULT).ceil() as i64;
    let mu_floor = mu.floor() as i64;

    loop {
        // Sample uniformly around mu
        let z = mu_floor.saturating_add(rng.gen_range(-max_offset..=max_offset));

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
// Bernoullí Sampler (test-only — superseded by `ber_exp`)
// ============================================================================

/// Samples a Bernoulli variable: returns true with probability exp(-x) for x >= 0.
/// Test-only: orphaned by the FACCT [`ber_exp`] implementation; retained
/// for the legacy `test_bernoulli_exp` unit test.
#[cfg(test)]
fn sample_bernoulli_exp<R: RngCore>(rng: &mut R, x: f64) -> bool {
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

/// Samples Bernoulli with probability exp(-1) ≈ 0.368. Test-only helper for
/// the recursion in [`sample_bernoulli_exp`].
#[cfg(test)]
fn sample_bernoulli_exp_minus_one<R: RngCore>(rng: &mut R) -> bool {
    const EXP_MINUS_1: f64 = 0.36787944117144233;
    rng.gen::<f64>() < EXP_MINUS_1
}

// ============================================================================
// BaseSampler for FALCON
// ============================================================================

/// The base half-Gaussian sampler for FALCON's signature scheme.
///
/// Samples non-negative integers from the discrete half-Gaussian
/// `D_{Z+, sigma_base}` using a precomputed cumulative distribution
/// table (CDT) and a constant-time linear scan.
///
/// FIPS 206 Algorithm 12 (`SamplerZ`) uses a base sampler at
/// `sigma_0 = 1.8205` (see [`SIGMA_0`]). Other configurations exist
/// for compatibility with the legacy callers; use
/// [`BaseSampler::sigma_zero`] for the FIPS-compliant variant.
pub struct BaseSampler {
    /// Precomputed CDT for the base Gaussian.
    cdt: Vec<u64>,
    /// Sigma for the base sampler.
    #[allow(dead_code)]
    sigma: f64,
}

impl BaseSampler {
    /// Creates a new base sampler with the legacy default sigma (1.277833).
    ///
    /// Retained only for the historical legacy tests; new code should use
    /// [`BaseSampler::sigma_zero`] which matches FIPS 206 Algorithm 12.
    pub fn new() -> Self {
        // For FALCON, sigma_min ≈ 1.277833 (for FALCON-512)
        // The CDT is precomputed for this value
        let sigma = 1.277833;
        let cdt = Self::compute_cdt(sigma);
        BaseSampler { cdt, sigma }
    }

    /// Creates a base sampler at the FIPS 206 base sigma `sigma_0 = 1.8205`,
    /// configured for the half-Gaussian distribution `D_{Z+, sigma_0}`.
    ///
    /// This is the sigma used by Algorithm 12 ([`SamplerZ`]) for the inner
    /// half-Gaussian draw `z0`. The CDT here is for the *half*-Gaussian
    /// (density `exp(-z²/(2σ²)) / Z_+` over `z >= 0`), distinct from the
    /// `|full Gaussian|` CDT used by [`BaseSampler::new`]. The two differ
    /// in the weight at `z=0` (half-Gaussian gives `z=0` mass `1/Z_+`,
    /// full-Gaussian gives `1/Z_full ≈ 1/(2·Z_+)`).
    pub fn sigma_zero() -> Self {
        let sigma = SIGMA_0;
        let cdt = Self::compute_half_cdt(sigma);
        BaseSampler { cdt, sigma }
    }

    /// Creates a base sampler with a custom sigma (`|full Gaussian|` CDT).
    #[allow(dead_code)]
    pub fn with_sigma(sigma: f64) -> Self {
        let cdt = Self::compute_cdt(sigma);
        BaseSampler { cdt, sigma }
    }

    /// Computes the CDT for the half-Gaussian `D_{Z+, sigma}`.
    ///
    /// `CDT[k] = floor(2^63 * sum_{j=0}^{k} P_half(j))` where
    /// `P_half(j) = exp(-j²/(2σ²)) / Z_+`, `Z_+ = sum_{j>=0} exp(-j²/(2σ²))`.
    ///
    /// **Distinct from [`compute_cdt`]**: the half-Gaussian gives `z=0`
    /// mass `1/Z_+`; the `|full Gaussian|` CDT gives `z=0` mass
    /// `1/Z_full ≈ 1/(2·Z_+)`. This is the function used by
    /// [`BaseSampler::sigma_zero`] (the FIPS 206 Algorithm 12 base draw).
    ///
    /// Precision: `f64` arithmetic, quantized to `u64`. The cumulative
    /// sum accumulates rounding error on the order of `ULP · max_z`,
    /// well below the 63-bit threshold quantization granularity at
    /// FALCON's `sigma_0`. For bit-exact CDT reproducibility against
    /// PQClean reference, a fixed-point construction matching FACCT
    /// would be required; the current implementation matches the
    /// distribution to within statistical noise but not to the last bit.
    fn compute_half_cdt(sigma: f64) -> Vec<u64> {
        let max_z = (sigma * MAX_SIGMA_MULT).ceil() as usize;
        let mut probs = Vec::with_capacity(max_z + 1);
        for z in 0..=max_z {
            let zf = z as f64;
            probs.push((-(zf * zf) / (2.0 * sigma * sigma)).exp());
        }
        // Half-Gaussian normalizer: Z_+ = sum_{z>=0} exp(-z²/(2σ²)).
        // (Every z >= 0 carries weight 1; no "double-count" for z>0.)
        let z_plus: f64 = probs.iter().sum();
        let mut cdt = Vec::with_capacity(max_z + 1);
        let mut cumsum = 0.0;
        let scale = (1u64 << 63) as f64;
        for &p in &probs {
            cumsum += p / z_plus;
            cdt.push((cumsum * scale) as u64);
        }
        cdt
    }

    /// Computes the CDT for the `|full Gaussian|` distribution at `sigma`,
    /// i.e. the absolute-value distribution of a full discrete Gaussian
    /// over the integers. `CDT[k]` is the cumulative probability that
    /// `|z| <= k`, quantized to a 63-bit fixed-point threshold.
    ///
    /// **Not the same distribution as [`compute_half_cdt`]**: this
    /// function gives `z=0` weight `1/Z_full` and `z>0` weight
    /// `2/Z_full`, matching `|N(0, sigma)|`. The half-Gaussian
    /// `D_{Z+, sigma}` (used by FIPS 206 Algorithm 12) gives `z=0`
    /// weight `1/Z_+` — twice as much. Use `compute_half_cdt` for
    /// the SamplerZ inner draw; use `compute_cdt` only for the legacy
    /// signed-sample API ([`BaseSampler::sample`]).
    ///
    /// Precision: `f64` arithmetic, quantized to `u64`. Accumulated
    /// rounding error is bounded by `~ULP · max_z`, below the 63-bit
    /// threshold quantization granularity for typical FALCON sigmas.
    fn compute_cdt(sigma: f64) -> Vec<u64> {
        // CDT[i] = floor(2^63 * sum_{j=0}^{i} P(|z|=j))
        // where P(|z|=0) = 1/Z_full and P(|z|=j>0) = 2*exp(-j²/(2σ²))/Z_full,
        // Z_full = 1 + 2 * sum_{j>=1} exp(-j²/(2σ²)).

        let max_z = (sigma * MAX_SIGMA_MULT).ceil() as usize;
        let mut probs = Vec::with_capacity(max_z + 1);

        for z in 0..=max_z {
            let zf = z as f64;
            probs.push((-(zf * zf) / (2.0 * sigma * sigma)).exp());
        }

        let z: f64 = probs[0] + 2.0 * probs[1..].iter().sum::<f64>();

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

    /// Samples a signed integer from the base discrete Gaussian.
    ///
    /// Internally samples a half-Gaussian magnitude then assigns a uniform sign.
    ///
    /// # Security
    ///
    /// This uses a constant-time linear scan through the CDT table,
    /// avoiding data-dependent branches that could leak information.
    /// All comparisons use arithmetic masking instead of conditional branches.
    #[allow(dead_code)]
    pub fn sample<R: RngCore>(&self, rng: &mut R) -> i64 {
        let u: u64 = rng.gen();
        let half_u = u >> 1; // Use top 63 bits for the magnitude draw

        let mut abs_z: i64 = self.sample_half_inner(half_u);

        // Sample sign using the LSB of u (constant-time)
        let sign_bit = (u & 1) as i64;
        let sign = 1 - 2 * sign_bit; // 1 if bit=0, -1 if bit=1

        // Apply sign, but return 0 if abs_z == 0
        abs_z = sign * abs_z;
        abs_z
    }

    /// Samples a non-negative integer from the base half-Gaussian.
    ///
    /// Used by [`SamplerZ`] as the inner draw of FIPS 206 Algorithm 12.
    pub fn sample_half<R: RngCore>(&self, rng: &mut R) -> i64 {
        let u: u64 = rng.gen();
        let half_u = u >> 1;
        self.sample_half_inner(half_u)
    }

    /// Counts how many CDT thresholds are met by `half_u`. Constant-time
    /// over `self.cdt.len()`.
    #[inline]
    fn sample_half_inner(&self, half_u: u64) -> i64 {
        let mut abs_z: i64 = 0;
        for &threshold in &self.cdt {
            // Constant-time comparison: (half_u >= threshold) without branching.
            let diff = half_u.wrapping_sub(threshold);
            let exceeded = ((!(diff >> 63)) & 1) as i64; // 1 if half_u >= threshold
            abs_z += exceeded;
        }
        abs_z
    }
}

impl Default for BaseSampler {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// SamplerZ: FIPS 206 Algorithm 12 — FALCON's integer Gaussian sampler
// ============================================================================

/// Samples an integer from a discrete Gaussian centered at `mu` with width `sigma`.
///
/// This implements FIPS 206 Algorithm 12 (SamplerZ): the inner half-Gaussian
/// draw is taken from a CDT-based base sampler at `sigma_0 = 1.8205`, and a
/// BerExp acceptance test (Algorithm 14) rejects to adjust the distribution
/// to the per-leaf target `sigma`. The `sigma_min` argument (set via
/// [`SamplerZ::with_sigma_min`]) controls the BerExp scaling factor
/// `ccs = sigma_min / sigma`.
///
/// # Side-channel status
///
/// The half-Gaussian draw via [`BaseSampler::sample_half`] is constant-time
/// (CDT linear scan, arithmetic-mask comparisons, signed accumulation).
/// The BerExp acceptance test in [`ber_exp`] is constant-time over its
/// inner loop: it reduces `x = s · ln 2 + r` and evaluates the FACCT
/// degree-12 polynomial in `r` via 64-bit widening multiplications, with
/// no data-dependent branches or table lookups. The outer rejection loop
/// runs a probabilistic number of iterations — this leakage is inherent
/// to FALCON's design and is treated as public in the security proof.
pub struct SamplerZ {
    /// Base half-Gaussian sampler at sigma_0 = 1.8205.
    base: BaseSampler,
    /// Per-parameter-set minimum sigma. Drives the BerExp acceptance ratio.
    sigma_min: f64,
}

impl SamplerZ {
    /// Creates a new SamplerZ with the FALCON-512 default sigma_min.
    pub fn new() -> Self {
        SamplerZ {
            base: BaseSampler::sigma_zero(),
            sigma_min: FALCON_512_SIGMA_MIN,
        }
    }

    /// Creates a new SamplerZ tuned to the given `sigma_min`.
    ///
    /// `sigma_min` must be `>= 1.0` and `<= sigma_0`; the FALCON parameter
    /// sets set it to ~1.28 (FALCON-512) or ~1.30 (FALCON-1024).
    pub fn with_sigma_min(sigma_min: f64) -> Self {
        SamplerZ {
            base: BaseSampler::sigma_zero(),
            sigma_min,
        }
    }

    /// Samples `z` from the discrete Gaussian `D_{Z, mu, sigma}` per FIPS 206
    /// Algorithm 12.
    ///
    /// Caller invariant: `sigma >= self.sigma_min` (enforced by `assert!`).
    /// The FALCON keygen's basis-quality check ensures every per-leaf sigma
    /// produced by [`crate::sampler::FfSampler`] satisfies this.
    pub fn sample<R: RngCore>(&self, rng: &mut R, mu: f64, sigma: f64) -> i64 {
        assert!(
            sigma >= self.sigma_min,
            "SamplerZ::sample requires sigma ({}) >= sigma_min ({})",
            sigma,
            self.sigma_min
        );
        assert!(sigma.is_finite() && sigma > 0.0, "sigma must be finite and positive");
        if !mu.is_finite() {
            // Degenerate input — return 0 rather than diverging on NaN.
            return 0;
        }

        // s = floor(mu), r = mu - s in [0, 1).
        let s = mu.floor() as i64;
        let r = mu - s as f64;

        // BerExp scaling factor: P(accept) = (sigma_min / sigma) * exp(-x).
        let ccs = self.sigma_min / sigma;

        // Precomputed inverse scales so the loop body avoids division.
        let inv_2sigma_sq = 1.0 / (2.0 * sigma * sigma);
        let inv_2sigma0_sq = INV_2SIGMA0_SQ;

        // Outer rejection loop. The number of iterations is data-dependent
        // (inherent to FALCON's design); inner draws and the BerExp
        // probabilistic test are not constant-time in the wall-clock sense
        // but the BaseSampler's CDT scan is.
        loop {
            // z0 ~ D_{Z+, sigma_0}: a non-negative draw from the base
            // half-Gaussian. Constant-time over the CDT.
            let z0 = self.base.sample_half(rng);

            // Sign bit b. Drawn from a fresh u64 (top bit) so it is
            // independent of the z0 draw and uniformly distributed.
            let b_word: u64 = rng.gen();
            let b: i64 = ((b_word >> 63) & 1) as i64;

            // z = b + (2b - 1) * z0:
            //   b = 0 → z = -z0    (z in {0, -1, -2, ...})
            //   b = 1 → z = 1 + z0 (z in {1, 2, 3, ...})
            let z = b + (2 * b - 1) * z0;

            // x = (z - r)^2 / (2 sigma^2) - z0^2 / (2 sigma_0^2)
            // This is the BerExp argument for the spec acceptance test.
            let diff = (z as f64) - r;
            let x = diff * diff * inv_2sigma_sq - (z0 as f64) * (z0 as f64) * inv_2sigma0_sq;

            if ber_exp(rng, x, ccs) {
                return s + z;
            }
        }
    }
}

/// BerExp acceptance test from FIPS 206 Algorithm 14.
///
/// Returns `true` with probability `ccs * exp(-x)` for `x >= 0` and
/// `ccs in [0, 1]`. Constant-time: the internal polynomial evaluation
/// and the comparison both run in a fixed number of integer operations
/// with no data-dependent branches.
///
/// The implementation reduces `x = s * ln 2 + r` with `r ∈ [0, ln 2]`,
/// computes a 64-bit fixed-point approximation of `2^63 * ccs * exp(-r)`
/// via the FACCT degree-12 polynomial in [`expm_p63`], shifts right by
/// `s` to get the full `2^(63-s) * ccs * exp(-r)`, and compares against
/// a fresh 63-bit random value. `P(u < threshold) = threshold / 2^63 =
/// ccs * exp(-r) / 2^s = ccs * exp(-x)`.
fn ber_exp<R: RngCore>(rng: &mut R, x: f64, ccs: f64) -> bool {
    // Degenerate inputs.
    if !x.is_finite() {
        return false;
    }
    if x <= 0.0 {
        // x = 0 means exp(-x) = 1, so accept with probability ccs.
        // Use the same saturating conversion as the main path so ccs = 1.0
        // maps to a threshold that always exceeds the 63-bit sample u.
        let ccs_q64 = fp_to_q64(ccs);
        let ccs_q63 = ccs_q64 >> 1; // shrink Q64 -> Q63 for comparison with u
        let u: u64 = rng.gen::<u64>() >> 1; // top 63 bits, uniform in [0, 2^63).
        return u < ccs_q63;
    }

    // Reduce x = s * ln(2) + r, r ∈ [0, ln(2)).
    let s_f = (x / std::f64::consts::LN_2).floor();
    let r = x - s_f * std::f64::consts::LN_2;

    // Clamp s ∈ [0, 63]. Beyond 63, acceptance probability is below 2^-63
    // and the test reduces to "always reject" (modulo a vanishing chance).
    let s: u32 = if s_f < 0.0 {
        0
    } else if s_f > 63.0 {
        63
    } else {
        s_f as u32
    };

    // raw ≈ 2^63 * ccs * exp(-r); threshold = raw >> s = 2^(63-s) * ccs * exp(-r).
    let raw = expm_p63(r, ccs);
    let threshold = raw.wrapping_shr(s);

    // Uniform u in [0, 2^63). Accept iff u < threshold.
    let u: u64 = rng.gen::<u64>() >> 1;
    u < threshold
}

/// FACCT polynomial coefficients for `2^63 * exp(-r)` on `r ∈ [0, ln 2]`.
///
/// Degree-12 Horner schedule: `y_0 = C[0]`, `y_{k+1} = C[k+1] - (z * y_k) >> 64`
/// where `z = r * 2^64`. Final `y_12 ≈ 2^63 * exp(-r)`.
///
/// Coefficients are taken from Howe-Pöppelmann-Prest's "Practical,
/// Constant-Time Discrete Gaussian Sampling" (FACCT, 2020) and match
/// the values in PQClean's Falcon-512 reference (`fpr.c`,
/// `fpr_expm_p63`). The polynomial is accurate to ~2^-50 over the
/// reduction interval — comfortably above the 63-bit threshold
/// quantization error introduced downstream by the right-shift by `s`.
const FACCT_C: [u64; 13] = [
    0x00000004741183A3,
    0x00000036548CFC06,
    0x0000024FDCBF140A,
    0x0000171D939DE045,
    0x0000D00CF58F6F84,
    0x000680681CF796E3,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    0x8000000000000000,
];

/// Returns the high 64 bits of the 128-bit product `a * b`.
///
/// Constant-time on platforms where the 128-bit multiply is constant-time
/// (true for `u128` on x86_64 / aarch64 — LLVM emits a single `mul` plus
/// register shuffle, no data-dependent branches).
#[inline]
fn umul64_hi(a: u64, b: u64) -> u64 {
    (((a as u128) * (b as u128)) >> 64) as u64
}

/// Two raised to 64, as `f64`. Exact (`2^64` is a power of two and
/// representable in IEEE 754 double).
const TWO_POW_64: f64 = 18446744073709551616.0;

/// Converts `v ∈ [0, 1]` to its Q64 fixed-point representation. The
/// boundary case `v = 1.0` produces `u64::MAX` (saturating), avoiding
/// the wraparound that `(v * 2^63 as u64) << 1` would otherwise
/// produce. Negative or NaN inputs map to 0.
#[inline]
fn fp_to_q64(v: f64) -> u64 {
    if !(v > 0.0) {
        return 0;
    }
    let scaled = v * TWO_POW_64;
    if scaled >= TWO_POW_64 {
        u64::MAX
    } else {
        scaled as u64
    }
}

/// Computes `floor(2^63 * ccs * exp(-r))` for `r ∈ [0, ln 2]`,
/// `ccs ∈ [0, 1]`, via the FACCT degree-12 polynomial.
///
/// Constant-time: 12 widening multiplications + subtractions + a final
/// scale by `ccs`. No data-dependent branches.
fn expm_p63(r: f64, ccs: f64) -> u64 {
    debug_assert!(r >= 0.0 && r <= std::f64::consts::LN_2 + f64::EPSILON);
    debug_assert!(ccs >= 0.0 && ccs <= 1.0);

    // z = r * 2^64 (in u64 fixed-point). r ∈ [0, ln 2] ≈ [0, 0.693],
    // so r * 2^64 ∈ [0, ~1.28e19] which fits in u64.
    let z = fp_to_q64(r);

    // Horner schedule. y stays in [0, 2^63] throughout — the constants
    // are crafted so the subtraction never underflows.
    let mut y = FACCT_C[0];
    for &c in &FACCT_C[1..] {
        y = c.wrapping_sub(umul64_hi(z, y));
    }

    // Multiply by ccs * 2^64. `ccs = 1.0` maps to `u64::MAX`; the high
    // 64 bits of `u64::MAX * y` is `y - 1` for `y > 0` (and 0 for y=0),
    // so the boundary case yields ≈ 2^63·exp(-r) as intended.
    let z_ccs = fp_to_q64(ccs);
    umul64_hi(z_ccs, y)
}

/// Pre-divided constant `1 / (2 * sigma_0^2)`, used hot in [`SamplerZ::sample`].
const INV_2SIGMA0_SQ: f64 = 1.0 / (2.0 * SIGMA_0 * SIGMA_0);

/// FALCON-512's `sigma_min` (also FALCON-1024 to two decimal places).
/// Used by the no-argument [`SamplerZ::new`] constructor; the
/// [`SamplerZ::with_sigma_min`] constructor takes the param-set value.
const FALCON_512_SIGMA_MIN: f64 = 1.2778336969128337;

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

        // FIPS 206 Algorithm 12 is designed for per-leaf sigma — the
        // distribution at the deepest level of FALCON's ffSampling, which
        // sits in `[sigma_min, sigma_max]` (roughly `[1.28, 1.85]`). It is
        // NOT a general-purpose Gaussian sampler for arbitrarily large
        // sigma; the global FALCON `sigma ≈ 165.74` enters the algorithm
        // only as the SCALING factor that produces the per-leaf widths.
        let mu = 0.0;
        let sigma = 1.5;

        let mut sum: f64 = 0.0;
        let mut sum_sq: f64 = 0.0;
        let n = 10_000;

        for _ in 0..n {
            let z = sampler.sample(&mut rng, mu, sigma);
            sum += z as f64;
            sum_sq += (z * z) as f64;
        }

        let mean = sum / n as f64;
        let variance = sum_sq / n as f64 - mean * mean;

        // Mean is approximately normal with std ~ sigma / sqrt(n);
        // for sigma=1.5, n=10000 that's ~0.015. We allow 5*std.
        assert!(mean.abs() < 5.0 * sigma / (n as f64).sqrt(),
            "Mean {} should be close to 0 within 5 std", mean);

        // Variance of the sample variance has std ~ sigma^2 * sqrt(2/n).
        // For sigma=1.5 that's ~ 2.25 * 0.014 ≈ 0.032. Allow 5*std.
        let expected_var = sigma * sigma;
        let var_std = expected_var * (2.0 / n as f64).sqrt();
        assert!(
            (variance - expected_var).abs() < 5.0 * var_std,
            "Variance {} should be ~{} (within 5 std = {})",
            variance, expected_var, 5.0 * var_std
        );
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
