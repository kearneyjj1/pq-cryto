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

/// 1 / (2 * sigma_0^2) for the BerExp rejection test.
/// Computed from SIGMA_0 to ensure internal consistency with the RCDT table.
const INV_2SIGMA0_SQ: f64 = 1.0 / (2.0 * SIGMA_0 * SIGMA_0);

/// Maximum deviation from the mean (in units of sigma) for test samplers.
/// Values beyond this are essentially zero probability.
#[cfg(test)]
const MAX_SIGMA_MULT: f64 = 10.0;

// ============================================================================
// Naive Gaussian Samplers (test-only)
// ============================================================================

/// Samples from a discrete Gaussian distribution N(0, sigma^2) over Z.
///
/// Uses rejection sampling. This is NOT constant-time and is only used
/// as a fallback for toy test parameters where sigma < sigma_min.
#[cfg(test)]
pub fn sample_z_gaussian<R: RngCore>(rng: &mut R, sigma: f64) -> i64 {
    if sigma <= 0.0 {
        return 0;
    }

    let sign: i64 = if rng.gen_bool(0.5) { 1 } else { -1 };
    let abs_z = sample_half_gaussian(rng, sigma);

    if abs_z == 0 {
        0
    } else {
        sign * abs_z
    }
}

/// Samples |z| from a half-Gaussian (z >= 0).
#[cfg(test)]
fn sample_half_gaussian<R: RngCore>(rng: &mut R, sigma: f64) -> i64 {
    let max_z = (sigma * MAX_SIGMA_MULT).ceil() as i64;

    loop {
        let z = rng.gen_range(0..=max_z);
        let prob = (-((z * z) as f64) / (2.0 * sigma * sigma)).exp();
        let u: f64 = rng.gen();
        if u < prob {
            return z;
        }
    }
}

/// Samples from a discrete Gaussian N(mu, sigma^2) centered at mu.
///
/// Uses rejection sampling. NOT constant-time — only for test builds.
#[cfg(test)]
pub fn sample_gaussian<R: RngCore>(rng: &mut R, mu: f64, sigma: f64) -> i64 {
    if sigma <= 0.0 {
        return mu.round() as i64;
    }

    let max_offset = (sigma * MAX_SIGMA_MULT).ceil() as i64;
    let mu_floor = mu.floor() as i64;

    loop {
        let z = mu_floor + rng.gen_range(-max_offset..=max_offset);
        let diff = z as f64 - mu;
        let prob = (-(diff * diff) / (2.0 * sigma * sigma)).exp();
        let u: f64 = rng.gen();
        if u < prob {
            return z;
        }
    }
}

// ============================================================================
// Keygen Gaussian Sampler
// ============================================================================

/// Samples from a discrete Gaussian N(0, sigma^2) for key generation.
///
/// Uses fixed-iteration rejection sampling with uniform proposal over [-4σ, 4σ].
/// All iterations consume the same RNG bytes and use branchless integer-based
/// acceptance to avoid timing side-channels on the sampled value.
///
/// The acceptance test uses `fpr_expm_p63` (Q63 integer polynomial) instead of
/// `f64::exp()` to avoid variable-time libm calls. The comparison uses unsigned
/// integer subtraction instead of floating-point comparison to avoid branches.
///
/// The iteration count is computed from the acceptance rate to guarantee failure
/// probability < 2^-128. Since sigma is derived from the public parameter n,
/// the iteration count is public and leaks no secret information.
///
/// The FIPS 206 SamplerZ (Algorithm 12) cannot be used for keygen because its
/// base sampler at sigma_0=1.82 has lighter tails than the target sigma≈4.05,
/// causing BerExp's acceptance probability to cap at 1.0 for large |z| and
/// producing incorrect (too-low) variance.
pub fn sample_keygen_gaussian<R: RngCore>(rng: &mut R, sigma: f64) -> i64 {
    if sigma <= 0.0 {
        return 0;
    }

    // Defense-in-depth: clamp sigma to a minimum of 0.5 to ensure the
    // dynamically computed iteration count always achieves < 2^-128 failure.
    // Keygen sigma is always >= 2.86 (for n=1024), so this never activates.
    let sigma = sigma.max(0.5);

    // 4σ tail bound: captures 99.994% of the Gaussian mass.
    // Variance bias from truncation is < 0.11%, negligible for keygen.
    let max_z = (sigma * 4.0).ceil() as i64;
    let two_sigma_sq = 2.0 * sigma * sigma;
    let range = (2 * max_z + 1) as u64;

    // Compute iterations for failure probability < 2^-128.
    // Acceptance rate ≈ sigma * sqrt(2π) / range ≈ 0.29–0.31 for typical keygen sigma.
    // Iterations = ceil(128 * ln(2) / -ln(1 - rate)).
    let acceptance_rate = (sigma * std::f64::consts::TAU.sqrt() / range as f64).min(0.99);
    let fixed_iters = (128.0_f64 * 2.0_f64.ln() / -(1.0 - acceptance_rate).ln())
        .ceil() as u32;
    let fixed_iters = fixed_iters.max(64);

    // Constants for decomposing x = s*ln(2) + r
    const LN2: f64 = 0.693_147_180_559_945_309_417_232_121_458;
    const INV_LN2: f64 = 1.442_695_040_888_963_407_359_924_681_002;

    let mut result = 0i64;
    let mut accepted = 0u64; // 0 = not yet, 1 = accepted

    for _ in 0..fixed_iters {
        // Constant-time uniform in [-max_z, max_z] via modular reduction.
        // Bias from u64 % range is < range / 2^64, negligible.
        let raw: u64 = rng.next_u64();
        let z = (raw % range) as i64 - max_z;

        // Integer-based acceptance: compute exp(-z^2/(2*sigma^2)) using Q63
        // polynomial, avoiding variable-time f64::exp() and float branches.
        // Decompose x = z^2/(2*sigma^2) into s*ln(2) + r, r in [0, ln(2)).
        let x = ((z * z) as f64) / two_sigma_sq;
        let s = (x * INV_LN2).floor() as u32;
        let r = x - (s as f64) * LN2;

        // z_q63 = floor(2^63 * exp(-r)) via integer polynomial
        let mut z_q63 = fpr_expm_p63(r, 1.0);

        // Branchless halving s times: z_q63 *= 2^(-s)
        // For keygen with 4σ tail, x <= 8 so s <= 12 — well within 63.
        let s_halve = s.min(63);
        for i in 0..63u32 {
            let is_active = ((i.wrapping_sub(s_halve)) >> 31) & 1;
            let mask = (is_active as u64).wrapping_neg();
            let halved = (z_q63 + 1) >> 1;
            z_q63 ^= (z_q63 ^ halved) & mask;
        }

        // Draw 63 random bits and accept via integer comparison (constant-time).
        // Unsigned subtraction: high bit set iff w < z_q63 (borrow).
        let w = rng.next_u64() >> 1;
        let accept = (w.wrapping_sub(z_q63)) >> 63; // 1 if w < z_q63

        let first_accept = accept & (1 - accepted);
        // Branchless select: update result only on first acceptance
        let mask = first_accept.wrapping_neg() as i64; // 0 or -1 (all-ones)
        result = (z & mask) | (result & !mask);
        accepted |= accept;
    }

    // Statistically unreachable: sigma >= 0.5 guarantees enough iterations
    // for < 2^-128 failure probability.
    assert!(accepted != 0, "keygen Gaussian sampling failed after {} iterations", fixed_iters);
    result
}

// ============================================================================
// Bernoulli Sampler (old floating-point version, test-only)
// ============================================================================

/// Samples a Bernoulli variable: returns true with probability exp(-x) for x >= 0.
///
/// Uses floating-point arithmetic. NOT constant-time — test builds only.
#[cfg(test)]
pub fn sample_bernoulli_exp<R: RngCore>(rng: &mut R, x: f64) -> bool {
    if x <= 0.0 {
        return true;
    }
    if x >= 1.0 {
        if !sample_bernoulli_exp_minus_one(rng) {
            return false;
        }
        return sample_bernoulli_exp(rng, x - 1.0);
    }
    let prob = (-x).exp();
    rng.gen::<f64>() < prob
}

/// Samples Bernoulli with probability exp(-1) ≈ 0.368 (test-only).
#[cfg(test)]
fn sample_bernoulli_exp_minus_one<R: RngCore>(rng: &mut R) -> bool {
    const EXP_MINUS_1: f64 = 0.36787944117144233;
    rng.gen::<f64>() < EXP_MINUS_1
}

// ============================================================================
// FIPS 206 RCDT Table and BaseSampler
// ============================================================================

/// FIPS 206 RCDT (Reverse Cumulative Distribution Table) for the half-Gaussian
/// at sigma_0 = 1.8205.
///
/// 18 entries, each stored as 3 × 24-bit limbs `[w2, w1, w0]`.
/// Full 72-bit value = `w2 * 2^48 + w1 * 2^24 + w0`.
/// Entry `i` represents `floor(2^72 * P(|z| >= i+1))`.
///
/// From the Falcon reference implementation (`sign.c` `dist[]` array).
const RCDT: [[u32; 3]; 18] = [
    [10745844, 3068844, 3741698],
    [5559083, 1580863, 8248194],
    [2260429, 13669192, 2736639],
    [708981, 4421575, 10046180],
    [169348, 7122675, 4136815],
    [30538, 13063405, 7650655],
    [4132, 14505003, 7826148],
    [417, 16768101, 11363290],
    [31, 8444042, 8086568],
    [1, 12844466, 265321],
    [0, 1232676, 13644283],
    [0, 38047, 9111839],
    [0, 870, 6138264],
    [0, 14, 12545723],
    [0, 0, 3104126],
    [0, 0, 28824],
    [0, 0, 198],
    [0, 0, 1],
];

/// The base sampler for FALCON's signature scheme (FIPS 206).
///
/// Samples from the half-Gaussian distribution at `sigma_0 = 1.8205`
/// using the RCDT method with 72-bit precision. The comparison is
/// constant-time using multi-limb subtraction with borrow propagation.
pub struct BaseSampler;

impl BaseSampler {
    /// Creates a new base sampler.
    pub fn new() -> Self {
        BaseSampler
    }

    /// Samples from the full discrete Gaussian at sigma_0.
    ///
    /// Returns a signed integer `z` with `P(z) ∝ exp(-z²/(2σ₀²))`.
    /// Uses 72-bit precision RCDT with constant-time comparison.
    #[cfg(test)]
    pub fn sample<R: RngCore>(&self, rng: &mut R) -> i64 {
        let z = self.sample_half(rng);

        // Random sign bit (constant-time)
        let mut sign_byte = [0u8; 1];
        rng.fill_bytes(&mut sign_byte);
        let b = (sign_byte[0] & 1) as i64;
        let sign = 1 - 2 * b; // +1 or -1

        // z=0 * any sign = 0, which is correct
        sign * z
    }

    /// Samples `z0 >= 0` from the half-Gaussian at `sigma_0 = 1.8205`.
    ///
    /// Draws 9 random bytes (72 bits), splits into 3 × 24-bit limbs,
    /// and performs a constant-time scan against the RCDT table.
    /// The result `z0` counts how many RCDT entries exceed the random value.
    pub fn sample_half<R: RngCore>(&self, rng: &mut R) -> i64 {
        // Draw 72 random bits
        let mut bytes = [0u8; 9];
        rng.fill_bytes(&mut bytes);

        // Split into 3 × 24-bit limbs (little-endian byte order)
        let v0 = (bytes[0] as u32) | ((bytes[1] as u32) << 8) | ((bytes[2] as u32) << 16);
        let v1 = (bytes[3] as u32) | ((bytes[4] as u32) << 8) | ((bytes[5] as u32) << 16);
        let v2 = (bytes[6] as u32) | ((bytes[7] as u32) << 8) | ((bytes[8] as u32) << 16);

        // Constant-time scan: count RCDT entries where (v2,v1,v0) < (w2,w1,w0).
        // Multi-limb subtraction: borrow propagates through bit 31 of each u32.
        let mut z: i64 = 0;
        for entry in &RCDT {
            let w2 = entry[0];
            let w1 = entry[1];
            let w0 = entry[2];
            let cc0 = (v0.wrapping_sub(w0)) >> 31;
            let cc1 = (v1.wrapping_sub(w1).wrapping_sub(cc0)) >> 31;
            let cc2 = (v2.wrapping_sub(w2).wrapping_sub(cc1)) >> 31;
            z += cc2 as i64;
        }
        z
    }
}

impl Default for BaseSampler {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Constant-time BerExp: Bernoulli exponential acceptance test
// ============================================================================

/// Q63 fixed-point coefficients for the exp(-x) Taylor polynomial.
/// `C[k] = floor(2^63 / k!)` for k = 0..12.
/// Used by the integer-only `fpr_expm_p63` implementation.
const EXP_COEFFS: [u64; 13] = [
    0x8000000000000000,  // C[ 0] = floor(2^63 / 0!)  = 9223372036854775808
    0x8000000000000000,  // C[ 1] = floor(2^63 / 1!)  = 9223372036854775808
    0x4000000000000000,  // C[ 2] = floor(2^63 / 2!)  = 4611686018427387904
    0x1555555555555555,  // C[ 3] = floor(2^63 / 3!)  = 1537228672809129301
    0x0555555555555555,  // C[ 4] = floor(2^63 / 4!)  = 384307168202282325
    0x0111111111111111,  // C[ 5] = floor(2^63 / 5!)  = 76861433640456465
    0x002D82D82D82D82D,  // C[ 6] = floor(2^63 / 6!)  = 12810238940076077
    0x0006806806806806,  // C[ 7] = floor(2^63 / 7!)  = 1830034134296582
    0x0000D00D00D00D00,  // C[ 8] = floor(2^63 / 8!)  = 228754266787072
    0x0000171DE3A556C7,  // C[ 9] = floor(2^63 / 9!)  = 25417140754119
    0x0000024FC9F6EF13,  // C[10] = floor(2^63 / 10!) = 2541714075411
    0x00000035CC8ACFEA,  // C[11] = floor(2^63 / 11!) = 231064915946
    0x000000047BB63BFE,  // C[12] = floor(2^63 / 12!) = 19255409662
];

/// Computes `(a * b) >> 63` using 128-bit intermediate arithmetic.
/// Both `a` and `b` are Q63 fixed-point values.
#[inline]
fn mulhi63(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) >> 63) as u64
}

/// Computes `floor(2^63 * ccs * exp(-x))` for `0 <= x < ln(2)`.
///
/// Uses pure integer fixed-point arithmetic (Q63 format) for the polynomial
/// evaluation, eliminating all floating-point operations in the inner loop.
/// Only two f64-to-integer conversions occur at the boundary (x and ccs).
///
/// The polynomial is a degree-12 Taylor approximation of exp(-x),
/// evaluated via Horner's method in Q63 fixed-point:
///   y = C[12]
///   y = C[k] - (x_q63 * y) >> 63   for k = 11 down to 0
///
/// Based on the Falcon reference implementation's `fpr_expm_p63` in `sign.c`.
fn fpr_expm_p63(x: f64, ccs: f64) -> u64 {
    // Convert x to Q63 fixed-point: x is in [0, ln(2)) ≈ [0, 0.693)
    // so x * 2^63 < 2^63, safe for u64.
    let x_q63 = (x * 9_223_372_036_854_775_808.0) as u64;

    // Convert ccs to Q63 fixed-point: ccs is in (0, 1] so ccs * 2^63 fits in u64.
    let ccs_q63 = (ccs * 9_223_372_036_854_775_808.0) as u64;

    // Horner evaluation in Q63 fixed-point, from highest to lowest coefficient.
    let mut y = EXP_COEFFS[12];
    y = EXP_COEFFS[11].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[10].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[9].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[8].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[7].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[6].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[5].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[4].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[3].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[2].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[1].wrapping_sub(mulhi63(x_q63, y));
    y = EXP_COEFFS[0].wrapping_sub(mulhi63(x_q63, y));

    // Final: multiply exp(-x) result by ccs (also in Q63)
    mulhi63(y, ccs_q63)
}

/// Original floating-point fpr_expm_p63, retained for comparison testing.
#[cfg(test)]
fn fpr_expm_p63_float(x: f64, ccs: f64) -> u64 {
    let mut y: f64 = 1.0 / 479_001_600.0;
    y = 1.0 / 39_916_800.0 - x * y;
    y = 1.0 / 3_628_800.0 - x * y;
    y = 1.0 / 362_880.0 - x * y;
    y = 1.0 / 40_320.0 - x * y;
    y = 1.0 / 5_040.0 - x * y;
    y = 1.0 / 720.0 - x * y;
    y = 1.0 / 120.0 - x * y;
    y = 1.0 / 24.0 - x * y;
    y = (1.0 / 6.0) - x * y;
    y = 0.5 - x * y;
    y = 1.0 - x * y;
    y = 1.0 - x * y;
    let z = y * ccs * 9_223_372_036_854_775_808.0;
    z as u64
}

/// Accepts with probability `ccs * exp(-x)` using branchless control flow.
///
/// Implements the Falcon reference BerExp algorithm with constant-time handling
/// of both positive and negative x:
///
/// 1. Decompose `x = s * ln(2) + r` where `0 <= r < ln(2)` (floor handles all x)
/// 2. Compute `z = floor(2^63 * ccs * exp(-r))` via integer polynomial
/// 3. For `s >= 0`: halve `z` exactly `s` times (constant-time masked loop)
/// 4. For `s < 0`: probability exceeds representable range, force accept via mask
/// 5. Accept if a random 63-bit value is less than `z`
///
/// The negative-x case (where `exp(-x) > 1`) occurs when sigma slightly exceeds
/// sigma_0 during signing. Rather than branching to floating-point `exp()`, we
/// use arithmetic masking to set z to its maximum value, ensuring acceptance.
fn ber_exp<R: RngCore>(rng: &mut R, x: f64, ccs: f64) -> bool {
    const LN2: f64 = 0.693_147_180_559_945_309_417_232_121_458;
    const INV_LN2: f64 = 1.442_695_040_888_963_407_359_924_681_002;

    // Decompose x = s * ln(2) + r where 0 <= r < ln(2).
    // floor() works for both positive and negative x.
    let s_float = (x * INV_LN2).floor();
    let r = x - s_float * LN2;
    let s = s_float as i64;

    // Compute z = floor(2^63 * ccs * exp(-r)) using integer polynomial.
    // r is always in [0, ln(2)) regardless of the sign of x.
    let mut z = fpr_expm_p63(r, ccs);

    // For s >= 0: halve z exactly s times.
    // exp(-x) = exp(-r) * 2^(-s), so each halving divides by 2.
    let s_halve = s.max(0).min(63) as u32;
    for i in 0..63u32 {
        let is_active = ((i.wrapping_sub(s_halve)) >> 31) & 1;
        let mask = (is_active as u64).wrapping_neg();
        let halved = (z + 1) >> 1;
        z ^= (z ^ halved) & mask;
    }

    // For s < 0: double z exactly |s| times (with saturation at 2^63-1).
    // exp(-x) = exp(-r) * 2^|s|, so each doubling multiplies by 2.
    // s_double = 0 when s >= 0 (loop is a no-op), |s| when s < 0.
    let s_double = ((-s).max(0).min(63)) as u32;
    for i in 0..63u32 {
        let is_active = ((i.wrapping_sub(s_double)) >> 31) & 1;
        let mask = (is_active as u64).wrapping_neg();
        // Saturating double: if z >= 2^62, cap at 2^63-1
        let cap = (z >> 62) & 1;
        let cap_mask = cap.wrapping_neg();
        let doubled = (z.wrapping_shl(1) & !cap_mask)
            | (0x7FFF_FFFF_FFFF_FFFF & cap_mask);
        z ^= (z ^ doubled) & mask;
    }

    // Draw random 63-bit value and compare
    let w: u64 = rng.gen::<u64>() >> 1;
    w < z
}

/// Original branching ber_exp for comparison testing.
#[cfg(test)]
fn ber_exp_branching<R: RngCore>(rng: &mut R, x: f64, ccs: f64) -> bool {
    const LN2: f64 = 0.693_147_180_559_945_309_417_232_121_458;
    const INV_LN2: f64 = 1.442_695_040_888_963_407_359_924_681_002;

    if x < 0.0 {
        let prob = ccs * (-x).exp();
        if prob >= 1.0 {
            let _w: u64 = rng.gen();
            return true;
        }
        let z = (prob * 9_223_372_036_854_775_808.0) as u64;
        let w: u64 = rng.gen::<u64>() >> 1;
        return w < z;
    }

    let s = (x * INV_LN2).floor() as u32;
    let r = x - (s as f64) * LN2;
    let mut z = fpr_expm_p63(r, ccs);
    for i in 0..63u32 {
        let is_active = ((i.wrapping_sub(s)) >> 31) & 1;
        let mask = (is_active as u64).wrapping_neg();
        let halved = (z + 1) >> 1;
        z ^= (z ^ halved) & mask;
    }
    let w: u64 = rng.gen::<u64>() >> 1;
    w < z
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
    /// Requires `sigma >= sigma_min`. In test builds, falls back to naive rejection
    /// sampling for toy parameters where sigma < sigma_min. In release builds,
    /// sigma is clamped to sigma_min with a debug_assert.
    pub fn sample<R: RngCore>(&self, rng: &mut R, mu: f64, sigma: f64) -> i64 {
        // Guard: sigma must be >= sigma_min for FIPS 206 correctness
        #[cfg(test)]
        if sigma < self.sigma_min {
            return sample_gaussian(rng, mu, sigma.max(0.01));
        }
        #[cfg(not(test))]
        let sigma = {
            debug_assert!(
                sigma >= self.sigma_min,
                "SamplerZ requires sigma ({}) >= sigma_min ({})",
                sigma, self.sigma_min
            );
            sigma.max(self.sigma_min)
        };

        let s = mu.floor() as i64;
        let r = mu - (s as f64);
        let dss = 1.0 / (2.0 * sigma * sigma);
        let ccs = self.sigma_min / sigma;

        // Note: this rejection loop is inherently variable-time (iteration count
        // depends on secret mu/sigma). This is fundamental to FALCON's design —
        // the reference implementation has the same property. The number of loop
        // iterations is considered public information in the security model.
        loop {
            // 1. Sample z0 >= 0 from half-Gaussian at sigma_0
            let z0 = self.base.sample_half(rng);

            // 2. Random sign bit (constant-time: avoid gen_bool which uses f64)
            let mut sign_byte = [0u8; 1];
            rng.fill_bytes(&mut sign_byte);
            let b: i64 = (sign_byte[0] & 1) as i64;

            // 3. z = b + (2b - 1) * z0
            //    b=1 → z = 1 + z0 (positive side)
            //    b=0 → z = -z0    (negative side)
            let z = b + (2 * b - 1) * z0;

            // 4. Compute exponent for rejection test
            //    x = (z - r)^2 / (2*sigma^2) - z0^2 / (2*sigma_0^2)
            let zr = (z as f64) - r;
            let x = zr * zr * dss - (z0 as f64) * (z0 as f64) * INV_2SIGMA0_SQ;

            // 5. BerExp acceptance test
            let accept_ber = ber_exp(rng, x, ccs);

            // 6. Constant-time rejection of (z=0, b=0) case.
            // z_nonzero_or_b: is 1 if z != 0 OR b != 0, i.e. the sample is valid.
            // For z as i64: z != 0 iff (z | z.wrapping_neg()) has bit 63 set.
            let z_nonzero = ((z as u64) | (z as u64).wrapping_neg()) >> 63;
            let valid = z_nonzero | (b as u64);

            // Accept iff BerExp accepted AND (z,b) is not the forbidden (0,0) case.
            // Use bitwise & (not short-circuit &&) to avoid timing leak on accept_ber.
            if accept_ber & (valid != 0) {
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

    #[test]
    fn test_ber_exp_negative_x() {
        // Verify ber_exp handles negative x correctly (ccs * exp(-x) < 1).
        // This is critical for SamplerZ when sigma > sigma_0 (keygen case).
        let mut rng = StdRng::seed_from_u64(222);
        let ccs: f64 = SIGMA_MIN / 4.05; // ≈ 0.315

        // For x = -0.5, correct probability = ccs * exp(0.5) ≈ 0.519
        let mut count = 0;
        let n = 20000;
        for _ in 0..n {
            if ber_exp(&mut rng, -0.5, ccs) {
                count += 1;
            }
        }
        let rate = (count as f64) / (n as f64);
        let expected = ccs * (0.5f64).exp();
        assert!(
            (rate - expected).abs() < 0.03,
            "BerExp(-0.5, {}): rate {} should be close to {}",
            ccs, rate, expected
        );

        // For x = -0.03, correct probability = ccs * exp(0.03) ≈ 0.325
        let mut count2 = 0;
        for _ in 0..n {
            if ber_exp(&mut rng, -0.03, ccs) {
                count2 += 1;
            }
        }
        let rate2 = (count2 as f64) / (n as f64);
        let expected2 = ccs * (0.03f64).exp();
        assert!(
            (rate2 - expected2).abs() < 0.03,
            "BerExp(-0.03, {}): rate {} should be close to {}",
            ccs, rate2, expected2
        );
    }

    #[test]
    fn test_fpr_expm_p63_integer_vs_float() {
        // Compare integer and floating-point implementations across [0, ln(2))
        let ccs_values = [1.0, 0.5, 0.315, 0.999, 0.01];
        for &ccs in &ccs_values {
            for i in 0..500 {
                let x = (i as f64) * 0.693 / 500.0;
                let int_result = fpr_expm_p63(x, ccs);
                let float_result = fpr_expm_p63_float(x, ccs);
                let diff = (int_result as i128 - float_result as i128).unsigned_abs();
                assert!(
                    diff <= 2048, // Allow small rounding differences from Q63 truncation
                    "Mismatch at x={:.6}, ccs={}: int={}, float={}, diff={}",
                    x, ccs, int_result, float_result, diff
                );
            }
        }
    }

    #[test]
    fn test_fpr_expm_p63_boundary_values() {
        // x = 0: exp(0) = 1, result should be ccs * 2^63
        let r = fpr_expm_p63(0.0, 1.0);
        assert!(r > 0x7FFFFFFFFF000000u64, "exp(0)*1.0*2^63 should be near 2^63, got {}", r);

        // x near ln(2): exp(-0.693) ≈ 0.50028, result should be near that * 2^63
        let r2 = fpr_expm_p63(0.693, 1.0);
        let expected = ((-0.693f64).exp() * 9_223_372_036_854_775_808.0) as u64;
        let diff = (r2 as i128 - expected as i128).unsigned_abs();
        assert!(diff < 1_000_000_000, "exp(-0.693) result too far from expected, diff={}", diff);
    }

    #[test]
    fn test_keygen_gaussian_variance() {
        // Verify sample_keygen_gaussian produces correct variance for sigma=4.05
        // (keygen sigma for FALCON-512).
        let mut rng = StdRng::seed_from_u64(42);
        let sigma = 4.05;

        let mut sum: i64 = 0;
        let mut sum_sq: i64 = 0;
        let n = 10000;

        for _ in 0..n {
            let z = sample_keygen_gaussian(&mut rng, sigma);
            sum += z;
            sum_sq += z * z;
        }

        let mean = (sum as f64) / (n as f64);
        let variance = (sum_sq as f64) / (n as f64) - mean * mean;
        let expected_var = sigma * sigma; // 16.4

        assert!(
            mean.abs() < 0.5,
            "Mean {} should be close to 0",
            mean
        );

        // Variance should be within 15% of sigma^2
        assert!(
            variance > expected_var * 0.85 && variance < expected_var * 1.15,
            "Variance {} should be close to {} (sigma={})",
            variance, expected_var, sigma
        );
    }

    #[test]
    fn test_ber_exp_branchless_vs_branching() {
        // Compare branchless ber_exp against the original branching version
        // using the same RNG seeds. Both should produce identical results.
        let ccs = 0.315;
        let x_values = [0.0, 0.1, 0.5, 1.0, 2.0, -0.01, -0.1, -0.5, -1.0];

        for &x in &x_values {
            let n = 5000;
            let mut rng1 = StdRng::seed_from_u64(777);
            let mut rng2 = StdRng::seed_from_u64(777);

            let mut count_bl = 0;
            let mut count_br = 0;
            for _ in 0..n {
                if ber_exp(&mut rng1, x, ccs) { count_bl += 1; }
                if ber_exp_branching(&mut rng2, x, ccs) { count_br += 1; }
            }

            let rate_bl = count_bl as f64 / n as f64;
            let rate_br = count_br as f64 / n as f64;

            // Rates should be statistically close (within 3%)
            assert!(
                (rate_bl - rate_br).abs() < 0.03,
                "BerExp mismatch at x={}: branchless={:.3}, branching={:.3}",
                x, rate_bl, rate_br
            );
        }
    }
}
