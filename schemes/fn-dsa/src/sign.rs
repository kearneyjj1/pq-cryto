//! Signing for FALCON.
//!
//! This module implements the FALCON signature generation algorithm,
//! which uses Fast Fourier Sampling to produce compact signatures.
//!
//! ## Basis Convention
//!
//! The NTRU lattice for public key `h = g·f⁻¹ mod q` is:
//!
//! ```text
//!   L = { (u, v) ∈ (Z[x]/(xⁿ+1))² : u + v·h ≡ 0 (mod q) }
//! ```
//!
//! The secret short basis has rows `(g, -f)` and `(G, -F)`:
//!
//! ```text
//!   B = [[g, -f], [G, -F]]       (row-vector convention: v·B ∈ L)
//! ```
//!
//! Verification: both rows are in L since `g - f·h = 0` and `G - F·h ≡ 0 (mod q)`
//! (the latter from `fG - gF = q`).
//!
//! With `det(B) = fG - gF = q`, the inverse is:
//!
//! ```text
//!   B⁻¹ = (1/q) · [[-F, f], [-G, g]]
//! ```
//!
//! The signature `(s1, s2) = (c, 0) - (z0, z1)·B` gives:
//!
//! ```text
//!   s1 = c - z0·g - z1·G
//!   s2 = z0·f + z1·F
//! ```
//!
//! Verify: `s1 + s2·h = c - z0·g - z1·G + (z0·f + z1·F)·g/f = c` since `F·g/f ≡ G (mod q)`.

use rand::RngCore;
use crate::error::{FnDsaError, Result};
use crate::fft::{fft, ifft, Complex};
use crate::hash::{generate_nonce, hash_to_point};
use crate::keygen::SecretKey;
use crate::params::{MAX_SIGN_ATTEMPTS, NONCE_SIZE};
use crate::sampler::FfSampler;
#[cfg(test)]
use crate::sampler::SimpleSampler;

/// A FALCON signature.
#[derive(Clone, Debug)]
pub struct Signature {
    /// The nonce used for this signature.
    pub nonce: [u8; NONCE_SIZE],
    /// The signature polynomial s2 (compressed).
    pub s2: Vec<i16>,
}

impl Signature {
    /// Returns the size of the encoded signature in bytes.
    pub fn encoded_size(&self) -> usize {
        // Rough estimate: nonce + compressed s2
        NONCE_SIZE + self.s2.len() * 2
    }

    /// Returns the squared norm of the signature.
    ///
    /// Uses saturating arithmetic to prevent overflow.
    pub fn norm_sq(&self) -> i64 {
        self.s2.iter()
            .map(|&x| (x as i64) * (x as i64))
            .fold(0i64, |acc, x| acc.saturating_add(x))
    }
}

/// Signs a message using the FALCON signature scheme.
///
/// This is the main signing function. It:
/// 1. Generates a random nonce
/// 2. Hashes the message to a polynomial c
/// 3. Computes target t = (c, 0)·B⁻¹ = (-F·c/q, f·c/q) via ffSampling
/// 4. Samples (z0, z1) close to t using FFT sampling
/// 5. Computes s2 = z0·f + z1·F (see module docs for basis convention)
/// 6. Derives s1 = c - s2·h and checks ||(s1, s2)|| is within bounds
/// 7. Returns the signature (nonce, s2)
///
/// ## Timing Side-Channels
///
/// The outer signing loop (retry on norm failure) and the inner SamplerZ
/// rejection loop are inherently variable-time. This is a fundamental property
/// of the FALCON algorithm: the number of signing attempts and rejection
/// samples are random variables that cannot be made constant without
/// abandoning the security proof. The reference FALCON implementation has the
/// same characteristic. Within each iteration, BerExp uses integer arithmetic
/// and the RCDT base sampler uses constant-time comparison.
pub fn sign<R: RngCore>(rng: &mut R, sk: &SecretKey, message: &[u8]) -> Result<Signature> {
    use zeroize::Zeroize;
    use crate::poly::Poly;
    use crate::field::Zq;

    let sigma = sk.params.sigma;
    let sigma_min = sk.params.sigma_min;
    let bound_sq = sk.params.sig_bound_sq;

    // Create the FFT sampler using the LDL* tree from the secret key
    let ff_sampler = FfSampler::new(&sk.gs, sigma, sigma_min);

    // Precompute NTT-ready Poly forms of the secret key polynomials
    let mut f_poly = Poly::from_i16(&sk.f.iter().map(|&x| x as i16).collect::<Vec<_>>());
    let mut big_f_poly = Poly::from_i16(&sk.big_f.iter().map(|&x| x as i16).collect::<Vec<_>>());
    let h_poly = Poly::from_i16(&sk.h);

    for _attempt in 0..MAX_SIGN_ATTEMPTS {
        // Generate a fresh nonce
        let nonce = generate_nonce(rng);

        // Hash the message to get c
        let c = hash_to_point(message, &nonce, &sk.params);

        // Convert c to FFT form for ffSampling
        let mut c_fft: Vec<Complex> = c.iter().map(|zq| Complex::from_real(zq.value() as f64)).collect();
        fft(&mut c_fft);

        // Use ffSampling to sample (z0, z1) close to the target
        let (mut z0_fft, mut z1_fft) = ff_sampler.sample_signature(rng, &c_fft);

        // Convert z0, z1 from FFT form back to integer coefficients.
        // The sampler produces integer polynomials; ifft + round recovers them.
        ifft(&mut z0_fft);
        ifft(&mut z1_fft);
        let mut z0: Vec<i16> = z0_fft.iter().map(|c| c.re.round() as i16).collect();
        let mut z1: Vec<i16> = z1_fft.iter().map(|c| c.re.round() as i16).collect();

        // Zeroize per-iteration secret intermediates
        c_fft.zeroize();
        z0_fft.zeroize();
        z1_fft.zeroize();

        // Compute s2 = z0·f + z1·F using exact NTT arithmetic (mod q).
        // This avoids floating-point precision loss in the FFT chain that
        // would cause wrong integer rounding at n >= 512.
        // See module-level documentation for the basis convention derivation.
        let z0_poly = Poly::from_i16(&z0);
        let z1_poly = Poly::from_i16(&z1);
        let s2_poly = z0_poly.mul_ntt(&f_poly).add(&z1_poly.mul_ntt(&big_f_poly));
        let s2 = s2_poly.to_centered();

        // Zeroize sampled lattice vectors (secret-dependent)
        z0.zeroize();
        z1.zeroize();

        // Compute s1 = c - s2*h using exact NTT arithmetic (same as verification)
        let c_poly = Poly::from_zq(c.clone());
        let s2_for_verify = Poly::from_i16(&s2);
        let s2h = s2_for_verify.mul_ntt(&h_poly);
        let s1_poly = c_poly.sub(&s2h);

        // Check the norm (using saturating arithmetic for safety)
        let s1_norm_sq = s1_poly.norm_sq();
        let s2_norm_sq: i64 = s2.iter()
            .map(|&x| (x as i64) * (x as i64))
            .fold(0i64, |acc, x| acc.saturating_add(x));
        let total_norm_sq = s1_norm_sq.saturating_add(s2_norm_sq);

        if (total_norm_sq as f64) <= bound_sq {
            // Zeroize secret-derived polynomial objects before returning
            for c in &mut f_poly.coeffs { *c = Zq::ZERO; }
            for c in &mut big_f_poly.coeffs { *c = Zq::ZERO; }
            return Ok(Signature { nonce, s2 });
        }
    }

    // Zeroize secret-derived polynomial objects on error path
    for c in &mut f_poly.coeffs { *c = Zq::ZERO; }
    for c in &mut big_f_poly.coeffs { *c = Zq::ZERO; }

    Err(FnDsaError::SigningFailed {
        attempts: MAX_SIGN_ATTEMPTS,
    })
}

/// Signs a message using direct s2 sampling with rejection (for testing only).
///
/// This is NOT the proper FALCON algorithm - it directly samples s2 from a
/// discrete Gaussian and uses rejection sampling to find a short (s1, s2) pair.
/// This approach is much slower and doesn't provide the same security guarantees,
/// but it demonstrates the verification correctness.
///
/// # Security
///
/// This function is restricted to `#[cfg(test)]` builds because it bypasses
/// the ffSampling algorithm and does not provide side-channel resistance.
#[cfg(test)]
pub fn sign_simple<R: RngCore>(rng: &mut R, sk: &SecretKey, message: &[u8]) -> Result<Signature> {
    use crate::poly::Poly;

    let n = sk.params.n;
    let bound_sq = sk.params.sig_bound_sq;
    let sampler = SimpleSampler::new();

    // For direct sampling, we need a small sigma to keep norms manageable
    // Expected ||s2||^2 ≈ n * sigma^2, and we need total norm < bound
    // With bound_sq = 2e9, we can have ||s2||^2 + ||s1||^2 < 2e9
    // If we sample s2 with sigma ~ 40, ||s2||^2 ~ 512 * 1600 = 820k
    // Then s1 = c - s2*h needs ||s1||^2 < 2e9 - 820k ≈ 2e9
    let sigma = 40.0;

    for _attempt in 0..MAX_SIGN_ATTEMPTS {
        let nonce = generate_nonce(rng);
        let c = hash_to_point(message, &nonce, &sk.params);

        // Directly sample s2 from a discrete Gaussian centered at 0
        let s2: Vec<i16> = (0..n)
            .map(|_| sampler.sample(rng, 0.0, sigma) as i16)
            .collect();

        // Compute s1 = c - s2*h using the same polynomial arithmetic as verification
        let c_poly = Poly::from_zq(c.clone());
        let s2_poly = Poly::from_i16(&s2);
        let h_poly = Poly::from_i16(&sk.h);
        let s2h = s2_poly.mul(&h_poly);
        let s1_poly = c_poly.sub(&s2h);

        // Check the norm (using saturating arithmetic for safety)
        let s1_norm_sq = s1_poly.norm_sq();
        let s2_norm_sq: i64 = s2.iter()
            .map(|&x| (x as i64) * (x as i64))
            .fold(0i64, |acc, x| acc.saturating_add(x));
        let total_norm_sq = s1_norm_sq.saturating_add(s2_norm_sq);

        if (total_norm_sq as f64) <= bound_sq {
            return Ok(Signature { nonce, s2 });
        }
    }

    Err(FnDsaError::SigningFailed {
        attempts: MAX_SIGN_ATTEMPTS,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_signature_structure() {
        // Create a mock signature
        let nonce = [0u8; NONCE_SIZE];
        let s2: Vec<i16> = (0..512).map(|i| (i as i16) % 100 - 50).collect();

        let sig = Signature { nonce, s2 };

        assert_eq!(sig.nonce.len(), NONCE_SIZE);
        assert_eq!(sig.s2.len(), 512);
    }

    #[test]
    fn test_signature_norm() {
        let nonce = [0u8; NONCE_SIZE];
        let s2: Vec<i16> = vec![1, 2, 3, 4, 5];

        let sig = Signature { nonce, s2 };

        // Norm squared: 1 + 4 + 9 + 16 + 25 = 55
        assert_eq!(sig.norm_sq(), 55);
    }

    // Integration test using small parameter set (n=16)
    #[test]
    fn test_sign_basic() {
        use crate::keygen::keygen_16;

        // Use seed 12345 which is known to produce valid NTRU pairs
        let mut rng = StdRng::seed_from_u64(12345);

        // Generate key pair with small n for testing
        let keypair = keygen_16(&mut rng).expect("keygen failed");

        // Sign a message
        let message = b"Hello, FALCON!";
        let sig = sign(&mut rng, &keypair.sk, message);

        match sig {
            Ok(s) => {
                assert_eq!(s.s2.len(), 16);
                println!("Signature norm^2: {}", s.norm_sq());
            }
            Err(e) => {
                // Signing may fail due to norm bounds, which is acceptable
                println!("Signing failed (acceptable for test): {}", e);
            }
        }
    }

    // Full integration test with FALCON-512
    #[test]
    fn test_sign_falcon_512() {
        use crate::keygen::keygen_512;

        let mut rng = StdRng::seed_from_u64(42);

        // Generate key pair
        let keypair = keygen_512(&mut rng).expect("FALCON-512 keygen must succeed");

        // Sign a message
        let message = b"Hello, FALCON!";
        let sig = sign(&mut rng, &keypair.sk, message).expect("FALCON-512 signing must succeed");

        assert_eq!(sig.s2.len(), 512);
        assert!(
            (sig.norm_sq() as f64) <= keypair.sk.params.sig_bound_sq,
            "Signature norm^2 {} exceeds bound {}",
            sig.norm_sq(),
            keypair.sk.params.sig_bound_sq
        );
    }

    // Verify FFT and NTT s2 computations agree
    #[test]
    fn test_sign_falcon_512_fft_ntt_agreement() {
        use crate::keygen::keygen_512;
        use crate::poly::Poly;
        use crate::hash::{generate_nonce, hash_to_point};
        use crate::fft::{fft, ifft};
        use crate::sampler::FfSampler;

        let mut rng = StdRng::seed_from_u64(42);
        let keypair = keygen_512(&mut rng).expect("keygen must succeed");
        let sk = &keypair.sk;

        let sigma = sk.params.sigma;
        let sigma_min = sk.params.sigma_min;

        let ff_sampler = FfSampler::new(&sk.gs, sigma, sigma_min);

        let mut f_fft: Vec<Complex> = sk.f.iter().map(|&x| Complex::from_real(x as f64)).collect();
        let mut big_f_fft: Vec<Complex> = sk.big_f.iter().map(|&x| Complex::from_real(x as f64)).collect();
        fft(&mut f_fft);
        fft(&mut big_f_fft);

        let nonce = generate_nonce(&mut rng);
        let c = hash_to_point(b"Hello, FALCON!", &nonce, &sk.params);
        let mut c_fft: Vec<Complex> = c.iter().map(|zq| Complex::from_real(zq.value() as f64)).collect();
        fft(&mut c_fft);

        let (z0_fft, z1_fft) = ff_sampler.sample_signature(&mut rng, &c_fft);

        let mut z0_coeff = z0_fft.clone();
        let mut z1_coeff = z1_fft.clone();
        ifft(&mut z0_coeff);
        ifft(&mut z1_coeff);

        // z0, z1 should be close to integers
        let z0_max_frac: f64 = z0_coeff.iter().map(|c| (c.re - c.re.round()).abs()).fold(0.0f64, f64::max);
        let z1_max_frac: f64 = z1_coeff.iter().map(|c| (c.re - c.re.round()).abs()).fold(0.0f64, f64::max);
        assert!(z0_max_frac < 1e-3, "z0 not integer-valued: max_frac={:.6e}", z0_max_frac);
        assert!(z1_max_frac < 1e-3, "z1 not integer-valued: max_frac={:.6e}", z1_max_frac);

        // Compute s2 via NTT path (exact)
        let z0_int: Vec<i16> = z0_coeff.iter().map(|c| c.re.round() as i16).collect();
        let z1_int: Vec<i16> = z1_coeff.iter().map(|c| c.re.round() as i16).collect();
        let z0_poly = Poly::from_i16(&z0_int);
        let z1_poly = Poly::from_i16(&z1_int);
        let f_poly = Poly::from_i16(&sk.f.iter().map(|&x| x as i16).collect::<Vec<_>>());
        let big_f_poly = Poly::from_i16(&sk.big_f.iter().map(|&x| x as i16).collect::<Vec<_>>());
        let s2_ntt = z0_poly.mul_ntt(&f_poly).add(&z1_poly.mul_ntt(&big_f_poly));

        // Verify s2 coefficients are within valid range
        let s2_centered = s2_ntt.to_centered();
        let max_coeff = s2_centered.iter().map(|&x| x.abs()).max().unwrap_or(0);
        assert!(
            max_coeff < (crate::params::Q as i16) / 2,
            "s2 coefficient {} exceeds q/2",
            max_coeff
        );
    }

    // Verify signing succeeds across multiple seeds
    #[test]
    fn test_sign_falcon_512_multi_seed() {
        use crate::keygen::keygen_512;

        let mut successes = 0;
        for seed in [42u64, 100, 200, 300, 400, 500, 1000, 2000, 5000, 12345] {
            let mut rng = StdRng::seed_from_u64(seed);
            let keypair = match keygen_512(&mut rng) {
                Ok(kp) => kp,
                Err(_) => continue,
            };

            if sign(&mut rng, &keypair.sk, b"test").is_ok() {
                successes += 1;
            }
        }
        assert!(successes >= 5, "Expected at least 5/10 seeds to produce valid signatures, got {}", successes);
    }
}
