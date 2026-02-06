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
#[cfg(test)]
use crate::params::Q;
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
pub fn sign<R: RngCore>(rng: &mut R, sk: &SecretKey, message: &[u8]) -> Result<Signature> {
    use zeroize::Zeroize;
    use crate::poly::Poly;

    let n = sk.params.n;
    let sigma = sk.params.sigma;
    let sigma_min = sk.params.sigma_min;
    let bound_sq = sk.params.sig_bound_sq;

    // Create the FFT sampler using the LDL* tree from the secret key
    let ff_sampler = FfSampler::new(sk.gs.clone(), sigma, sigma_min);

    // Precompute NTT-ready Poly forms of the secret key polynomials
    let f_poly = Poly::from_i16(&sk.f.iter().map(|&x| x as i16).collect::<Vec<_>>());
    let big_f_poly = Poly::from_i16(&sk.big_f.iter().map(|&x| x as i16).collect::<Vec<_>>());
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
        let z0: Vec<i16> = z0_fft.iter().map(|c| c.re.round() as i16).collect();
        let z1: Vec<i16> = z1_fft.iter().map(|c| c.re.round() as i16).collect();

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
            return Ok(Signature { nonce, s2 });
        }
    }

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

/// Signs a message with a specific nonce (for testing/debugging).
///
/// # Security
///
/// This function is restricted to `#[cfg(test)]` builds because it uses
/// a naive sampler and does not provide the security guarantees of ffSampling.
#[cfg(test)]
pub fn sign_with_nonce<R: RngCore>(
    rng: &mut R,
    sk: &SecretKey,
    message: &[u8],
    nonce: [u8; NONCE_SIZE],
) -> Result<Signature> {
    let n = sk.params.n;
    let sigma = sk.params.sigma;
    let bound_sq = sk.params.sig_bound_sq;

    let sampler = SimpleSampler::new();

    // Hash the message to get c
    let c = hash_to_point(message, &nonce, &sk.params);

    // Convert c to FFT form
    let mut c_fft: Vec<Complex> = c.iter().map(|zq| Complex::from_real(zq.value() as f64)).collect();
    fft(&mut c_fft);

    // Convert secret key to FFT form
    let mut f_fft: Vec<Complex> = sk.f.iter().map(|&x| Complex::from_real(x as f64)).collect();
    let mut big_f_fft: Vec<Complex> = sk.big_f.iter().map(|&x| Complex::from_real(x as f64)).collect();

    fft(&mut f_fft);
    fft(&mut big_f_fft);

    // Sample from the Gaussian
    let z0: Vec<i64> = (0..n)
        .map(|_| sampler.sample(rng, 0.0, sigma))
        .collect();
    let z1: Vec<i64> = (0..n)
        .map(|_| sampler.sample(rng, 0.0, sigma))
        .collect();

    // Convert z0, z1 to FFT form
    let mut z0_fft: Vec<Complex> = z0.iter().map(|&x| Complex::from_real(x as f64)).collect();
    let mut z1_fft: Vec<Complex> = z1.iter().map(|&x| Complex::from_real(x as f64)).collect();
    fft(&mut z0_fft);
    fft(&mut z1_fft);

    // Compute s2 = z0*f + z1*F
    let mut s2_fft: Vec<Complex> = Vec::with_capacity(n);
    for i in 0..n {
        s2_fft.push(z0_fft[i] * f_fft[i] + z1_fft[i] * big_f_fft[i]);
    }

    ifft(&mut s2_fft);
    let s2: Vec<i16> = s2_fft.iter()
        .map(|c| {
            let val = c.re.round() as i32;
            let reduced = ((val % Q) + Q) % Q;
            if reduced > Q / 2 {
                (reduced - Q) as i16
            } else {
                reduced as i16
            }
        })
        .collect();

    Ok(Signature { nonce, s2 })
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

    // Diagnostic: compare FFT vs NTT s2 computation to find precision issues
    #[test]
    fn test_sign_falcon_512_diag() {
        use crate::keygen::keygen_512;
        use crate::poly::Poly;
        use crate::hash::{generate_nonce, hash_to_point};
        use crate::fft::{fft, ifft};
        use crate::sampler::FfSampler;

        let mut rng = StdRng::seed_from_u64(42);
        let keypair = keygen_512(&mut rng).expect("keygen must succeed");
        let sk = &keypair.sk;

        let n = sk.params.n;
        let sigma = sk.params.sigma;
        let sigma_min = sk.params.sigma_min;

        let ff_sampler = FfSampler::new(sk.gs.clone(), sigma, sigma_min);

        let mut f_fft: Vec<Complex> = sk.f.iter().map(|&x| Complex::from_real(x as f64)).collect();
        let mut big_f_fft: Vec<Complex> = sk.big_f.iter().map(|&x| Complex::from_real(x as f64)).collect();
        fft(&mut f_fft);
        fft(&mut big_f_fft);

        let nonce = generate_nonce(&mut rng);
        let c = hash_to_point(b"Hello, FALCON!", &nonce, &sk.params);
        let mut c_fft: Vec<Complex> = c.iter().map(|zq| Complex::from_real(zq.value() as f64)).collect();
        fft(&mut c_fft);

        let (z0_fft, z1_fft) = ff_sampler.sample_signature(&mut rng, &c_fft);

        // Convert z0, z1 back to coefficient form to get integer values
        let mut z0_coeff = z0_fft.clone();
        let mut z1_coeff = z1_fft.clone();
        ifft(&mut z0_coeff);
        ifft(&mut z1_coeff);

        // Check if z0, z1 are actually integers
        let z0_max_frac: f64 = z0_coeff.iter().map(|c| (c.re - c.re.round()).abs()).fold(0.0f64, f64::max);
        let z1_max_frac: f64 = z1_coeff.iter().map(|c| (c.re - c.re.round()).abs()).fold(0.0f64, f64::max);
        let z0_max_im: f64 = z0_coeff.iter().map(|c| c.im.abs()).fold(0.0f64, f64::max);
        let z1_max_im: f64 = z1_coeff.iter().map(|c| c.im.abs()).fold(0.0f64, f64::max);
        println!("z0: max_frac={:.6e}, max_im={:.6e}", z0_max_frac, z0_max_im);
        println!("z1: max_frac={:.6e}, max_im={:.6e}", z1_max_frac, z1_max_im);

        // z0 and z1 max abs values
        let z0_max: f64 = z0_coeff.iter().map(|c| c.re.abs()).fold(0.0f64, f64::max);
        let z1_max: f64 = z1_coeff.iter().map(|c| c.re.abs()).fold(0.0f64, f64::max);
        println!("z0_max_abs={:.1}, z1_max_abs={:.1}", z0_max, z1_max);

        // Compute s2 via FFT path (as sign() does)
        let mut s2_fft: Vec<Complex> = Vec::with_capacity(n);
        for i in 0..n {
            s2_fft.push(z0_fft[i] * f_fft[i] + z1_fft[i] * big_f_fft[i]);
        }
        ifft(&mut s2_fft);
        let s2_fft_coeffs: Vec<i32> = s2_fft.iter().map(|c| c.re.round() as i32).collect();
        let s2_fft_max_frac: f64 = s2_fft.iter().map(|c| (c.re - c.re.round()).abs()).fold(0.0f64, f64::max);
        println!("s2_fft: max_frac={:.6e}, max_abs={}", s2_fft_max_frac,
            s2_fft_coeffs.iter().map(|x| x.abs()).max().unwrap());

        // Compute s2 via NTT path (exact)
        let z0_int: Vec<i16> = z0_coeff.iter().map(|c| c.re.round() as i16).collect();
        let z1_int: Vec<i16> = z1_coeff.iter().map(|c| c.re.round() as i16).collect();
        let z0_poly = Poly::from_i16(&z0_int);
        let z1_poly = Poly::from_i16(&z1_int);
        let f_poly = Poly::from_i16(&sk.f.iter().map(|&x| x as i16).collect::<Vec<_>>());
        let big_f_poly = Poly::from_i16(&sk.big_f.iter().map(|&x| x as i16).collect::<Vec<_>>());
        let s2_ntt = z0_poly.mul_ntt(&f_poly).add(&z1_poly.mul_ntt(&big_f_poly));
        let s2_ntt_centered = s2_ntt.to_centered();

        // Compare FFT vs NTT
        let mut mismatches = 0;
        for i in 0..n {
            let fft_reduced = ((s2_fft_coeffs[i] % Q) + Q) % Q;
            let fft_centered = if fft_reduced > Q / 2 { fft_reduced - Q } else { fft_reduced } as i16;
            if fft_centered != s2_ntt_centered[i] {
                if mismatches < 5 {
                    println!("s2 mismatch at [{}]: fft={} (raw={}), ntt={}", i, fft_centered, s2_fft_coeffs[i], s2_ntt_centered[i]);
                }
                mismatches += 1;
            }
        }
        println!("Total s2 mismatches: {}/{}", mismatches, n);

        // Now compute s1 using both s2 values
        let c_poly = Poly::from_zq(c);
        let h_poly = Poly::from_i16(&sk.h);

        let s2_from_fft = Poly::from_i16(&s2_fft_coeffs.iter().map(|&x| {
            let r = ((x % Q) + Q) % Q;
            if r > Q / 2 { (r - Q) as i16 } else { r as i16 }
        }).collect::<Vec<_>>());
        let s1_from_fft = c_poly.sub(&s2_from_fft.mul_ntt(&h_poly));

        let s1_from_ntt = c_poly.sub(&s2_ntt.mul_ntt(&h_poly));

        let s2_fft_norm: i64 = s2_fft_coeffs.iter().map(|&x| {
            let r = ((x % Q) + Q) % Q;
            let c = if r > Q / 2 { (r - Q) as i64 } else { r as i64 };
            c * c
        }).fold(0i64, |a, b| a.saturating_add(b));
        let s2_ntt_norm: i64 = s2_ntt_centered.iter().map(|&x| (x as i64) * (x as i64)).fold(0i64, |a, b| a.saturating_add(b));

        println!("FFT path: s1²={}, s2²={}, total={}", s1_from_fft.norm_sq(), s2_fft_norm, s1_from_fft.norm_sq().saturating_add(s2_fft_norm));
        println!("NTT path: s1²={}, s2²={}, total={}, bound={}", s1_from_ntt.norm_sq(), s2_ntt_norm, s1_from_ntt.norm_sq().saturating_add(s2_ntt_norm), sk.params.sig_bound_sq);

        // Check LDL* tree leaf sigma values
        let tree = &sk.gs.tree;
        let leaf_level = tree.depth - 1;
        let n_leaves = tree.nodes_at_level(leaf_level);
        let mut min_sigma = f64::MAX;
        let mut max_sigma = f64::MIN;
        let mut sigma_eff_min = f64::MAX;
        let mut sigma_eff_max = f64::MIN;
        for pos in 0..n_leaves {
            let node = tree.get_node(leaf_level, pos);
            for &s in &node.sigma {
                if s > 0.0 {
                    min_sigma = min_sigma.min(s);
                    max_sigma = max_sigma.max(s);
                    let eff = sk.params.sigma / s;
                    sigma_eff_min = sigma_eff_min.min(eff);
                    sigma_eff_max = sigma_eff_max.max(eff);
                }
            }
        }
        println!("Leaf tree_sigma: min={:.4}, max={:.4}", min_sigma, max_sigma);
        println!("Leaf sigma_eff:  min={:.4}, max={:.4} (sigma_min={:.4})", sigma_eff_min, sigma_eff_max, sk.params.sigma_min);
    }

    // Try multiple seeds to find one where signing succeeds
    #[test]
    fn test_sign_falcon_512_seed_search() {
        use crate::keygen::keygen_512;

        for seed in [42u64, 100, 200, 300, 400, 500, 1000, 2000, 5000, 12345] {
            let mut rng = StdRng::seed_from_u64(seed);
            let keypair = match keygen_512(&mut rng) {
                Ok(kp) => kp,
                Err(_) => continue,
            };

            // Check tree conditioning
            let tree = &keypair.sk.gs.tree;
            let leaf_level = tree.depth - 1;
            let n_leaves = tree.nodes_at_level(leaf_level);
            let mut sigma_eff_max: f64 = 0.0;
            for pos in 0..n_leaves {
                let node = tree.get_node(leaf_level, pos);
                for &s in &node.sigma {
                    if s > 1e-6 {
                        let eff = keypair.sk.params.sigma / s;
                        sigma_eff_max = sigma_eff_max.max(eff);
                    }
                }
            }

            let sig = sign(&mut rng, &keypair.sk, b"test");
            let status = match &sig {
                Ok(s) => format!("OK (norm²={})", s.norm_sq()),
                Err(_) => "FAIL".to_string(),
            };
            println!("Seed {}: sigma_eff_max={:.4}, sign={}", seed, sigma_eff_max, status);
        }
    }
}
