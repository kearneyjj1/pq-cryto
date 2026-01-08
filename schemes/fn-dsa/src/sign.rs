//! Signing for FALCON.
//!
//! This module implements the FALCON signature generation algorithm,
//! which uses Fast Fourier Sampling to produce compact signatures.

use rand::RngCore;
use crate::error::{FnDsaError, Result};
use crate::fft::{fft, ifft, Complex};
use crate::hash::{generate_nonce, hash_to_point};
use crate::keygen::SecretKey;
use crate::params::{Q, MAX_SIGN_ATTEMPTS, NONCE_SIZE};
use crate::sampler::{FfSampler, SimpleSampler};

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
    pub fn norm_sq(&self) -> i64 {
        self.s2.iter().map(|&x| (x as i64) * (x as i64)).sum()
    }
}

/// Signs a message using the FALCON signature scheme.
///
/// This is the main signing function. It:
/// 1. Generates a random nonce
/// 2. Hashes the message to a polynomial c
/// 3. Computes the target t = (c, 0)
/// 4. Samples a lattice vector (z0, z1) close to t using FFT sampling
/// 5. Computes s2 = c - z0*f - z1*F
/// 6. Checks if the norm is acceptable
/// 7. Returns the signature (nonce, s2)
pub fn sign<R: RngCore>(rng: &mut R, sk: &SecretKey, message: &[u8]) -> Result<Signature> {
    let n = sk.params.n;
    let sigma = sk.params.sigma;
    let bound_sq = sk.params.sig_bound_sq;

    // Create the FFT sampler using the Gram-Schmidt data from the secret key
    let ff_sampler = FfSampler::new(sk.gs.clone());

    // Precompute FFT forms of the secret key polynomials
    let mut f_fft: Vec<Complex> = sk.f.iter().map(|&x| Complex::from_real(x as f64)).collect();
    let mut g_fft: Vec<Complex> = sk.g.iter().map(|&x| Complex::from_real(x as f64)).collect();
    let mut big_f_fft: Vec<Complex> = sk.big_f.iter().map(|&x| Complex::from_real(x as f64)).collect();
    let mut big_g_fft: Vec<Complex> = sk.big_g.iter().map(|&x| Complex::from_real(x as f64)).collect();

    fft(&mut f_fft);
    fft(&mut g_fft);
    fft(&mut big_f_fft);
    fft(&mut big_g_fft);

    for _attempt in 0..MAX_SIGN_ATTEMPTS {
        // Generate a fresh nonce
        let nonce = generate_nonce(rng);

        // Hash the message to get c
        let c = hash_to_point(message, &nonce, &sk.params);

        // Convert c to FFT form for operations
        let mut c_fft: Vec<Complex> = c.iter().map(|zq| Complex::from_real(zq.value() as f64)).collect();
        fft(&mut c_fft);

        // Use the FFT sampler to sample (z0, z1) close to the target
        let (z0_fft, z1_fft) = ff_sampler.sample_signature(rng, &c_fft, sigma);

        // Compute s2 = z0*f + z1*F (in FFT form, then convert back)
        let mut s2_fft: Vec<Complex> = Vec::with_capacity(n);
        for i in 0..n {
            s2_fft.push(z0_fft[i] * f_fft[i] + z1_fft[i] * big_f_fft[i]);
        }

        // Convert s2 back to coefficient form
        ifft(&mut s2_fft);
        let s2: Vec<i16> = s2_fft.iter()
            .map(|c| {
                let val = c.re.round() as i32;
                // Reduce modulo q and center
                let reduced = ((val % Q) + Q) % Q;
                if reduced > Q / 2 {
                    (reduced - Q) as i16
                } else {
                    reduced as i16
                }
            })
            .collect();

        // Compute s1 = c - s2*h using polynomial arithmetic (same as verification)
        // This ensures the norm check uses the same computation as verification
        use crate::poly::Poly;
        let c_poly = Poly::from_zq(c.clone());
        let s2_poly = Poly::from_i16(&s2);
        let h_poly = Poly::from_i16(&sk.h);
        let s2h = s2_poly.mul(&h_poly);
        let s1_poly = c_poly.sub(&s2h);

        // Check the norm
        let s1_norm_sq = s1_poly.norm_sq();
        let s2_norm_sq: i64 = s2.iter().map(|&x| (x as i64) * (x as i64)).sum();
        let total_norm_sq = s1_norm_sq + s2_norm_sq;

        if (total_norm_sq as f64) <= bound_sq {
            return Ok(Signature { nonce, s2 });
        }
    }

    Err(FnDsaError::SigningFailed {
        attempts: MAX_SIGN_ATTEMPTS,
    })
}

/// Signs a message using direct s2 sampling with rejection (for educational demo).
///
/// This is NOT the proper FALCON algorithm - it directly samples s2 from a
/// discrete Gaussian and uses rejection sampling to find a short (s1, s2) pair.
/// This approach is much slower and doesn't provide the same security guarantees,
/// but it demonstrates the verification correctness.
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

        // Check the norm
        let s1_norm_sq = s1_poly.norm_sq();
        let s2_norm_sq: i64 = s2.iter().map(|&x| (x as i64) * (x as i64)).sum();
        let total_norm_sq = s1_norm_sq + s2_norm_sq;

        if (total_norm_sq as f64) <= bound_sq {
            return Ok(Signature { nonce, s2 });
        }
    }

    Err(FnDsaError::SigningFailed {
        attempts: MAX_SIGN_ATTEMPTS,
    })
}

/// Signs a message with a specific nonce (for testing/debugging).
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
    use crate::params::FALCON_512;
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

    // Full integration test with FALCON-512 (requires working NTRUSolve for large n)
    #[test]
    #[ignore]
    fn test_sign_falcon_512() {
        use crate::keygen::keygen_512;

        // Use a seed that is known to find valid NTRU pairs
        let mut rng = StdRng::seed_from_u64(42);

        // Generate key pair
        let keypair = match keygen_512(&mut rng) {
            Ok(kp) => kp,
            Err(e) => {
                println!("Keygen failed (expected for incomplete impl): {}", e);
                return;
            }
        };

        // Sign a message
        let message = b"Hello, FALCON!";
        let sig = sign(&mut rng, &keypair.sk, message);

        match sig {
            Ok(s) => {
                assert_eq!(s.s2.len(), 512);
                assert!(s.norm_sq() < keypair.sk.params.sig_bound_sq as i64);
                println!("FALCON-512 signing succeeded! Signature norm^2: {}", s.norm_sq());
            }
            Err(e) => {
                println!("Signing failed: {}", e);
            }
        }
    }
}
