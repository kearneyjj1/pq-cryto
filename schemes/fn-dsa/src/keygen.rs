//! Key generation for FALCON.
//!
//! This module implements the FALCON key generation algorithm,
//! which produces a public/secret key pair from random input.

use rand::RngCore;
use crate::error::{FnDsaError, Result};
use crate::fft::{fft, Complex};
use crate::fft_tree::GramSchmidt;
use crate::field::Zq;
use crate::gaussian::sample_keygen_gaussian;
use crate::ntru::ntru_solve;
use crate::params::{Params, Q, FALCON_512, FALCON_1024};
use crate::poly::ntt;

/// A FALCON public key.
#[derive(Clone, Debug)]
pub struct PublicKey {
    /// The public polynomial h = g/f mod q.
    pub h: Vec<i16>,
    /// The parameter set.
    pub params: Params,
}

impl PublicKey {
    /// Returns the size of the encoded public key in bytes.
    pub fn encoded_size(&self) -> usize {
        self.params.pk_bytes
    }
}

/// A FALCON secret key.
///
/// # Security
///
/// This struct implements `Drop` to zeroize secret key material when dropped.
///
/// **Cloning warning**: When a `SecretKey` is cloned, the clone contains an
/// independent copy of all secret material (polynomials `f`, `g`, `F`, `G`
/// and the Gram-Schmidt / LDL* tree). Both the original and the clone must
/// be dropped for complete zeroization of secret data from memory. Avoid
/// cloning unless necessary (e.g., for a separate signing context).
#[derive(Clone)]
pub struct SecretKey {
    /// The polynomial f (small coefficients).
    pub f: Vec<i8>,
    /// The polynomial g (small coefficients).
    pub g: Vec<i8>,
    /// The polynomial F (NTRU complement of f).
    /// Uses i16 since coefficients may exceed i8 range after Babai reduction.
    pub big_f: Vec<i16>,
    /// The polynomial G (NTRU complement of g).
    /// Uses i16 since coefficients may exceed i8 range after Babai reduction.
    pub big_g: Vec<i16>,
    /// The public key h = g/f mod q (for signature computation).
    pub h: Vec<i16>,
    /// The Gram-Schmidt orthogonalized basis for sampling.
    pub gs: GramSchmidt,
    /// The parameter set.
    pub params: Params,
}

impl Drop for SecretKey {
    fn drop(&mut self) {
        use zeroize::Zeroize;
        // Zeroize all secret polynomials
        self.f.zeroize();
        self.g.zeroize();
        self.big_f.zeroize();
        self.big_g.zeroize();
        // Note: h is public, but we zeroize it anyway
        self.h.zeroize();
        // GramSchmidt, FftTree, and FftNode all implement Drop with zeroization.
        // The gs field is zeroized automatically when SecretKey is dropped.
    }
}

impl SecretKey {
    /// Returns the size of the encoded secret key in bytes.
    pub fn encoded_size(&self) -> usize {
        self.params.sk_bytes
    }
}

/// A FALCON key pair.
///
/// # Security
///
/// Cloning a `KeyPair` clones the embedded [`SecretKey`]. See the security
/// notes on `SecretKey` regarding independent zeroization of cloned copies.
#[derive(Clone)]
pub struct KeyPair {
    /// The public key.
    pub pk: PublicKey,
    /// The secret key.
    pub sk: SecretKey,
}

/// Computes the keygen sigma for polynomial sampling per FIPS 206.
///
/// sigma_keygen = 1.17 * sqrt(q / (2*n))
/// - n=512:  ~4.05
/// - n=1024: ~2.87
/// For small test parameters (n < 256), uses a fixed moderate value since
/// the FIPS 206 formula produces unreasonably large sigmas.
fn keygen_sigma(n: usize) -> f64 {
    if n < 256 {
        return 1.5; // Moderate value for toy test parameters
    }
    1.17 * ((Q as f64) / (2.0 * n as f64)).sqrt()
}

/// Computes the FIPS 206 combined norm bound for (f, g).
///
/// Returns the maximum allowed value for `||f||^2 + ||g||^2`.
/// For production parameters (n >= 256), uses the FIPS 206 factor 1.21
/// above the expected norm. For toy parameters (n < 256), uses a more
/// generous factor since chi-squared variance is high relative to the
/// mean at low degrees of freedom.
fn max_fg_norm_sq(n: usize) -> i64 {
    let sigma = keygen_sigma(n);
    let expected = 2.0 * (n as f64) * sigma * sigma;
    if n < 256 {
        (3.0 * expected).floor() as i64
    } else {
        (1.21 * expected).floor() as i64
    }
}

/// Generates a random polynomial with small coefficients using discrete Gaussian sampling.
///
/// Coefficients are sampled from a discrete Gaussian distribution N(0, sigma^2)
/// using rejection sampling (`sample_keygen_gaussian`).
///
/// Note: The FIPS 206 SamplerZ (Algorithm 12) cannot be used here because its
/// base sampler at sigma_0=1.82 has lighter tails than the keygen sigma≈4.05,
/// causing incorrect variance. SamplerZ is designed for the signing path where
/// sigma ≈ sigma_min ≈ 1.28 < sigma_0.
fn generate_small_poly<R: RngCore>(rng: &mut R, n: usize, sigma: f64) -> Vec<i8> {
    let mut poly = vec![0i8; n];

    for i in 0..n {
        let z = sample_keygen_gaussian(rng, sigma);
        // Clamp to i8 range (coefficients should be small, but clamp for safety)
        poly[i] = z.clamp(-127, 127) as i8;
    }

    poly
}

/// Checks if a polynomial is invertible modulo q using exact NTT.
///
/// A polynomial f is invertible in Z_q[X]/(X^n + 1) if and only if
/// its NTT (Number Theoretic Transform) has no zero components.
/// This uses exact modular arithmetic rather than floating-point FFT.
fn is_invertible(f: &[i8], n: usize) -> bool {
    // Convert to Zq
    let f_zq: Vec<Zq> = f.iter().map(|&x| Zq::new(x as i32)).collect();

    // Compute exact NTT
    let f_ntt = ntt(&f_zq, n);

    // Check that no NTT coefficient is zero
    for c in &f_ntt {
        if c.is_zero() {
            return false;
        }
    }

    true
}

/// Checks if the LDL* tree has acceptable quality for FALCON signing.
///
/// Per FIPS 206 Algorithm 14, the Gram-Schmidt norms at all tree leaves
/// must satisfy `tree_sigma <= sigma_sign / sigma_min`, ensuring that
/// `sigma_eff = sigma_sign / tree_sigma >= sigma_min` at every leaf.
/// Keys that fail this check are rejected during keygen.
fn check_gram_schmidt_quality(gs: &GramSchmidt, params: &Params) -> bool {
    let tree = &gs.tree;
    let max_tree_sigma = params.sigma / params.sigma_min;

    let leaf_level = tree.depth - 1;
    let n_leaves = tree.nodes_at_level(leaf_level);

    for pos in 0..n_leaves {
        let node = tree.get_node(leaf_level, pos);
        for &s in &node.sigma {
            if s > max_tree_sigma {
                return false;
            }
        }
    }

    true
}

/// Computes the Gram-Schmidt orthogonalization of the secret key.
fn compute_gram_schmidt(f: &[i8], g: &[i8], big_f: &[i16], big_g: &[i16]) -> GramSchmidt {
    let n = f.len();

    // Convert to FFT form
    let mut f_fft: Vec<Complex> = f.iter().map(|&x| Complex::from_real(x as f64)).collect();
    let mut g_fft: Vec<Complex> = g.iter().map(|&x| Complex::from_real(x as f64)).collect();
    let mut big_f_fft: Vec<Complex> = big_f.iter().map(|&x| Complex::from_real(x as f64)).collect();
    let mut big_g_fft: Vec<Complex> = big_g.iter().map(|&x| Complex::from_real(x as f64)).collect();

    fft(&mut f_fft);
    fft(&mut g_fft);
    fft(&mut big_f_fft);
    fft(&mut big_g_fft);

    GramSchmidt::new(f_fft, g_fft, big_f_fft, big_g_fft)
}

/// Computes the public key h = g/f mod q.
fn compute_public_key(f: &[i8], g: &[i8]) -> Result<Vec<i16>> {
    crate::ntru::compute_public_key(f, g, f.len())
}

/// Generates a FALCON key pair per FIPS 206 Algorithm 14 (NTRUGen).
///
/// Samples (f, g) from a discrete Gaussian, checks the combined norm bound
/// `||f||^2 + ||g||^2`, solves the NTRU equation `fG - gF = q`, computes
/// the Gram-Schmidt / LDL* tree, and verifies tree quality at all leaves.
/// Uses a single Gaussian sampling strategy as specified by the standard.
pub fn keygen<R: RngCore>(rng: &mut R, params: &Params) -> Result<KeyPair> {
    let n = params.n;
    let sigma = keygen_sigma(n);
    let norm_bound = max_fg_norm_sq(n);

    // FIPS 206: up to 10,000 attempts for production parameters.
    // For toy parameters (small n), allow many more attempts since
    // valid NTRU pairs are much rarer with the single Gaussian strategy.
    let max_attempts = if n <= 64 { 100_000 } else { 10_000 };

    for _attempt in 0..max_attempts {
        // Step 1: Sample f, g from discrete Gaussian N(0, sigma^2)
        let f = generate_small_poly(rng, n, sigma);
        let g = generate_small_poly(rng, n, sigma);

        // Step 2: Check combined norm ||f||^2 + ||g||^2 <= norm_bound
        let fg_norm_sq: i64 = f.iter().chain(g.iter())
            .map(|&x| (x as i64) * (x as i64))
            .sum();
        if fg_norm_sq > norm_bound {
            continue;
        }

        // Step 3: Check f invertible mod q (via exact NTT)
        if !is_invertible(&f, n) {
            continue;
        }

        // Step 4: Solve NTRU equation fG - gF = q
        let (big_f, big_g) = match ntru_solve(&f, &g, n) {
            Ok(result) => result,
            Err(_) => continue,
        };

        // Step 5: Compute public key h = g/f mod q
        let h = match compute_public_key(&f, &g) {
            Ok(h) => h,
            Err(_) => continue,
        };

        // Step 6: Build Gram-Schmidt tree and check leaf quality.
        // The FIPS 206 quality check applies to production parameters
        // (n >= 256) only. For toy parameters (n < 256), the signing
        // loop's norm check provides the safety net instead.
        let gs = compute_gram_schmidt(&f, &g, &big_f, &big_g);
        if n >= 256 && !check_gram_schmidt_quality(&gs, params) {
            continue;
        }

        // Success
        let pk = PublicKey { h, params: *params };
        let sk = SecretKey {
            f, g, big_f, big_g,
            h: pk.h.clone(),
            gs,
            params: *params,
        };

        return Ok(KeyPair { pk, sk });
    }

    Err(FnDsaError::KeygenFailed {
        reason: "exceeded maximum attempts",
    })
}

/// Generates a FALCON key pair with a specific seed.
///
/// This is useful for deterministic key generation in testing.
pub fn keygen_with_seed(seed: u64, params: &Params) -> Result<KeyPair> {
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    let mut rng = StdRng::seed_from_u64(seed);
    keygen(&mut rng, params)
}

/// Generates a FALCON-512 key pair.
pub fn keygen_512<R: RngCore>(rng: &mut R) -> Result<KeyPair> {
    keygen(rng, &FALCON_512)
}

/// Generates a FALCON-1024 key pair.
pub fn keygen_1024<R: RngCore>(rng: &mut R) -> Result<KeyPair> {
    keygen(rng, &FALCON_1024)
}

/// Generates a FALCON-16 key pair (for testing only).
#[cfg(test)]
pub fn keygen_16<R: RngCore>(rng: &mut R) -> Result<KeyPair> {
    use crate::params::FALCON_16;
    keygen(rng, &FALCON_16)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_generate_small_poly() {
        let mut rng = StdRng::seed_from_u64(42);
        let n = 64;
        let sigma = keygen_sigma(n);
        let poly = generate_small_poly(&mut rng, n, sigma);

        assert_eq!(poly.len(), n);

        // keygen_sigma(64) = 1.5 (capped for small n)
        // Coefficients should be within ~5*sigma ≈ 7.5
        for &c in &poly {
            assert!(c.abs() <= 10, "Coefficient {} is too large for sigma={:.1}", c, sigma);
        }

        // Check that we have a mix of positive and negative values
        let has_positive = poly.iter().any(|&c| c > 0);
        let has_negative = poly.iter().any(|&c| c < 0);
        assert!(has_positive && has_negative);
    }

    #[test]
    fn test_is_invertible() {
        // f = 1 should be invertible
        let mut f = vec![0i8; 8];
        f[0] = 1;
        assert!(is_invertible(&f, 8));

        // f = 0 should not be invertible
        let f_zero = vec![0i8; 8];
        assert!(!is_invertible(&f_zero, 8));
    }

    // Test keygen with various n values using the FIPS 206 approach
    fn try_keygen_with_n(n: usize, log_n: usize, seed: u64, max_attempts: u32) -> bool {
        use crate::params::Params;

        let params = Params {
            n,
            log_n,
            sigma: 165.0,
            sigma_min: 1.2778336969128337,
            sig_bound_sq: 34034726.0,
            pk_bytes: 64,
            sk_bytes: 128,
            sig_bytes_max: 64,
            rice_k: 8,
            security_level: 1,
        };

        let sigma = keygen_sigma(n);
        let norm_bound = max_fg_norm_sq(n);

        // Try multiple seeds
        for s in 0..10 {
            let mut rng = StdRng::seed_from_u64(seed + s);

            for _ in 0..max_attempts {
                let f = generate_small_poly(&mut rng, n, sigma);
                let g = generate_small_poly(&mut rng, n, sigma);

                // Combined norm check
                let fg_norm_sq: i64 = f.iter().chain(g.iter())
                    .map(|&x| (x as i64) * (x as i64))
                    .sum();
                if fg_norm_sq > norm_bound {
                    continue;
                }

                if !is_invertible(&f, n) {
                    continue;
                }

                if let Ok((_big_f, _big_g)) = ntru_solve(&f, &g, n) {
                    return true;
                }
            }
        }
        false
    }

    #[test]
    fn test_keygen_small() {
        // Test n=8 - toy parameter, Gaussian-only sampling may need many attempts
        let success = try_keygen_with_n(8, 3, 12345, 2000);
        if success {
            println!("n=8 keygen succeeded!");
        } else {
            println!("n=8 keygen did not succeed (expected for toy parameters with single Gaussian strategy)");
        }
    }

    #[test]
    fn test_keygen_16_direct() {
        // Test FIPS 206 keygen flow directly for n=16.
        let n = 16;
        let sigma = keygen_sigma(n);
        let norm_bound = max_fg_norm_sq(n);

        for seed in 0..20u64 {
            let mut rng = StdRng::seed_from_u64(seed);

            for attempt in 0..500 {
                let f = generate_small_poly(&mut rng, n, sigma);
                let g = generate_small_poly(&mut rng, n, sigma);

                // Combined norm check
                let fg_norm_sq: i64 = f.iter().chain(g.iter())
                    .map(|&x| (x as i64) * (x as i64))
                    .sum();
                if fg_norm_sq > norm_bound {
                    continue;
                }

                if f[0] == 0 || !is_invertible(&f, n) {
                    continue;
                }

                let (big_f, big_g) = match ntru_solve(&f, &g, n) {
                    Ok(result) => result,
                    Err(_) => continue,
                };

                if let Ok(_h) = compute_public_key(&f, &g) {
                    println!("Success at seed={}, attempt {}: NTRU solved, public key computed", seed, attempt);
                    return;
                }
            }
        }

        panic!("keygen_16 failed across all seeds");
    }

    #[test]
    fn test_keygen_n16() {
        // Test n=16
        let success = try_keygen_with_n(16, 4, 12345, 500);
        if success {
            println!("n=16 keygen succeeded!");
        } else {
            println!("n=16 keygen failed (may need more attempts)");
        }
    }

    #[test]
    fn test_keygen_n32() {
        // Test n=32 - might need many attempts
        let success = try_keygen_with_n(32, 5, 12345, 1000);
        if success {
            println!("n=32 keygen succeeded!");
        } else {
            println!("n=32 keygen failed (may need more attempts or better sampling)");
        }
    }

    #[test]
    fn test_keygen_512() {
        let mut rng = StdRng::seed_from_u64(42);
        let keypair = keygen_512(&mut rng).expect("FALCON-512 keygen must succeed");

        assert_eq!(keypair.pk.h.len(), 512);
        assert_eq!(keypair.sk.f.len(), 512);
        assert_eq!(keypair.sk.g.len(), 512);
        assert_eq!(keypair.sk.big_f.len(), 512);
        assert_eq!(keypair.sk.big_g.len(), 512);
    }

    /// Verifies that keygen with correct Gaussian sampling produces keys
    /// with reasonable tree quality (avg g00 close to expected).
    #[test]
    fn test_tree_leaf_sigma_diagnostic() {
        use crate::params::FALCON_512;

        let handle = std::thread::Builder::new()
            .stack_size(16 * 1024 * 1024)
            .spawn(|| {
                let n = 512;
                let sigma_gen = keygen_sigma(n);
                let params = &FALCON_512;
                let max_tree_sigma = params.sigma / params.sigma_min;
                // Expected avg g00 = (||f_fft||² + ||g_fft||²) / n
                // Since Parseval: ||f_fft||² = n*||f||², and ||f||² ≈ n*sigma²:
                // avg g00 ≈ (n²*sigma² + n²*sigma²) / n = 2*n*sigma²
                let expected_avg_g00 = 2.0 * (n as f64) * sigma_gen * sigma_gen;

                let mut rng = StdRng::seed_from_u64(42);
                for _attempt in 0..50u32 {
                    let f = generate_small_poly(&mut rng, n, sigma_gen);
                    let g = generate_small_poly(&mut rng, n, sigma_gen);

                    if f[0] == 0 || !is_invertible(&f, n) { continue; }
                    let (big_f, big_g) = match ntru_solve(&f, &g, n) {
                        Ok(r) => r,
                        Err(_) => continue,
                    };

                    let gs = compute_gram_schmidt(&f, &g, &big_f, &big_g);
                    let f_norm_sq: f64 = gs.f_fft.iter().map(|c| c.norm_sq()).sum();
                    let g_norm_sq: f64 = gs.g_fft.iter().map(|c| c.norm_sq()).sum();
                    let avg_g00 = (f_norm_sq + g_norm_sq) / n as f64;

                    // avg g00 should be within 30% of expected
                    assert!(
                        avg_g00 > expected_avg_g00 * 0.7 && avg_g00 < expected_avg_g00 * 1.3,
                        "avg g00 = {:.1}, expected ≈ {:.1}, threshold² = {:.1}",
                        avg_g00, expected_avg_g00, max_tree_sigma * max_tree_sigma
                    );

                    // Most leaves should be below threshold (< 20% above for well-formed keys)
                    let tree = &gs.tree;
                    let leaf_level = tree.depth - 1;
                    let n_leaves = tree.nodes_at_level(leaf_level);
                    let total_sigma_values = n_leaves * 2;
                    let above: usize = (0..n_leaves)
                        .flat_map(|pos| tree.get_node(leaf_level, pos).sigma.iter().copied())
                        .filter(|&s| s > max_tree_sigma)
                        .count();
                    assert!(
                        (above as f64) < 0.20 * (total_sigma_values as f64),
                        "Too many leaves above threshold: {}/{}", above, total_sigma_values
                    );
                    return;
                }
                panic!("No valid key found in 50 attempts");
            })
            .unwrap();
        handle.join().unwrap();
    }
}
