//! Sampling functions for ML-DSA.
//!
//! This module implements the deterministic sampling functions used in ML-DSA,
//! all based on SHAKE128 or SHAKE256 as the underlying extendable-output function.
//!
//! Functions:
//! - ExpandA: Expands seed ρ to public matrix A (FIPS 204 Algorithm 26)
//! - ExpandS: Expands seed ρ' to secret vectors s1, s2 (FIPS 204 Algorithm 27)
//! - ExpandMask: Expands seed to masking polynomial y (FIPS 204 Algorithm 28)
//! - SampleInBall: Samples challenge polynomial c (FIPS 204 Algorithm 29)

use crate::params::{Params, N, Q};
use crate::poly::Poly;
use crate::polyvec::{PolyMatrix, PolyVec};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128, Shake256,
};

/// Rejection bound for uniform sampling (q = 8380417 < 2^23).
/// We sample 3 bytes at a time and reject if >= Q.
const SAMPLE_BOUND: i32 = Q;

/// Expands seed ρ into the public matrix A.
///
/// A is a k×l matrix where each element is a polynomial with coefficients
/// uniformly distributed in [0, q).
///
/// Uses SHAKE128 for each matrix entry with domain separator.
///
/// FIPS 204 Algorithm 26: ExpandA
pub fn expand_a(rho: &[u8; 32], params: &Params) -> PolyMatrix {
    let k = params.k;
    let l = params.l;
    let mut matrix = PolyMatrix::zero(k, l);

    for i in 0..k {
        for j in 0..l {
            // Domain separator: ρ || j || i (little-endian bytes)
            let mut hasher = Shake128::default();
            hasher.update(rho);
            hasher.update(&[j as u8, i as u8]);
            let mut reader = hasher.finalize_xof();

            // Sample 256 coefficients using rejection sampling
            let poly = &mut matrix.rows[i].polys[j];
            sample_uniform_poly(&mut reader, poly);
        }
    }

    matrix
}

/// Samples a polynomial with uniform coefficients in [0, q) using rejection sampling.
fn sample_uniform_poly<R: XofReader>(reader: &mut R, poly: &mut Poly) {
    let mut coeffs_sampled = 0;
    let mut buf = [0u8; 3];

    while coeffs_sampled < N {
        reader.read(&mut buf);

        // Extract 23-bit value from 3 bytes
        let val = (buf[0] as i32) | ((buf[1] as i32) << 8) | (((buf[2] & 0x7F) as i32) << 16);

        if val < SAMPLE_BOUND {
            poly.coeffs[coeffs_sampled] = val;
            coeffs_sampled += 1;
        }
    }
}

/// Expands seed ρ' into secret vectors s1 and s2.
///
/// - s1 is a vector of l polynomials with coefficients in [-η, η]
/// - s2 is a vector of k polynomials with coefficients in [-η, η]
///
/// Uses SHAKE256 for each polynomial with domain separator.
///
/// FIPS 204 Algorithm 27: ExpandS
pub fn expand_s(rho_prime: &[u8; 64], params: &Params) -> (PolyVec, PolyVec) {
    let k = params.k;
    let l = params.l;
    let eta = params.eta;

    let mut s1 = PolyVec::zero(l);
    let mut s2 = PolyVec::zero(k);

    // Sample s1 (l polynomials)
    for i in 0..l {
        let mut hasher = Shake256::default();
        hasher.update(rho_prime);
        hasher.update(&(i as u16).to_le_bytes());
        let mut reader = hasher.finalize_xof();
        sample_eta_poly(&mut reader, &mut s1.polys[i], eta);
    }

    // Sample s2 (k polynomials)
    for i in 0..k {
        let mut hasher = Shake256::default();
        hasher.update(rho_prime);
        hasher.update(&((l + i) as u16).to_le_bytes());
        let mut reader = hasher.finalize_xof();
        sample_eta_poly(&mut reader, &mut s2.polys[i], eta);
    }

    (s1, s2)
}

/// Samples a polynomial with coefficients in [-η, η] using rejection sampling.
fn sample_eta_poly<R: XofReader>(reader: &mut R, poly: &mut Poly, eta: usize) {
    if eta == 2 {
        sample_eta2_poly(reader, poly);
    } else if eta == 4 {
        sample_eta4_poly(reader, poly);
    } else {
        panic!("Unsupported eta value: {}", eta);
    }
}

/// Samples a polynomial with coefficients in [-2, 2].
///
/// Each byte gives two coefficients using rejection from [0, 14] to [-2, 2].
fn sample_eta2_poly<R: XofReader>(reader: &mut R, poly: &mut Poly) {
    let mut coeffs_sampled = 0;
    let mut byte = [0u8; 1];

    while coeffs_sampled < N {
        reader.read(&mut byte);
        let b = byte[0];

        // Low nibble
        let t0 = b & 0x0F;
        if t0 < 15 && coeffs_sampled < N {
            // t0 mod 5 gives [0, 4], subtract 2 to get [-2, 2]
            poly.coeffs[coeffs_sampled] = (t0 % 5) as i32 - 2;
            coeffs_sampled += 1;
        }

        // High nibble
        let t1 = b >> 4;
        if t1 < 15 && coeffs_sampled < N {
            poly.coeffs[coeffs_sampled] = (t1 % 5) as i32 - 2;
            coeffs_sampled += 1;
        }
    }
}

/// Samples a polynomial with coefficients in [-4, 4].
///
/// Each byte gives two coefficients using rejection from [0, 8] to [-4, 4].
fn sample_eta4_poly<R: XofReader>(reader: &mut R, poly: &mut Poly) {
    let mut coeffs_sampled = 0;
    let mut byte = [0u8; 1];

    while coeffs_sampled < N {
        reader.read(&mut byte);
        let b = byte[0];

        // Low nibble
        let t0 = b & 0x0F;
        if t0 < 9 && coeffs_sampled < N {
            poly.coeffs[coeffs_sampled] = t0 as i32 - 4;
            coeffs_sampled += 1;
        }

        // High nibble
        let t1 = b >> 4;
        if t1 < 9 && coeffs_sampled < N {
            poly.coeffs[coeffs_sampled] = t1 as i32 - 4;
            coeffs_sampled += 1;
        }
    }
}

/// Expands seed into masking polynomial y with coefficients in [-γ1+1, γ1].
///
/// Uses SHAKE256 for deterministic expansion.
///
/// FIPS 204 Algorithm 28: ExpandMask
pub fn expand_mask(seed: &[u8], nonce: u16, gamma1: i32) -> Poly {
    let mut hasher = Shake256::default();
    hasher.update(seed);
    hasher.update(&nonce.to_le_bytes());
    let mut reader = hasher.finalize_xof();

    let mut poly = Poly::zero();

    if gamma1 == (1 << 17) {
        // γ1 = 2^17, need 18 bits per coefficient
        sample_gamma1_17(&mut reader, &mut poly);
    } else if gamma1 == (1 << 19) {
        // γ1 = 2^19, need 20 bits per coefficient
        sample_gamma1_19(&mut reader, &mut poly);
    } else {
        panic!("Unsupported gamma1 value: {}", gamma1);
    }

    poly
}

/// Samples coefficients in [-γ1+1, γ1] for γ1 = 2^17.
fn sample_gamma1_17<R: XofReader>(reader: &mut R, poly: &mut Poly) {
    let gamma1 = 1 << 17;
    let mut buf = [0u8; 576]; // 256 * 18 / 8 = 576 bytes
    reader.read(&mut buf);

    for i in 0..128 {
        // Each 9 bytes gives 4 coefficients (18 bits each)
        let base = i * 9 / 2;
        let offset = (i * 9) % 2;

        let (c0, c1) = if offset == 0 {
            let c0 = (buf[base] as i32)
                | ((buf[base + 1] as i32) << 8)
                | (((buf[base + 2] & 0x03) as i32) << 16);
            let c1 = ((buf[base + 2] as i32) >> 2)
                | ((buf[base + 3] as i32) << 6)
                | (((buf[base + 4] & 0x0F) as i32) << 14);
            (c0, c1)
        } else {
            let c0 = ((buf[base] as i32) >> 4)
                | ((buf[base + 1] as i32) << 4)
                | (((buf[base + 2] & 0x3F) as i32) << 12);
            let c1 = ((buf[base + 2] as i32) >> 6)
                | ((buf[base + 3] as i32) << 2)
                | ((buf[base + 4] as i32) << 10);
            (c0, c1)
        };

        poly.coeffs[2 * i] = gamma1 - c0;
        poly.coeffs[2 * i + 1] = gamma1 - c1;
    }
}

/// Samples coefficients in [-γ1+1, γ1] for γ1 = 2^19.
fn sample_gamma1_19<R: XofReader>(reader: &mut R, poly: &mut Poly) {
    let gamma1 = 1 << 19;
    let mut buf = [0u8; 640]; // 256 * 20 / 8 = 640 bytes
    reader.read(&mut buf);

    for i in 0..128 {
        // Each 5 bytes gives 2 coefficients (20 bits each)
        let base = i * 5;

        let c0 = (buf[base] as i32)
            | ((buf[base + 1] as i32) << 8)
            | (((buf[base + 2] & 0x0F) as i32) << 16);

        let c1 = ((buf[base + 2] as i32) >> 4)
            | ((buf[base + 3] as i32) << 4)
            | ((buf[base + 4] as i32) << 12);

        poly.coeffs[2 * i] = gamma1 - c0;
        poly.coeffs[2 * i + 1] = gamma1 - c1;
    }
}

/// Expands masking vector y (l polynomials).
pub fn expand_mask_vec(seed: &[u8], kappa: u16, l: usize, gamma1: i32) -> PolyVec {
    let mut y = PolyVec::zero(l);
    for i in 0..l {
        y.polys[i] = expand_mask(seed, kappa + i as u16, gamma1);
    }
    y
}

/// Samples challenge polynomial c with exactly τ non-zero coefficients in {-1, 1}.
///
/// Uses SHAKE256 for deterministic sampling.
///
/// FIPS 204 Algorithm 29: SampleInBall
pub fn sample_in_ball(seed: &[u8], tau: usize) -> Poly {
    let mut hasher = Shake256::default();
    hasher.update(seed);
    let mut reader = hasher.finalize_xof();

    let mut poly = Poly::zero();

    // First 8 bytes determine the signs
    let mut sign_bytes = [0u8; 8];
    reader.read(&mut sign_bytes);
    let signs = u64::from_le_bytes(sign_bytes);

    // Use rejection sampling to get τ distinct positions
    let mut byte = [0u8; 1];
    for i in (N - tau)..N {
        // Sample j ∈ [0, i]
        loop {
            reader.read(&mut byte);
            let j = byte[0] as usize;
            if j <= i {
                // Swap: c[i] = c[j], c[j] = ±1
                poly.coeffs[i] = poly.coeffs[j];

                let sign_bit = (signs >> (i - (N - tau))) & 1;
                poly.coeffs[j] = if sign_bit == 0 { 1 } else { Q - 1 }; // +1 or -1

                break;
            }
        }
    }

    poly
}

/// Hashes the public key to get tr (used in signing).
///
/// tr = H(ρ || t1_packed) with 64-byte output.
pub fn hash_public_key(pk_bytes: &[u8]) -> [u8; 64] {
    let mut hasher = Shake256::default();
    hasher.update(pk_bytes);
    let mut reader = hasher.finalize_xof();
    let mut tr = [0u8; 64];
    reader.read(&mut tr);
    tr
}

/// Hashes message with tr for signing.
///
/// μ = H(tr || M) with 64-byte output.
pub fn hash_message(tr: &[u8; 64], message: &[u8]) -> [u8; 64] {
    let mut hasher = Shake256::default();
    hasher.update(tr);
    hasher.update(message);
    let mut reader = hasher.finalize_xof();
    let mut mu = [0u8; 64];
    reader.read(&mut mu);
    mu
}

/// Computes challenge hash c_tilde from commitment and message hash.
///
/// c_tilde = H(μ || w1_packed) with λ/4 bytes output.
pub fn hash_challenge(mu: &[u8; 64], w1_bytes: &[u8], lambda: usize) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(mu);
    hasher.update(w1_bytes);
    let mut reader = hasher.finalize_xof();

    let len = lambda / 4;
    let mut c_tilde = vec![0u8; len];
    reader.read(&mut c_tilde);
    c_tilde
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ML_DSA_44, ML_DSA_65};

    #[test]
    fn test_expand_a_dimensions() {
        let rho = [0u8; 32];
        let params = &ML_DSA_44;
        let a = expand_a(&rho, params);

        assert_eq!(a.nrows, params.k);
        assert_eq!(a.ncols, params.l);

        for row in &a.rows {
            assert_eq!(row.len(), params.l);
            for poly in &row.polys {
                // All coefficients should be in [0, Q)
                for &c in &poly.coeffs {
                    assert!(c >= 0 && c < Q);
                }
            }
        }
    }

    #[test]
    fn test_expand_a_deterministic() {
        let rho = [42u8; 32];
        let params = &ML_DSA_44;

        let a1 = expand_a(&rho, params);
        let a2 = expand_a(&rho, params);

        // Same seed should produce same matrix
        for i in 0..params.k {
            for j in 0..params.l {
                assert_eq!(a1.rows[i].polys[j].coeffs, a2.rows[i].polys[j].coeffs);
            }
        }
    }

    #[test]
    fn test_expand_s_dimensions() {
        let rho_prime = [0u8; 64];
        let params = &ML_DSA_44;

        let (s1, s2) = expand_s(&rho_prime, params);

        assert_eq!(s1.len(), params.l);
        assert_eq!(s2.len(), params.k);
    }

    #[test]
    fn test_expand_s_eta_bounds() {
        let rho_prime = [99u8; 64];
        let params = &ML_DSA_44;
        let eta = params.eta as i32;

        let (s1, s2) = expand_s(&rho_prime, params);

        // All coefficients should be in [-η, η]
        for poly in &s1.polys {
            for &c in &poly.coeffs {
                assert!(c >= -eta && c <= eta, "s1 coeff {} out of bounds", c);
            }
        }
        for poly in &s2.polys {
            for &c in &poly.coeffs {
                assert!(c >= -eta && c <= eta, "s2 coeff {} out of bounds", c);
            }
        }
    }

    #[test]
    fn test_expand_s_eta4() {
        let rho_prime = [123u8; 64];
        let params = &ML_DSA_65; // Uses η = 4

        let (s1, s2) = expand_s(&rho_prime, params);

        for poly in s1.polys.iter().chain(s2.polys.iter()) {
            for &c in &poly.coeffs {
                assert!(c >= -4 && c <= 4, "coeff {} out of [-4, 4]", c);
            }
        }
    }

    #[test]
    fn test_sample_in_ball_weight() {
        let seed = [0u8; 32];
        let tau = 39;

        let c = sample_in_ball(&seed, tau);

        // Count non-zero coefficients
        let weight: usize = c.coeffs.iter().filter(|&&x| x != 0).count();
        assert_eq!(weight, tau);

        // All non-zero coefficients should be ±1
        for &coeff in &c.coeffs {
            assert!(coeff == 0 || coeff == 1 || coeff == Q - 1);
        }
    }

    #[test]
    fn test_sample_in_ball_deterministic() {
        let seed = [42u8; 32];
        let tau = 49;

        let c1 = sample_in_ball(&seed, tau);
        let c2 = sample_in_ball(&seed, tau);

        assert_eq!(c1.coeffs, c2.coeffs);
    }

    #[test]
    fn test_expand_mask_bounds() {
        let seed = [0u8; 66];
        let gamma1 = 1 << 17;

        let y = expand_mask(&seed, 0, gamma1);

        // Coefficients should be in [-γ1+1, γ1]
        for &c in &y.coeffs {
            assert!(c >= -gamma1 + 1 && c <= gamma1, "coeff {} out of bounds", c);
        }
    }

    #[test]
    fn test_hash_functions_output_sizes() {
        let data = b"test data";

        let tr = hash_public_key(data);
        assert_eq!(tr.len(), 64);

        let mu = hash_message(&tr, data);
        assert_eq!(mu.len(), 64);

        let c_tilde_128 = hash_challenge(&mu, data, 128);
        assert_eq!(c_tilde_128.len(), 32);

        let c_tilde_256 = hash_challenge(&mu, data, 256);
        assert_eq!(c_tilde_256.len(), 64);
    }
}
