//! Bit packing and serialization for ML-DSA.
//!
//! This module implements the packing functions for serializing and deserializing
//! public keys, secret keys, and signatures.
//!
//! Packing formats:
//! - t1: 10 bits per coefficient (coefficients in [0, 2^10))
//! - t0: 13 bits per coefficient (coefficients in [-(2^12-1), 2^12])
//! - s1/s2: 3 bits (η=2) or 4 bits (η=4) per coefficient
//! - z: 18 bits (γ1=2^17) or 20 bits (γ1=2^19) per coefficient
//! - w1: 6 bits (γ2=95232) or 4 bits (γ2=261888) per coefficient
//! - hints: sparse encoding with indices

use crate::error::{MlDsaError, Result};
use crate::params::{D, N};
use crate::poly::Poly;
use crate::polyvec::PolyVec;

/// Packs t1 polynomial (10 bits per coefficient).
///
/// t1 coefficients are in [0, 2^10).
pub fn pack_t1(poly: &Poly) -> Vec<u8> {
    let mut bytes = vec![0u8; N * 10 / 8]; // 320 bytes per polynomial

    for i in 0..N / 4 {
        let c0 = poly.coeffs[4 * i] as u32;
        let c1 = poly.coeffs[4 * i + 1] as u32;
        let c2 = poly.coeffs[4 * i + 2] as u32;
        let c3 = poly.coeffs[4 * i + 3] as u32;

        bytes[5 * i] = c0 as u8;
        bytes[5 * i + 1] = ((c0 >> 8) | (c1 << 2)) as u8;
        bytes[5 * i + 2] = ((c1 >> 6) | (c2 << 4)) as u8;
        bytes[5 * i + 3] = ((c2 >> 4) | (c3 << 6)) as u8;
        bytes[5 * i + 4] = (c3 >> 2) as u8;
    }

    bytes
}

/// Unpacks t1 polynomial (10 bits per coefficient).
///
/// Validates that all coefficients are in the valid range [0, 2^10).
pub fn unpack_t1(bytes: &[u8]) -> Result<Poly> {
    if bytes.len() != N * 10 / 8 {
        return Err(MlDsaError::DecodingError { context: "t1" });
    }

    let mut poly = Poly::zero();

    for i in 0..N / 4 {
        poly.coeffs[4 * i] =
            (bytes[5 * i] as i32) | (((bytes[5 * i + 1] & 0x03) as i32) << 8);
        poly.coeffs[4 * i + 1] =
            ((bytes[5 * i + 1] >> 2) as i32) | (((bytes[5 * i + 2] & 0x0F) as i32) << 6);
        poly.coeffs[4 * i + 2] =
            ((bytes[5 * i + 2] >> 4) as i32) | (((bytes[5 * i + 3] & 0x3F) as i32) << 4);
        poly.coeffs[4 * i + 3] =
            ((bytes[5 * i + 3] >> 6) as i32) | ((bytes[5 * i + 4] as i32) << 2);
    }

    // Validate coefficient range: t1 must be in [0, 2^10)
    let max_t1 = 1 << 10;
    for coeff in &poly.coeffs {
        if *coeff < 0 || *coeff >= max_t1 {
            return Err(MlDsaError::DecodingError { context: "t1 coefficient out of range" });
        }
    }

    Ok(poly)
}

/// Packs t0 polynomial (13 bits per coefficient, centered).
///
/// t0 coefficients are in [-(2^(D-1)-1), 2^(D-1)].
pub fn pack_t0(poly: &Poly) -> Vec<u8> {
    let mut bytes = vec![0u8; N * D / 8]; // 416 bytes per polynomial
    let center = 1 << (D - 1); // 2^12 = 4096

    for i in 0..N / 8 {
        // Pack 8 coefficients into 13 bytes
        let mut coeffs = [0u32; 8];
        for j in 0..8 {
            // Convert from centered to unsigned: add 2^(D-1)
            coeffs[j] = (center - poly.coeffs[8 * i + j]) as u32;
        }

        bytes[13 * i] = coeffs[0] as u8;
        bytes[13 * i + 1] = ((coeffs[0] >> 8) | (coeffs[1] << 5)) as u8;
        bytes[13 * i + 2] = (coeffs[1] >> 3) as u8;
        bytes[13 * i + 3] = ((coeffs[1] >> 11) | (coeffs[2] << 2)) as u8;
        bytes[13 * i + 4] = ((coeffs[2] >> 6) | (coeffs[3] << 7)) as u8;
        bytes[13 * i + 5] = (coeffs[3] >> 1) as u8;
        bytes[13 * i + 6] = ((coeffs[3] >> 9) | (coeffs[4] << 4)) as u8;
        bytes[13 * i + 7] = (coeffs[4] >> 4) as u8;
        bytes[13 * i + 8] = ((coeffs[4] >> 12) | (coeffs[5] << 1)) as u8;
        bytes[13 * i + 9] = ((coeffs[5] >> 7) | (coeffs[6] << 6)) as u8;
        bytes[13 * i + 10] = (coeffs[6] >> 2) as u8;
        bytes[13 * i + 11] = ((coeffs[6] >> 10) | (coeffs[7] << 3)) as u8;
        bytes[13 * i + 12] = (coeffs[7] >> 5) as u8;
    }

    bytes
}

/// Unpacks t0 polynomial (13 bits per coefficient).
///
/// Validates that all coefficients are in the valid range [-(2^(D-1)-1), 2^(D-1)].
pub fn unpack_t0(bytes: &[u8]) -> Result<Poly> {
    if bytes.len() != N * D / 8 {
        return Err(MlDsaError::DecodingError { context: "t0" });
    }

    let mut poly = Poly::zero();
    let center = 1 << (D - 1);

    for i in 0..N / 8 {
        let b = &bytes[13 * i..13 * i + 13];

        let c0 = (b[0] as u32) | (((b[1] & 0x1F) as u32) << 8);
        let c1 = ((b[1] >> 5) as u32) | ((b[2] as u32) << 3) | (((b[3] & 0x03) as u32) << 11);
        let c2 = ((b[3] >> 2) as u32) | (((b[4] & 0x7F) as u32) << 6);
        let c3 = ((b[4] >> 7) as u32) | ((b[5] as u32) << 1) | (((b[6] & 0x0F) as u32) << 9);
        let c4 = ((b[6] >> 4) as u32) | ((b[7] as u32) << 4) | (((b[8] & 0x01) as u32) << 12);
        let c5 = ((b[8] >> 1) as u32) | (((b[9] & 0x3F) as u32) << 7);
        let c6 = ((b[9] >> 6) as u32) | ((b[10] as u32) << 2) | (((b[11] & 0x07) as u32) << 10);
        let c7 = ((b[11] >> 3) as u32) | ((b[12] as u32) << 5);

        poly.coeffs[8 * i] = (center as i32) - (c0 as i32);
        poly.coeffs[8 * i + 1] = (center as i32) - (c1 as i32);
        poly.coeffs[8 * i + 2] = (center as i32) - (c2 as i32);
        poly.coeffs[8 * i + 3] = (center as i32) - (c3 as i32);
        poly.coeffs[8 * i + 4] = (center as i32) - (c4 as i32);
        poly.coeffs[8 * i + 5] = (center as i32) - (c5 as i32);
        poly.coeffs[8 * i + 6] = (center as i32) - (c6 as i32);
        poly.coeffs[8 * i + 7] = (center as i32) - (c7 as i32);
    }

    // Validate coefficient range: t0 must be in [-(2^(D-1)-1), 2^(D-1)]
    let min_t0 = -(center as i32) + 1;
    let max_t0 = center as i32;
    for coeff in &poly.coeffs {
        if *coeff < min_t0 || *coeff > max_t0 {
            return Err(MlDsaError::DecodingError { context: "t0 coefficient out of range" });
        }
    }

    Ok(poly)
}

/// Packs eta polynomial (η=2: 3 bits, η=4: 4 bits per coefficient).
///
/// Returns an error if eta is not 2 or 4.
pub fn pack_eta(poly: &Poly, eta: usize) -> Result<Vec<u8>> {
    match eta {
        2 => Ok(pack_eta2(poly)),
        4 => Ok(pack_eta4(poly)),
        _ => Err(MlDsaError::InvalidParams {
            reason: "eta must be 2 or 4",
        }),
    }
}

/// Packs polynomial with η=2 (3 bits per coefficient).
fn pack_eta2(poly: &Poly) -> Vec<u8> {
    let mut bytes = vec![0u8; N * 3 / 8]; // 96 bytes

    for i in 0..N / 8 {
        // Convert from [-2, 2] to [0, 4]
        let mut t = [0u8; 8];
        for j in 0..8 {
            t[j] = (2 - poly.coeffs[8 * i + j]) as u8;
        }

        bytes[3 * i] = t[0] | (t[1] << 3) | (t[2] << 6);
        bytes[3 * i + 1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
        bytes[3 * i + 2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    }

    bytes
}

/// Unpacks polynomial with η=2.
///
/// Validates that all coefficients are in the valid range [-2, 2].
pub fn unpack_eta2(bytes: &[u8]) -> Result<Poly> {
    if bytes.len() != N * 3 / 8 {
        return Err(MlDsaError::DecodingError { context: "eta2" });
    }

    let mut poly = Poly::zero();

    for i in 0..N / 8 {
        let b = &bytes[3 * i..3 * i + 3];

        let t0 = b[0] & 0x07;
        let t1 = (b[0] >> 3) & 0x07;
        let t2 = ((b[0] >> 6) | (b[1] << 2)) & 0x07;
        let t3 = (b[1] >> 1) & 0x07;
        let t4 = (b[1] >> 4) & 0x07;
        let t5 = ((b[1] >> 7) | (b[2] << 1)) & 0x07;
        let t6 = (b[2] >> 2) & 0x07;
        let t7 = b[2] >> 5;

        // Validate raw values are in [0, 4] before conversion
        if t0 > 4 || t1 > 4 || t2 > 4 || t3 > 4 || t4 > 4 || t5 > 4 || t6 > 4 || t7 > 4 {
            return Err(MlDsaError::DecodingError { context: "eta2 coefficient out of range" });
        }

        poly.coeffs[8 * i] = 2 - t0 as i32;
        poly.coeffs[8 * i + 1] = 2 - t1 as i32;
        poly.coeffs[8 * i + 2] = 2 - t2 as i32;
        poly.coeffs[8 * i + 3] = 2 - t3 as i32;
        poly.coeffs[8 * i + 4] = 2 - t4 as i32;
        poly.coeffs[8 * i + 5] = 2 - t5 as i32;
        poly.coeffs[8 * i + 6] = 2 - t6 as i32;
        poly.coeffs[8 * i + 7] = 2 - t7 as i32;
    }

    Ok(poly)
}

/// Packs polynomial with η=4 (4 bits per coefficient).
fn pack_eta4(poly: &Poly) -> Vec<u8> {
    let mut bytes = vec![0u8; N / 2]; // 128 bytes

    for i in 0..N / 2 {
        // Convert from [-4, 4] to [0, 8]
        let t0 = (4 - poly.coeffs[2 * i]) as u8;
        let t1 = (4 - poly.coeffs[2 * i + 1]) as u8;
        bytes[i] = t0 | (t1 << 4);
    }

    bytes
}

/// Unpacks polynomial with η=4.
///
/// Validates that all coefficients are in the valid range [-4, 4].
pub fn unpack_eta4(bytes: &[u8]) -> Result<Poly> {
    if bytes.len() != N / 2 {
        return Err(MlDsaError::DecodingError { context: "eta4" });
    }

    let mut poly = Poly::zero();

    for i in 0..N / 2 {
        let t0 = bytes[i] & 0x0F;
        let t1 = bytes[i] >> 4;

        // Validate raw values are in [0, 8] before conversion
        if t0 > 8 || t1 > 8 {
            return Err(MlDsaError::DecodingError { context: "eta4 coefficient out of range" });
        }

        poly.coeffs[2 * i] = 4 - t0 as i32;
        poly.coeffs[2 * i + 1] = 4 - t1 as i32;
    }

    Ok(poly)
}

/// Packs z polynomial (γ1 determines bit width).
///
/// Returns an error if gamma1 is not 2^17 or 2^19.
pub fn pack_z(poly: &Poly, gamma1: i32) -> Result<Vec<u8>> {
    match gamma1 {
        g if g == (1 << 17) => Ok(pack_z_17(poly)),
        g if g == (1 << 19) => Ok(pack_z_19(poly)),
        _ => Err(MlDsaError::InvalidParams {
            reason: "gamma1 must be 2^17 or 2^19",
        }),
    }
}

/// Packs z with γ1=2^17 (18 bits per coefficient).
fn pack_z_17(poly: &Poly) -> Vec<u8> {
    let gamma1 = 1 << 17;
    let mut bytes = vec![0u8; N * 18 / 8]; // 576 bytes

    for i in 0..N / 4 {
        // Convert from [-(γ1-1), γ1] to [0, 2γ1-1]
        let c0 = (gamma1 - poly.coeffs[4 * i]) as u32;
        let c1 = (gamma1 - poly.coeffs[4 * i + 1]) as u32;
        let c2 = (gamma1 - poly.coeffs[4 * i + 2]) as u32;
        let c3 = (gamma1 - poly.coeffs[4 * i + 3]) as u32;

        bytes[9 * i] = c0 as u8;
        bytes[9 * i + 1] = (c0 >> 8) as u8;
        bytes[9 * i + 2] = ((c0 >> 16) | (c1 << 2)) as u8;
        bytes[9 * i + 3] = (c1 >> 6) as u8;
        bytes[9 * i + 4] = ((c1 >> 14) | (c2 << 4)) as u8;
        bytes[9 * i + 5] = (c2 >> 4) as u8;
        bytes[9 * i + 6] = ((c2 >> 12) | (c3 << 6)) as u8;
        bytes[9 * i + 7] = (c3 >> 2) as u8;
        bytes[9 * i + 8] = (c3 >> 10) as u8;
    }

    bytes
}

/// Packs z with γ1=2^19 (20 bits per coefficient).
fn pack_z_19(poly: &Poly) -> Vec<u8> {
    let gamma1 = 1 << 19;
    let mut bytes = vec![0u8; N * 20 / 8]; // 640 bytes

    for i in 0..N / 4 {
        let c0 = (gamma1 - poly.coeffs[4 * i]) as u32;
        let c1 = (gamma1 - poly.coeffs[4 * i + 1]) as u32;
        let c2 = (gamma1 - poly.coeffs[4 * i + 2]) as u32;
        let c3 = (gamma1 - poly.coeffs[4 * i + 3]) as u32;

        bytes[10 * i] = c0 as u8;
        bytes[10 * i + 1] = (c0 >> 8) as u8;
        bytes[10 * i + 2] = ((c0 >> 16) | (c1 << 4)) as u8;
        bytes[10 * i + 3] = (c1 >> 4) as u8;
        bytes[10 * i + 4] = (c1 >> 12) as u8;
        bytes[10 * i + 5] = c2 as u8;
        bytes[10 * i + 6] = (c2 >> 8) as u8;
        bytes[10 * i + 7] = ((c2 >> 16) | (c3 << 4)) as u8;
        bytes[10 * i + 8] = (c3 >> 4) as u8;
        bytes[10 * i + 9] = (c3 >> 12) as u8;
    }

    bytes
}

/// Unpacks z polynomial.
pub fn unpack_z(bytes: &[u8], gamma1: i32) -> Result<Poly> {
    if gamma1 == (1 << 17) {
        unpack_z_17(bytes)
    } else if gamma1 == (1 << 19) {
        unpack_z_19(bytes)
    } else {
        Err(MlDsaError::InvalidParams {
            reason: "unsupported gamma1",
        })
    }
}

fn unpack_z_17(bytes: &[u8]) -> Result<Poly> {
    if bytes.len() != N * 18 / 8 {
        return Err(MlDsaError::DecodingError { context: "z" });
    }

    let gamma1 = 1 << 17;
    let max_encoded = (2 * gamma1 - 1) as u32; // 2γ1 - 1
    let mut poly = Poly::zero();

    for i in 0..N / 4 {
        let b = &bytes[9 * i..9 * i + 9];

        let c0 = (b[0] as u32) | ((b[1] as u32) << 8) | (((b[2] & 0x03) as u32) << 16);
        let c1 = ((b[2] >> 2) as u32) | ((b[3] as u32) << 6) | (((b[4] & 0x0F) as u32) << 14);
        let c2 = ((b[4] >> 4) as u32) | ((b[5] as u32) << 4) | (((b[6] & 0x3F) as u32) << 12);
        let c3 = ((b[6] >> 6) as u32) | ((b[7] as u32) << 2) | ((b[8] as u32) << 10);

        // Validate encoded values are in [0, 2γ1 - 1]
        if c0 > max_encoded || c1 > max_encoded || c2 > max_encoded || c3 > max_encoded {
            return Err(MlDsaError::DecodingError { context: "z coefficient out of range" });
        }

        poly.coeffs[4 * i] = gamma1 - c0 as i32;
        poly.coeffs[4 * i + 1] = gamma1 - c1 as i32;
        poly.coeffs[4 * i + 2] = gamma1 - c2 as i32;
        poly.coeffs[4 * i + 3] = gamma1 - c3 as i32;
    }

    Ok(poly)
}

fn unpack_z_19(bytes: &[u8]) -> Result<Poly> {
    if bytes.len() != N * 20 / 8 {
        return Err(MlDsaError::DecodingError { context: "z" });
    }

    let gamma1 = 1 << 19;
    let max_encoded = (2 * gamma1 - 1) as u32; // 2γ1 - 1
    let mut poly = Poly::zero();

    for i in 0..N / 4 {
        let b = &bytes[10 * i..10 * i + 10];

        let c0 = (b[0] as u32) | ((b[1] as u32) << 8) | (((b[2] & 0x0F) as u32) << 16);
        let c1 = ((b[2] >> 4) as u32) | ((b[3] as u32) << 4) | ((b[4] as u32) << 12);
        let c2 = (b[5] as u32) | ((b[6] as u32) << 8) | (((b[7] & 0x0F) as u32) << 16);
        let c3 = ((b[7] >> 4) as u32) | ((b[8] as u32) << 4) | ((b[9] as u32) << 12);

        // Validate encoded values are in [0, 2γ1 - 1]
        if c0 > max_encoded || c1 > max_encoded || c2 > max_encoded || c3 > max_encoded {
            return Err(MlDsaError::DecodingError { context: "z coefficient out of range" });
        }

        poly.coeffs[4 * i] = gamma1 - c0 as i32;
        poly.coeffs[4 * i + 1] = gamma1 - c1 as i32;
        poly.coeffs[4 * i + 2] = gamma1 - c2 as i32;
        poly.coeffs[4 * i + 3] = gamma1 - c3 as i32;
    }

    Ok(poly)
}

/// Packs w1 polynomial for commitment.
///
/// Returns an error if gamma2 is not 95232 or 261888.
pub fn pack_w1(poly: &Poly, gamma2: i32) -> Result<Vec<u8>> {
    match gamma2 {
        95232 => Ok(pack_w1_6bit(poly)),    // (q-1)/88, 6 bits per coefficient
        261888 => Ok(pack_w1_4bit(poly)),   // (q-1)/32, 4 bits per coefficient
        _ => Err(MlDsaError::InvalidParams {
            reason: "gamma2 must be 95232 or 261888",
        }),
    }
}

fn pack_w1_6bit(poly: &Poly) -> Vec<u8> {
    let mut bytes = vec![0u8; N * 6 / 8]; // 192 bytes

    for i in 0..N / 4 {
        bytes[3 * i] = (poly.coeffs[4 * i] as u8) | ((poly.coeffs[4 * i + 1] as u8) << 6);
        bytes[3 * i + 1] =
            ((poly.coeffs[4 * i + 1] >> 2) as u8) | ((poly.coeffs[4 * i + 2] as u8) << 4);
        bytes[3 * i + 2] =
            ((poly.coeffs[4 * i + 2] >> 4) as u8) | ((poly.coeffs[4 * i + 3] as u8) << 2);
    }

    bytes
}

fn pack_w1_4bit(poly: &Poly) -> Vec<u8> {
    let mut bytes = vec![0u8; N / 2]; // 128 bytes

    for i in 0..N / 2 {
        bytes[i] = (poly.coeffs[2 * i] as u8) | ((poly.coeffs[2 * i + 1] as u8) << 4);
    }

    bytes
}

/// Packs hint vector using sparse encoding.
///
/// Format: For each polynomial, store indices where hint[i] = 1,
/// followed by count of hints for that polynomial.
pub fn pack_hints(hints: &[Vec<bool>], omega: usize, k: usize) -> Vec<u8> {
    let mut bytes = vec![0u8; omega + k];
    let mut index = 0;

    for i in 0..k {
        for j in 0..N {
            if hints[i][j] {
                bytes[index] = j as u8;
                index += 1;
            }
        }
        bytes[omega + i] = index as u8;
    }

    bytes
}

/// Unpacks hint vector.
pub fn unpack_hints(bytes: &[u8], omega: usize, k: usize) -> Result<Vec<Vec<bool>>> {
    if bytes.len() != omega + k {
        return Err(MlDsaError::DecodingError { context: "hints" });
    }

    let mut hints = vec![vec![false; N]; k];
    let mut index = 0;

    for i in 0..k {
        let end = bytes[omega + i] as usize;
        if end < index || end > omega {
            return Err(MlDsaError::InvalidHint);
        }

        while index < end {
            let j = bytes[index] as usize;
            if j >= N {
                return Err(MlDsaError::InvalidHint);
            }
            if hints[i][j] {
                return Err(MlDsaError::InvalidHint); // Duplicate
            }
            hints[i][j] = true;
            index += 1;
        }
    }

    // Remaining positions should be zero
    for i in index..omega {
        if bytes[i] != 0 {
            return Err(MlDsaError::InvalidHint);
        }
    }

    Ok(hints)
}

/// Packs a vector of t1 polynomials.
pub fn pack_t1_vec(v: &PolyVec) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(v.len() * N * 10 / 8);
    for poly in &v.polys {
        bytes.extend(pack_t1(poly));
    }
    bytes
}

/// Packs a vector of t0 polynomials.
pub fn pack_t0_vec(v: &PolyVec) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(v.len() * N * D / 8);
    for poly in &v.polys {
        bytes.extend(pack_t0(poly));
    }
    bytes
}

/// Packs a vector of eta polynomials.
pub fn pack_eta_vec(v: &PolyVec, eta: usize) -> Result<Vec<u8>> {
    let bytes_per_poly = if eta == 2 { N * 3 / 8 } else { N / 2 };
    let mut bytes = Vec::with_capacity(v.len() * bytes_per_poly);
    for poly in &v.polys {
        bytes.extend(pack_eta(poly, eta)?);
    }
    Ok(bytes)
}

/// Packs a vector of z polynomials.
pub fn pack_z_vec(v: &PolyVec, gamma1: i32) -> Result<Vec<u8>> {
    let bits = if gamma1 == (1 << 17) { 18 } else { 20 };
    let bytes_per_poly = N * bits / 8;
    let mut bytes = Vec::with_capacity(v.len() * bytes_per_poly);
    for poly in &v.polys {
        bytes.extend(pack_z(poly, gamma1)?);
    }
    Ok(bytes)
}

/// Packs a vector of w1 polynomials.
pub fn pack_w1_vec(v: &PolyVec, gamma2: i32) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    for poly in &v.polys {
        bytes.extend(pack_w1(poly, gamma2)?);
    }
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_unpack_t1() {
        let mut poly = Poly::zero();
        for i in 0..N {
            poly.coeffs[i] = (i * 3) as i32 % (1 << 10);
        }

        let packed = pack_t1(&poly);
        let unpacked = unpack_t1(&packed).unwrap();

        for i in 0..N {
            assert_eq!(poly.coeffs[i], unpacked.coeffs[i], "mismatch at {}", i);
        }
    }

    #[test]
    fn test_pack_unpack_t0() {
        let mut poly = Poly::zero();
        let center = 1 << (D - 1); // 4096

        // t0 coefficients are in (-2^(D-1), 2^(D-1)] = (-4096, 4096]
        // Using range [-4095, 4096] for testing
        for i in 0..N {
            // Generate values in [-4095, 4095] (avoiding the excluded boundary -4096)
            poly.coeffs[i] = ((i as i32 * 7) % (2 * center - 1)) - center + 1;
        }

        let packed = pack_t0(&poly);
        let unpacked = unpack_t0(&packed).unwrap();

        for i in 0..N {
            assert_eq!(poly.coeffs[i], unpacked.coeffs[i], "mismatch at {}", i);
        }
    }

    #[test]
    fn test_pack_unpack_eta2() {
        let mut poly = Poly::zero();
        for i in 0..N {
            poly.coeffs[i] = (i as i32 % 5) - 2; // [-2, 2]
        }

        let packed = pack_eta(&poly, 2).unwrap();
        let unpacked = unpack_eta2(&packed).unwrap();

        for i in 0..N {
            assert_eq!(poly.coeffs[i], unpacked.coeffs[i], "mismatch at {}", i);
        }
    }

    #[test]
    fn test_pack_unpack_eta4() {
        let mut poly = Poly::zero();
        for i in 0..N {
            poly.coeffs[i] = (i as i32 % 9) - 4; // [-4, 4]
        }

        let packed = pack_eta(&poly, 4).unwrap();
        let unpacked = unpack_eta4(&packed).unwrap();

        for i in 0..N {
            assert_eq!(poly.coeffs[i], unpacked.coeffs[i], "mismatch at {}", i);
        }
    }

    #[test]
    fn test_pack_unpack_z_17() {
        let gamma1 = 1 << 17;
        let mut poly = Poly::zero();
        for i in 0..N {
            poly.coeffs[i] = ((i as i32 * 1000) % (2 * gamma1)) - gamma1 + 1;
        }

        let packed = pack_z(&poly, gamma1).unwrap();
        let unpacked = unpack_z(&packed, gamma1).unwrap();

        for i in 0..N {
            assert_eq!(poly.coeffs[i], unpacked.coeffs[i], "mismatch at {}", i);
        }
    }

    #[test]
    fn test_pack_unpack_z_19() {
        let gamma1 = 1 << 19;
        let mut poly = Poly::zero();
        for i in 0..N {
            poly.coeffs[i] = ((i as i32 * 2000) % (2 * gamma1)) - gamma1 + 1;
        }

        let packed = pack_z(&poly, gamma1).unwrap();
        let unpacked = unpack_z(&packed, gamma1).unwrap();

        for i in 0..N {
            assert_eq!(poly.coeffs[i], unpacked.coeffs[i], "mismatch at {}", i);
        }
    }

    #[test]
    fn test_pack_unpack_hints() {
        let k = 4;
        let omega = 80;

        let mut hints = vec![vec![false; N]; k];
        hints[0][0] = true;
        hints[0][100] = true;
        hints[1][50] = true;
        hints[2][200] = true;
        hints[3][255] = true;

        let packed = pack_hints(&hints, omega, k);
        let unpacked = unpack_hints(&packed, omega, k).unwrap();

        for i in 0..k {
            for j in 0..N {
                assert_eq!(hints[i][j], unpacked[i][j], "mismatch at [{}, {}]", i, j);
            }
        }
    }

    #[test]
    fn test_hints_invalid_count() {
        let k = 4;
        let omega = 80;

        // Create invalid hint encoding
        let mut bytes = vec![0u8; omega + k];
        bytes[omega] = 100; // Count exceeds omega

        assert!(unpack_hints(&bytes, omega, k).is_err());
    }
}
