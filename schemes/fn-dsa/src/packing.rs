//! Serialization and deserialization for FALCON keys and signatures.
//!
//! This module implements the encoding/decoding of:
//! - Public keys (h polynomial)
//! - Secret keys (f, g, F, G polynomials)
//! - Signatures (nonce and compressed s2)
//!
//! The encoding follows a simplified format suitable for educational purposes.
//! Production implementations should use the FALCON specification's encoding.

use crate::error::{FnDsaError, Result};
use crate::fft::Complex;
use crate::fft_tree::GramSchmidt;
use crate::keygen::{KeyPair, PublicKey, SecretKey};
use crate::params::{FALCON_512, FALCON_1024, NONCE_SIZE};
use crate::sign::Signature;

// ============================================================================
// Constants
// ============================================================================

/// Magic bytes for FALCON public key.
const PK_MAGIC: [u8; 4] = *b"FPK1";

/// Magic bytes for FALCON secret key.
const SK_MAGIC: [u8; 4] = *b"FSK1";

/// Magic bytes for FALCON signature.
const SIG_MAGIC: [u8; 4] = *b"FSG1";

// ============================================================================
// Public Key Encoding
// ============================================================================

/// Encodes a public key to bytes.
///
/// Format:
/// - 4 bytes: magic "FPK1"
/// - 2 bytes: n (little-endian)
/// - n * 2 bytes: h coefficients (little-endian i16)
pub fn encode_public_key(pk: &PublicKey) -> Vec<u8> {
    let n = pk.params.n;
    let mut bytes = Vec::with_capacity(4 + 2 + n * 2);

    // Magic
    bytes.extend_from_slice(&PK_MAGIC);

    // n
    bytes.extend_from_slice(&(n as u16).to_le_bytes());

    // h coefficients
    for &coeff in &pk.h {
        bytes.extend_from_slice(&coeff.to_le_bytes());
    }

    bytes
}

/// Decodes a public key from bytes.
pub fn decode_public_key(bytes: &[u8]) -> Result<PublicKey> {
    if bytes.len() < 6 {
        return Err(FnDsaError::InvalidInput {
            field: "public_key",
            reason: "too short",
        });
    }

    // Check magic
    if &bytes[0..4] != &PK_MAGIC {
        return Err(FnDsaError::InvalidInput {
            field: "public_key",
            reason: "invalid magic bytes",
        });
    }

    // Read n
    let n = u16::from_le_bytes([bytes[4], bytes[5]]) as usize;

    // Determine params
    let params = match n {
        512 => FALCON_512,
        1024 => FALCON_1024,
        _ => {
            return Err(FnDsaError::InvalidInput {
                field: "public_key",
                reason: "unsupported n value",
            });
        }
    };

    // Check length
    let expected_len = 6 + n * 2;
    if bytes.len() < expected_len {
        return Err(FnDsaError::InvalidInput {
            field: "public_key",
            reason: "incomplete h coefficients",
        });
    }

    // Read h coefficients
    let mut h = Vec::with_capacity(n);
    for i in 0..n {
        let offset = 6 + i * 2;
        let coeff = i16::from_le_bytes([bytes[offset], bytes[offset + 1]]);
        h.push(coeff);
    }

    // Validate h coefficients are in valid range [-(Q-1)/2, (Q-1)/2] (centered) or [0, Q) (positive)
    // Public keys can use either representation, so accept both
    let q = crate::params::Q as i16;
    for &coeff in &h {
        if coeff < -(q - 1) || coeff >= q {
            return Err(FnDsaError::InvalidInput {
                field: "public_key",
                reason: "h coefficient out of range",
            });
        }
    }

    Ok(PublicKey { h, params })
}

// ============================================================================
// Secret Key Encoding
// ============================================================================

/// Encodes a secret key to bytes.
///
/// Format:
/// - 4 bytes: magic "FSK1"
/// - 2 bytes: n (little-endian)
/// - n bytes: f coefficients (i8)
/// - n bytes: g coefficients (i8)
/// - n * 2 bytes: F coefficients (little-endian i16)
/// - n * 2 bytes: G coefficients (little-endian i16)
///
/// Note: The Gram-Schmidt data is recomputed on decode rather than stored.
pub fn encode_secret_key(sk: &SecretKey) -> Vec<u8> {
    let n = sk.params.n;
    let mut bytes = Vec::with_capacity(4 + 2 + n + n + n * 2 + n * 2);

    // Magic
    bytes.extend_from_slice(&SK_MAGIC);

    // n
    bytes.extend_from_slice(&(n as u16).to_le_bytes());

    // f coefficients (i8)
    for &coeff in &sk.f {
        bytes.push(coeff as u8);
    }

    // g coefficients (i8)
    for &coeff in &sk.g {
        bytes.push(coeff as u8);
    }

    // F coefficients (i16)
    for &coeff in &sk.big_f {
        bytes.extend_from_slice(&coeff.to_le_bytes());
    }

    // G coefficients (i16)
    for &coeff in &sk.big_g {
        bytes.extend_from_slice(&coeff.to_le_bytes());
    }

    bytes
}

/// Decodes a secret key from bytes.
pub fn decode_secret_key(bytes: &[u8]) -> Result<SecretKey> {
    if bytes.len() < 6 {
        return Err(FnDsaError::InvalidInput {
            field: "secret_key",
            reason: "too short",
        });
    }

    // Check magic
    if &bytes[0..4] != &SK_MAGIC {
        return Err(FnDsaError::InvalidInput {
            field: "secret_key",
            reason: "invalid magic bytes",
        });
    }

    // Read n
    let n = u16::from_le_bytes([bytes[4], bytes[5]]) as usize;

    // Determine params
    let params = match n {
        512 => FALCON_512,
        1024 => FALCON_1024,
        _ => {
            return Err(FnDsaError::InvalidInput {
                field: "secret_key",
                reason: "unsupported n value",
            });
        }
    };

    // Check length
    let expected_len = 6 + n + n + n * 2 + n * 2;
    if bytes.len() < expected_len {
        return Err(FnDsaError::InvalidInput {
            field: "secret_key",
            reason: "incomplete data",
        });
    }

    // Read f coefficients
    let mut f = Vec::with_capacity(n);
    for i in 0..n {
        f.push(bytes[6 + i] as i8);
    }

    // Read g coefficients
    let mut g = Vec::with_capacity(n);
    for i in 0..n {
        g.push(bytes[6 + n + i] as i8);
    }

    // Read F coefficients
    let mut big_f = Vec::with_capacity(n);
    let f_offset = 6 + n + n;
    for i in 0..n {
        let offset = f_offset + i * 2;
        let coeff = i16::from_le_bytes([bytes[offset], bytes[offset + 1]]);
        big_f.push(coeff);
    }

    // Read G coefficients
    let mut big_g = Vec::with_capacity(n);
    let g_offset = f_offset + n * 2;
    for i in 0..n {
        let offset = g_offset + i * 2;
        let coeff = i16::from_le_bytes([bytes[offset], bytes[offset + 1]]);
        big_g.push(coeff);
    }

    // Recompute Gram-Schmidt data
    let gs = compute_gram_schmidt(&f, &g, &big_f, &big_g);

    // Compute public key h = g * f^(-1) mod q
    let h = crate::ntru::compute_public_key(&f, &g, n)
        .map_err(|_| FnDsaError::InvalidKey { reason: "failed to compute h from f,g" })?;

    Ok(SecretKey {
        f,
        g,
        big_f,
        big_g,
        h,
        gs,
        params,
    })
}

/// Helper function to compute Gram-Schmidt from polynomials.
fn compute_gram_schmidt(f: &[i8], g: &[i8], big_f: &[i16], big_g: &[i16]) -> GramSchmidt {
    use crate::fft::fft;

    let n = f.len();

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

// ============================================================================
// Signature Encoding
// ============================================================================

/// Encodes a signature to bytes.
///
/// Format:
/// - 4 bytes: magic "FSG1"
/// - 2 bytes: n (little-endian)
/// - 40 bytes: nonce
/// - n * 2 bytes: s2 coefficients (little-endian i16)
pub fn encode_signature(sig: &Signature, n: usize) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(4 + 2 + NONCE_SIZE + n * 2);

    // Magic
    bytes.extend_from_slice(&SIG_MAGIC);

    // n
    bytes.extend_from_slice(&(n as u16).to_le_bytes());

    // Nonce
    bytes.extend_from_slice(&sig.nonce);

    // s2 coefficients
    for &coeff in &sig.s2 {
        bytes.extend_from_slice(&coeff.to_le_bytes());
    }

    bytes
}

/// Encodes a signature with Golomb-Rice compression.
///
/// Per FIPS 206, each coefficient is encoded as:
/// - Zero: 1 bit (the value 0)
/// - Nonzero: 1 (flag) + 1 (sign) + k (low bits) + unary(high)
///
/// Where `low = (|val| - 1) mod 2^k`, `high = (|val| - 1) / 2^k`,
/// and `unary(h) = h zero-bits followed by a 1-bit.
///
/// The parameter `k` (rice_k) controls the split: k = 8 for FALCON-512,
/// k = 9 for FALCON-1024.
pub fn encode_signature_compressed(sig: &Signature, n: usize) -> Vec<u8> {
    let k = rice_k_for_n(n);
    let k_mask = (1u64 << k) - 1;

    let mut bytes = Vec::with_capacity(4 + 2 + NONCE_SIZE + n);

    // Magic (different for compressed)
    bytes.extend_from_slice(b"FSC1");

    // n
    bytes.extend_from_slice(&(n as u16).to_le_bytes());

    // Nonce
    bytes.extend_from_slice(&sig.nonce);

    // Bit-packing state
    let mut bit_buffer: u64 = 0;
    let mut bit_count: usize = 0;

    for &coeff in &sig.s2 {
        if coeff == 0 {
            // Zero flag: single 0 bit
            // bit_buffer |= 0 << bit_count; (no-op)
            bit_count += 1;
        } else {
            // Nonzero flag: 1
            bit_buffer |= 1u64 << bit_count;
            bit_count += 1;

            // Sign bit: 0 = positive, 1 = negative
            let sign = if coeff < 0 { 1u64 } else { 0u64 };
            bit_buffer |= sign << bit_count;
            bit_count += 1;

            // Encode abs_val - 1
            let m = (coeff.unsigned_abs() - 1) as u64;
            let low = m & k_mask;
            let high = m >> k;

            // Low bits (k bits, LSB-first)
            bit_buffer |= low << bit_count;
            bit_count += k;

            // Unary: high zeros followed by a 1
            bit_count += high as usize; // high zero bits (buffer already 0)
            bit_buffer |= 1u64 << bit_count;
            bit_count += 1;
        }

        // Flush full bytes
        while bit_count >= 8 {
            bytes.push((bit_buffer & 0xFF) as u8);
            bit_buffer >>= 8;
            bit_count -= 8;
        }
    }

    // Flush remaining bits
    if bit_count > 0 {
        bytes.push((bit_buffer & 0xFF) as u8);
    }

    bytes
}

/// Returns the Golomb-Rice parameter k for a given polynomial degree n.
fn rice_k_for_n(n: usize) -> usize {
    match n {
        1024 => 9,
        _ => 8, // FALCON-512 and test params
    }
}

/// Decodes a signature from bytes.
pub fn decode_signature(bytes: &[u8]) -> Result<(Signature, usize)> {
    if bytes.len() < 6 {
        return Err(FnDsaError::InvalidInput {
            field: "signature",
            reason: "too short",
        });
    }

    // Check magic
    let magic = &bytes[0..4];
    let compressed = magic == b"FSC1";

    if magic != &SIG_MAGIC && !compressed {
        return Err(FnDsaError::InvalidInput {
            field: "signature",
            reason: "invalid magic bytes",
        });
    }

    // Read n
    let n = u16::from_le_bytes([bytes[4], bytes[5]]) as usize;

    if compressed {
        return decode_signature_compressed(bytes, n);
    }

    // Check length for uncompressed
    let expected_len = 6 + NONCE_SIZE + n * 2;
    if bytes.len() < expected_len {
        return Err(FnDsaError::InvalidInput {
            field: "signature",
            reason: "incomplete data",
        });
    }

    // Read nonce
    let mut nonce = [0u8; NONCE_SIZE];
    nonce.copy_from_slice(&bytes[6..6 + NONCE_SIZE]);

    // Read s2 coefficients
    let mut s2 = Vec::with_capacity(n);
    let s2_offset = 6 + NONCE_SIZE;
    for i in 0..n {
        let offset = s2_offset + i * 2;
        let coeff = i16::from_le_bytes([bytes[offset], bytes[offset + 1]]);
        s2.push(coeff);
    }

    // Validate coefficient range: s2 should be in centered representation [-Q/2, Q/2]
    let max_coeff = (crate::params::Q / 2) as i16;
    for &coeff in &s2 {
        if coeff < -max_coeff || coeff > max_coeff {
            return Err(FnDsaError::InvalidInput {
                field: "signature",
                reason: "s2 coefficient out of range",
            });
        }
    }

    Ok((Signature { nonce, s2 }, n))
}

/// Decodes a Golomb-Rice compressed signature.
fn decode_signature_compressed(bytes: &[u8], n: usize) -> Result<(Signature, usize)> {
    let k = rice_k_for_n(n);
    let k_mask = (1u64 << k) - 1;
    let header_len = 6 + NONCE_SIZE;

    if bytes.len() < header_len {
        return Err(FnDsaError::InvalidInput {
            field: "signature",
            reason: "incomplete compressed data",
        });
    }

    // Read nonce
    let mut nonce = [0u8; NONCE_SIZE];
    nonce.copy_from_slice(&bytes[6..6 + NONCE_SIZE]);

    // Decode s2 using Golomb-Rice bit unpacking
    let mut s2 = Vec::with_capacity(n);
    let data = &bytes[header_len..];

    let mut bit_buffer: u64 = 0;
    let mut bit_count: usize = 0;
    let mut byte_idx: usize = 0;

    // Helper: ensure at least `need` bits are buffered
    let refill = |buf: &mut u64, bc: &mut usize, bi: &mut usize, need: usize| {
        while *bc < need && *bi < data.len() {
            *buf |= (data[*bi] as u64) << *bc;
            *bc += 8;
            *bi += 1;
        }
    };

    let max_coeff = (crate::params::Q / 2) as i16;

    for _ in 0..n {
        // Read zero flag (1 bit)
        refill(&mut bit_buffer, &mut bit_count, &mut byte_idx, 1);
        if bit_count < 1 {
            return Err(FnDsaError::InvalidInput {
                field: "signature",
                reason: "truncated compressed data",
            });
        }
        let nonzero = bit_buffer & 1;
        bit_buffer >>= 1;
        bit_count -= 1;

        if nonzero == 0 {
            s2.push(0i16);
            continue;
        }

        // Read sign (1 bit) + low (k bits) = 1 + k bits
        refill(&mut bit_buffer, &mut bit_count, &mut byte_idx, 1 + k);
        if bit_count < 1 + k {
            return Err(FnDsaError::InvalidInput {
                field: "signature",
                reason: "truncated compressed data",
            });
        }
        let sign = bit_buffer & 1;
        bit_buffer >>= 1;
        bit_count -= 1;

        let low = bit_buffer & k_mask;
        bit_buffer >>= k;
        bit_count -= k;

        // Read unary: count zeros until a 1
        let mut high: u64 = 0;
        loop {
            refill(&mut bit_buffer, &mut bit_count, &mut byte_idx, 1);
            if bit_count < 1 {
                return Err(FnDsaError::InvalidInput {
                    field: "signature",
                    reason: "truncated unary in compressed data",
                });
            }
            let bit = bit_buffer & 1;
            bit_buffer >>= 1;
            bit_count -= 1;
            if bit == 1 {
                break;
            }
            high += 1;
            // Safety limit: high should never exceed a few bits for valid signatures
            if high > 64 {
                return Err(FnDsaError::InvalidInput {
                    field: "signature",
                    reason: "unary overflow in compressed data",
                });
            }
        }

        let abs_val = ((high << k) | low) as i16 + 1;
        if abs_val > max_coeff || abs_val < 0 {
            return Err(FnDsaError::InvalidInput {
                field: "signature",
                reason: "s2 coefficient out of range",
            });
        }
        let coeff = if sign == 1 { -abs_val } else { abs_val };
        s2.push(coeff);
    }

    Ok((Signature { nonce, s2 }, n))
}

// ============================================================================
// Key Pair Encoding
// ============================================================================

/// Encodes a full key pair (public + secret key).
pub fn encode_keypair(kp: &KeyPair) -> Vec<u8> {
    let pk_bytes = encode_public_key(&kp.pk);
    let sk_bytes = encode_secret_key(&kp.sk);

    let mut bytes = Vec::with_capacity(4 + pk_bytes.len() + sk_bytes.len());
    bytes.extend_from_slice(&(pk_bytes.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&pk_bytes);
    bytes.extend_from_slice(&sk_bytes);

    bytes
}

/// Decodes a full key pair.
pub fn decode_keypair(bytes: &[u8]) -> Result<KeyPair> {
    if bytes.len() < 4 {
        return Err(FnDsaError::InvalidInput {
            field: "keypair",
            reason: "too short",
        });
    }

    let pk_len = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as usize;

    if bytes.len() < 4 + pk_len {
        return Err(FnDsaError::InvalidInput {
            field: "keypair",
            reason: "incomplete public key",
        });
    }

    let pk = decode_public_key(&bytes[4..4 + pk_len])?;
    let sk = decode_secret_key(&bytes[4 + pk_len..])?;

    Ok(KeyPair { pk, sk })
}

// ============================================================================
// Hex encoding utilities
// ============================================================================

/// Encodes bytes to a hexadecimal string.
pub fn to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{:02x}", b)).collect()
}

/// Decodes a hexadecimal string to bytes.
pub fn from_hex(hex: &str) -> Result<Vec<u8>> {
    if hex.len() % 2 != 0 {
        return Err(FnDsaError::InvalidInput {
            field: "hex",
            reason: "odd length",
        });
    }

    let mut bytes = Vec::with_capacity(hex.len() / 2);
    for i in (0..hex.len()).step_by(2) {
        let byte = u8::from_str_radix(&hex[i..i + 2], 16).map_err(|_| {
            FnDsaError::InvalidInput {
                field: "hex",
                reason: "invalid hex character",
            }
        })?;
        bytes.push(byte);
    }

    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_public_key_roundtrip() {
        let pk = PublicKey {
            h: (0..512).map(|i| (i % 1000) as i16 - 500).collect(),
            params: FALCON_512,
        };

        let encoded = encode_public_key(&pk);
        let decoded = decode_public_key(&encoded).unwrap();

        assert_eq!(pk.h, decoded.h);
        assert_eq!(pk.params.n, decoded.params.n);
    }

    #[test]
    fn test_signature_roundtrip() {
        let sig = Signature {
            nonce: [42u8; NONCE_SIZE],
            s2: (0..512).map(|i| (i % 200) as i16 - 100).collect(),
        };

        let encoded = encode_signature(&sig, 512);
        let (decoded, n) = decode_signature(&encoded).unwrap();

        assert_eq!(n, 512);
        assert_eq!(sig.nonce, decoded.nonce);
        assert_eq!(sig.s2, decoded.s2);
    }

    #[test]
    fn test_signature_compressed_roundtrip() {
        let sig = Signature {
            nonce: [17u8; NONCE_SIZE],
            s2: (0..512).map(|i| ((i as i16 * 7) % 400) - 200).collect(),
        };

        let encoded = encode_signature_compressed(&sig, 512);
        let (decoded, n) = decode_signature(&encoded).unwrap();

        assert_eq!(n, 512);
        assert_eq!(sig.nonce, decoded.nonce);
        assert_eq!(sig.s2, decoded.s2);

        // Golomb-Rice compressed should be smaller than uncompressed (n*2 bytes for coeffs)
        let uncompressed = encode_signature(&sig, 512);
        assert!(
            encoded.len() < uncompressed.len(),
            "Compressed {} should be smaller than uncompressed {}",
            encoded.len(),
            uncompressed.len()
        );
    }

    #[test]
    fn test_signature_compressed_zeros() {
        // All-zero signature: maximally compressible (1 bit per coeff)
        let sig = Signature {
            nonce: [0u8; NONCE_SIZE],
            s2: vec![0i16; 512],
        };

        let encoded = encode_signature_compressed(&sig, 512);
        let (decoded, n) = decode_signature(&encoded).unwrap();

        assert_eq!(n, 512);
        assert_eq!(sig.s2, decoded.s2);

        // 512 coefficients * 1 bit each = 64 bytes for coeff data + header
        let header = 4 + 2 + NONCE_SIZE; // 46
        assert!(
            encoded.len() <= header + 65, // 64 bytes of bits + 1 rounding
            "All-zero compressed size {} should be near {} bytes",
            encoded.len(),
            header + 64
        );
    }

    #[test]
    fn test_signature_compressed_small_values() {
        // Small values typical of real FALCON signatures
        let sig = Signature {
            nonce: [42u8; NONCE_SIZE],
            s2: (0..512)
                .map(|i| {
                    let v = ((i as i16 * 13 + 7) % 30) - 15; // range [-15, 14]
                    v
                })
                .collect(),
        };

        let encoded = encode_signature_compressed(&sig, 512);
        let (decoded, _) = decode_signature(&encoded).unwrap();
        assert_eq!(sig.s2, decoded.s2);

        let uncompressed = encode_signature(&sig, 512);
        assert!(
            encoded.len() < uncompressed.len(),
            "Small-value compressed {} < uncompressed {}",
            encoded.len(),
            uncompressed.len()
        );
    }

    #[test]
    fn test_hex_roundtrip() {
        let original = vec![0xDE, 0xAD, 0xBE, 0xEF, 0x00, 0xFF];
        let hex = to_hex(&original);
        assert_eq!(hex, "deadbeef00ff");

        let decoded = from_hex(&hex).unwrap();
        assert_eq!(original, decoded);
    }

    #[test]
    fn test_invalid_magic() {
        let bad_bytes = b"XXXX\x00\x02";
        assert!(decode_public_key(bad_bytes).is_err());
        assert!(decode_secret_key(bad_bytes).is_err());
        assert!(decode_signature(bad_bytes).is_err());
    }
}
