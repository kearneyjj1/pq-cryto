//! Signing for SQI-SIGN.
//!
//! SQI-SIGN uses a commitment-challenge-response structure similar to
//! Fiat-Shamir identification schemes.
//!
//! # Signing Algorithm
//!
//! 1. **Commitment**: Sample random isogeny ψ: E₀ → E₁, publish E₁
//! 2. **Challenge**: Hash (E_A, E₁, message) to get challenge scalar c
//! 3. **Response**: Compute isogeny σ: E₁ → E₂ such that φ ∘ σ has
//!    specific properties related to c
//!
//! The response uses KLPT and quaternion algebra computations to find
//! a "short" isogeny that satisfies the verification equation.

use crate::curve::Curve;
use crate::error::Result;
use crate::fp2::Fp2;
use crate::ideal::{ideal_to_isogeny, isogeny_to_ideal, klpt, LeftIdeal};
use crate::isogeny::Isogeny;
use crate::keygen::SecretKey;
use crate::params::Params;
use crate::quaternion::{MaximalOrder, Quaternion, Rational};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use rand::{CryptoRng, RngCore};
use sha3::{Digest, Sha3_256};

/// Maximum signing attempts before giving up.
const MAX_ATTEMPTS: u32 = 100;

/// Converts a BigInt to a fixed-size byte array in little-endian format.
///
/// If the value requires fewer bytes, the array is zero-padded.
/// If the value requires more bytes, it is truncated (caller must ensure this doesn't happen).
fn bigint_to_fixed_bytes(value: &BigInt, size: usize) -> Vec<u8> {
    let (_, mut bytes) = value.to_bytes_le();

    // Pad with zeros if needed
    if bytes.len() < size {
        bytes.resize(size, 0);
    }

    // Truncate if too long (shouldn't happen if size is chosen correctly)
    bytes.truncate(size);

    bytes
}

/// Derives a field element from a seed with full entropy.
///
/// Uses SHA3-256 with domain separation to derive a coefficient.
/// The full 256-bit hash output is used, then reduced modulo p.
/// For primes larger than 256 bits, multiple rounds are used.
fn derive_commitment_coefficient(seed: &[u8], domain: &[u8], modulus: &BigInt) -> BigInt {
    // Calculate how many bytes we need (prime size + security margin)
    let prime_bits = modulus.bits();
    let bytes_needed = ((prime_bits + 128) / 8) as usize; // Extra 128 bits for uniform distribution

    let mut output = Vec::with_capacity(bytes_needed);
    let mut counter: u32 = 0;

    // Use counter mode to generate enough bytes
    while output.len() < bytes_needed {
        let mut hasher = Sha3_256::new();
        hasher.update(seed);
        hasher.update(domain);
        hasher.update(&counter.to_le_bytes());
        let hash = hasher.finalize();
        output.extend_from_slice(&hash);
        counter += 1;
    }

    output.truncate(bytes_needed);

    // Convert to BigInt and reduce modulo p
    BigInt::from_bytes_le(num_bigint::Sign::Plus, &output) % modulus
}

/// A SQI-SIGN signature.
#[derive(Clone, Debug)]
pub struct Signature {
    /// The commitment curve E₁.
    pub commitment: Curve,
    /// The response isogeny (compressed representation).
    pub response: CompressedIsogeny,
}

/// Compressed representation of a response isogeny.
///
/// The full isogeny can be reconstructed during verification.
#[derive(Clone, Debug)]
pub struct CompressedIsogeny {
    /// Compressed data for the isogeny (degree and auxiliary data).
    pub data: Vec<u8>,
    /// Degree of the isogeny.
    pub degree: BigInt,
}

impl Signature {
    /// Serializes the signature to bytes using canonical encoding.
    ///
    /// # Encoding Format
    ///
    /// The signature uses fixed-size encoding to ensure canonical representation:
    /// - Commitment A coefficient (real part): `field_bytes` bytes, little-endian
    /// - Commitment A coefficient (imag part): `field_bytes` bytes, little-endian
    /// - Response degree: 8 bytes, little-endian
    /// - Response data length: 4 bytes, little-endian
    /// - Response data: variable length
    ///
    /// All field elements are reduced to canonical form [0, p-1] before encoding.
    pub fn to_bytes(&self) -> Vec<u8> {
        let modulus = &self.commitment.modulus;
        let field_bytes = (modulus.bits() as usize + 7) / 8;

        let mut bytes = Vec::new();

        // Serialize commitment curve A coefficient with fixed-size encoding
        // Ensure canonical form (value in [0, p-1])
        let re_canonical = &self.commitment.a.re.value % modulus;
        let im_canonical = &self.commitment.a.im.value % modulus;

        // Convert to fixed-size byte arrays
        let re_bytes = bigint_to_fixed_bytes(&re_canonical, field_bytes);
        let im_bytes = bigint_to_fixed_bytes(&im_canonical, field_bytes);

        bytes.extend(&re_bytes);
        bytes.extend(&im_bytes);

        // Serialize response degree (fixed 8 bytes for up to 64-bit degrees)
        let degree_bytes = bigint_to_fixed_bytes(&self.response.degree, 8);
        bytes.extend(&degree_bytes);

        // Serialize response data with length prefix
        bytes.extend(&(self.response.data.len() as u32).to_le_bytes());
        bytes.extend(&self.response.data);

        bytes
    }

    /// Serializes the signature using the legacy variable-length format.
    /// Provided for backwards compatibility only.
    #[deprecated(note = "Use to_bytes() for canonical encoding")]
    pub fn to_bytes_legacy(&self) -> Vec<u8> {
        let mut bytes = Vec::new();

        // Serialize commitment curve (A coefficient)
        let re_bytes = self.commitment.a.re.value.to_bytes_le().1;
        let im_bytes = self.commitment.a.im.value.to_bytes_le().1;

        // Length prefix for re
        bytes.extend(&(re_bytes.len() as u32).to_le_bytes());
        bytes.extend(&re_bytes);

        // Length prefix for im
        bytes.extend(&(im_bytes.len() as u32).to_le_bytes());
        bytes.extend(&im_bytes);

        // Serialize response
        let degree_bytes = self.response.degree.to_bytes_le().1;
        bytes.extend(&(degree_bytes.len() as u32).to_le_bytes());
        bytes.extend(&degree_bytes);

        bytes.extend(&(self.response.data.len() as u32).to_le_bytes());
        bytes.extend(&self.response.data);

        bytes
    }

    /// Deserializes a signature from canonical encoding.
    ///
    /// Validates that field elements are in canonical form [0, p-1].
    /// Rejects signatures with non-canonical representations to prevent malleability.
    pub fn from_bytes(bytes: &[u8], modulus: &BigInt) -> Result<Self> {
        let field_bytes = (modulus.bits() as usize + 7) / 8;
        let mut offset = 0;

        // Read re (fixed-size)
        if bytes.len() < offset + field_bytes {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }
        let re = BigInt::from_bytes_le(num_bigint::Sign::Plus, &bytes[offset..offset + field_bytes]);
        offset += field_bytes;

        // Validate canonical form: re must be < modulus
        if &re >= modulus {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }

        // Read im (fixed-size)
        if bytes.len() < offset + field_bytes {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }
        let im = BigInt::from_bytes_le(num_bigint::Sign::Plus, &bytes[offset..offset + field_bytes]);
        offset += field_bytes;

        // Validate canonical form: im must be < modulus
        if &im >= modulus {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }

        // Construct commitment curve
        let a = Fp2::new(
            crate::field::Fp::new(re, modulus.clone()),
            crate::field::Fp::new(im, modulus.clone()),
        );
        let commitment = Curve::new(a);

        // Read degree (fixed 8 bytes)
        if bytes.len() < offset + 8 {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }
        let degree = BigInt::from_bytes_le(num_bigint::Sign::Plus, &bytes[offset..offset + 8]);
        offset += 8;

        // Validate degree is positive
        if degree <= BigInt::zero() {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }

        // Read data length (4 bytes)
        if bytes.len() < offset + 4 {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }
        let data_len = u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap()) as usize;
        offset += 4;

        // Validate data_len is reasonable (prevent DoS via large allocations)
        if data_len > 1024 * 1024 {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }

        if bytes.len() < offset + data_len {
            return Err(crate::error::SqiSignError::InvalidSignature);
        }
        let data = bytes[offset..offset + data_len].to_vec();

        Ok(Self {
            commitment,
            response: CompressedIsogeny { data, degree },
        })
    }

    /// Returns the size of the signature in bytes.
    pub fn size(&self) -> usize {
        self.to_bytes().len()
    }
}

/// Signs a message using SQI-SIGN.
///
/// # Algorithm
///
/// 1. Sample random commitment isogeny ψ: E₀ → E₁
/// 2. Compute challenge c = H(E_A || E₁ || message)
/// 3. Use KLPT to find ideal I_σ corresponding to response
/// 4. Translate I_σ to response isogeny σ
/// 5. Compress σ for the signature
///
/// # Arguments
///
/// * `rng` - Cryptographically secure random number generator
/// * `sk` - Secret key
/// * `message` - Message to sign
///
/// # Returns
///
/// A signature on success, or an error if signing fails.
pub fn sign<R: RngCore + CryptoRng>(
    rng: &mut R,
    sk: &SecretKey,
    message: &[u8],
) -> Result<Signature> {
    for _attempt in 0..MAX_ATTEMPTS {
        // Step 1: Generate random commitment
        let (commitment_curve, commitment_isogeny) = sample_commitment(rng, &sk.params)?;

        // Step 2: Compute challenge
        let challenge = compute_challenge(&sk.public_key().curve, &commitment_curve, message);

        // Step 3: Compute response using KLPT
        if let Some(response) = compute_response(sk, &commitment_isogeny, &challenge) {
            // Step 4: Compress response
            let compressed = compress_isogeny(&response);

            return Ok(Signature {
                commitment: commitment_curve,
                response: compressed,
            });
        }
        // If response computation failed, try again with new commitment
    }

    Err(crate::error::SqiSignError::SigningFailed {
        attempts: MAX_ATTEMPTS,
    })
}

/// Samples a random commitment isogeny ψ: E₀ → E₁.
///
/// The commitment is a random isogeny of smooth degree from E₀.
fn sample_commitment<R: RngCore + CryptoRng>(
    rng: &mut R,
    params: &Params,
) -> Result<(Curve, Isogeny)> {
    let p = params.prime().clone();
    let o0 = MaximalOrder::standard(p.clone());

    // Generate random seed with full entropy
    let mut seed = [0u8; 32];
    rng.fill_bytes(&mut seed);

    // Derive coefficients with full entropy using separate domain separators
    let a = derive_commitment_coefficient(&seed, b"SQI-SIGN-COMMIT-A", &p) + BigInt::one();
    let b = derive_commitment_coefficient(&seed, b"SQI-SIGN-COMMIT-B", &p);
    let c = derive_commitment_coefficient(&seed, b"SQI-SIGN-COMMIT-C", &p);
    let d = derive_commitment_coefficient(&seed, b"SQI-SIGN-COMMIT-D", &p);

    // Create generator quaternion with all four components for better security
    let generator = Quaternion::new(
        Rational::from_int(a),
        Rational::from_int(b),
        Rational::from_int(c),
        Rational::from_int(d),
        p.clone(),
    );

    // Create commitment ideal
    let commitment_ideal = LeftIdeal::principal(&o0, generator);

    // Translate ideal to isogeny
    let commitment_isogeny = ideal_to_isogeny(&commitment_ideal);
    let commitment_curve = commitment_isogeny.codomain.clone();

    Ok((commitment_curve, commitment_isogeny))
}

/// Computes the challenge hash.
///
/// c = H(E_A || E₁ || message) where:
/// - E_A is the public key curve
/// - E₁ is the commitment curve
/// - message is the message being signed
fn compute_challenge(public_curve: &Curve, commitment: &Curve, message: &[u8]) -> BigInt {
    let mut hasher = Sha3_256::new();

    // Serialize public curve A coefficient
    let pk_re = public_curve.a.re.value.to_bytes_le().1;
    let pk_im = public_curve.a.im.value.to_bytes_le().1;
    hasher.update(&pk_re);
    hasher.update(&pk_im);

    // Serialize commitment curve A coefficient
    let cm_re = commitment.a.re.value.to_bytes_le().1;
    let cm_im = commitment.a.im.value.to_bytes_le().1;
    hasher.update(&cm_re);
    hasher.update(&cm_im);

    // Hash the message
    hasher.update(message);

    let hash = hasher.finalize();

    // Convert hash to challenge scalar
    BigInt::from_bytes_be(num_bigint::Sign::Plus, &hash)
}

/// Computes the response isogeny using KLPT.
///
/// The response computation finds an isogeny σ from the commitment curve E₁
/// such that the verification equation is satisfied.
///
/// Specifically, we need to find an ideal I_σ such that:
/// - I_σ has the correct norm (determined by challenge)
/// - I_σ corresponds to an isogeny starting at E₁
fn compute_response(
    sk: &SecretKey,
    commitment: &Isogeny,
    challenge: &BigInt,
) -> Option<Isogeny> {
    let p = &sk.endo_ring.p;

    // Convert commitment isogeny to ideal
    let commitment_ideal = isogeny_to_ideal(commitment, &sk.endo_ring);

    // Compute target norm based on challenge
    // The target should be smooth (product of small primes) for efficient evaluation
    let target_norm = compute_target_norm(challenge, p);

    // Use KLPT to find connecting element (returns a quaternion)
    let response_element = klpt(&commitment_ideal, &target_norm)?;

    // Create the response ideal from the KLPT element
    let response_ideal = LeftIdeal::principal(&sk.endo_ring, response_element);

    // Translate response ideal to isogeny
    let response = ideal_to_isogeny(&response_ideal);

    Some(response)
}

/// Computes the target norm for the response isogeny.
///
/// The target norm is derived from the challenge and must be smooth.
fn compute_target_norm(challenge: &BigInt, p: &BigInt) -> BigInt {
    // Use challenge to derive a smooth norm
    // Target norm should be of form 2^a * 3^b for efficient evaluation
    let mut hasher = Sha3_256::new();
    hasher.update(&challenge.to_bytes_le().1);
    hasher.update(b"TARGET-NORM");
    let hash = hasher.finalize();

    // Extract exponents from hash
    let a = (hash[0] % 32) as u32 + 16; // 2^a where 16 <= a <= 47
    let b = (hash[1] % 16) as u32 + 8;  // 3^b where 8 <= b <= 23

    // Compute target norm: 2^a * 3^b (mod p if needed)
    let two = BigInt::from(2);
    let three = BigInt::from(3);

    let pow2 = two.pow(a);
    let pow3 = three.pow(b);

    &pow2 * &pow3 % p
}

/// Compresses an isogeny for signature transmission.
///
/// The compressed form includes:
/// - Degree of the isogeny
/// - Auxiliary data for decompression
fn compress_isogeny(isogeny: &Isogeny) -> CompressedIsogeny {
    let mut data = Vec::new();

    // Store kernel polynomial coefficients (x-coordinates of kernel points)
    for point in isogeny.kernel_points() {
        if let Some(norm) = point.normalize() {
            let x_re = norm.x.re.value.to_bytes_le().1;
            let x_im = norm.x.im.value.to_bytes_le().1;

            data.extend(&(x_re.len() as u16).to_le_bytes());
            data.extend(&x_re);
            data.extend(&(x_im.len() as u16).to_le_bytes());
            data.extend(&x_im);
        }
    }

    CompressedIsogeny {
        data,
        degree: isogeny.degree.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;

    #[test]
    fn test_compute_challenge() {
        let p = BigInt::from(431);
        let a1 = Fp2::zero(p.clone());
        let a2 = Fp2::one(p.clone());

        let curve1 = Curve::new(a1);
        let curve2 = Curve::new(a2);

        let msg1 = b"test message 1";
        let msg2 = b"test message 2";

        let c1 = compute_challenge(&curve1, &curve2, msg1);
        let c2 = compute_challenge(&curve1, &curve2, msg2);

        // Different messages should produce different challenges
        assert_ne!(c1, c2);
    }

    #[test]
    fn test_signature_serialization() {
        let p = BigInt::from(431);
        let a = Fp2::zero(p.clone());
        let curve = Curve::new(a);

        let sig = Signature {
            commitment: curve,
            response: CompressedIsogeny {
                data: vec![1, 2, 3, 4],
                degree: BigInt::from(16),
            },
        };

        let bytes = sig.to_bytes();
        let sig2 = Signature::from_bytes(&bytes, &p).unwrap();

        assert_eq!(sig.response.degree, sig2.response.degree);
        assert_eq!(sig.response.data, sig2.response.data);
    }

    #[test]
    fn test_compute_target_norm() {
        let p = BigInt::from(431);
        let challenge = BigInt::from(12345);

        let target = compute_target_norm(&challenge, &p);

        // Target should be non-zero and less than p
        assert!(!target.is_zero());
    }
}
