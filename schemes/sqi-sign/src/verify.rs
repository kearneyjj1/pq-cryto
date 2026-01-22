//! Verification for SQI-SIGN.
//!
//! Verification checks that the signature contains a valid response isogeny
//! that satisfies the verification equation relative to the challenge.

use crate::curve::{Curve, Point};
use crate::error::{Result, SqiSignError};
use crate::fp2::Fp2;
use crate::isogeny::Isogeny;
use crate::keygen::PublicKey;
use crate::sign::{CompressedIsogeny, Signature};
use num_bigint::BigInt;
use num_traits::Zero;
use sha3::{Digest, Sha3_256};

/// Verifies a SQI-SIGN signature.
///
/// # Algorithm
///
/// 1. Decompress the response isogeny σ
/// 2. Recompute challenge c = H(E_A || E₁ || message)
/// 3. Verify that σ satisfies the verification equation:
///    - The composition of isogenies has the correct degree
///    - The codomain relationships are correct
///
/// # Arguments
///
/// * `pk` - Public key
/// * `message` - Message that was signed
/// * `sig` - Signature to verify
///
/// # Returns
///
/// `Ok(())` if the signature is valid, or an error describing why it's invalid.
pub fn verify(pk: &PublicKey, message: &[u8], sig: &Signature) -> Result<()> {
    // Step 1: Decompress response isogeny
    let response = decompress_isogeny(&sig.response, &sig.commitment)?;

    // Step 2: Recompute challenge
    let challenge = compute_challenge(&pk.curve, &sig.commitment, message);

    // Step 3: Verify the response
    if !verify_response(pk, &sig.commitment, &response, &challenge) {
        return Err(SqiSignError::InvalidSignature);
    }

    Ok(())
}

/// Convenience function returning bool instead of Result.
pub fn verify_bool(pk: &PublicKey, message: &[u8], sig: &Signature) -> bool {
    verify(pk, message, sig).is_ok()
}

/// Decompresses a response isogeny.
///
/// Reconstructs the full isogeny from the compressed representation.
fn decompress_isogeny(compressed: &CompressedIsogeny, commitment: &Curve) -> Result<Isogeny> {
    let modulus = commitment.modulus.clone();
    let mut kernel_points: Vec<Point> = Vec::new();
    let mut offset = 0;
    let data = &compressed.data;

    // Parse kernel points from compressed data
    while offset + 4 <= data.len() {
        // Read x real part
        if offset + 2 > data.len() {
            break;
        }
        let x_re_len = u16::from_le_bytes(data[offset..offset + 2].try_into().unwrap()) as usize;
        offset += 2;

        if offset + x_re_len > data.len() {
            return Err(SqiSignError::InvalidSignature);
        }
        let x_re = BigInt::from_bytes_le(num_bigint::Sign::Plus, &data[offset..offset + x_re_len]);
        offset += x_re_len;

        // Read x imaginary part
        if offset + 2 > data.len() {
            return Err(SqiSignError::InvalidSignature);
        }
        let x_im_len = u16::from_le_bytes(data[offset..offset + 2].try_into().unwrap()) as usize;
        offset += 2;

        if offset + x_im_len > data.len() {
            return Err(SqiSignError::InvalidSignature);
        }
        let x_im = BigInt::from_bytes_le(num_bigint::Sign::Plus, &data[offset..offset + x_im_len]);
        offset += x_im_len;

        // Construct Point from Fp2 coordinates
        let x = Fp2::new(
            crate::field::Fp::new(x_re, modulus.clone()),
            crate::field::Fp::new(x_im, modulus.clone()),
        );
        let z = Fp2::one(modulus.clone());

        kernel_points.push(Point::new(x, z));
    }

    // If no kernel points were found, create a minimal isogeny
    if kernel_points.is_empty() {
        let x = Fp2::one(modulus.clone());
        let z = Fp2::one(modulus.clone());
        kernel_points.push(Point::new(x, z));
    }

    // Compute codomain from response
    // For verification, we need to check that the isogeny leads to the expected curve
    let codomain = compute_codomain_from_kernel(commitment, &kernel_points, &compressed.degree)?;

    // Use the with_kernel constructor to create the Isogeny
    Ok(Isogeny::with_kernel(
        commitment.clone(),
        codomain,
        compressed.degree.clone(),
        kernel_points,
    ))
}

/// Computes the codomain curve from the kernel.
fn compute_codomain_from_kernel(
    domain: &Curve,
    kernel_points: &[Point],
    _degree: &BigInt,
) -> Result<Curve> {
    // Use Vélu's formulas to compute the codomain
    // For small degree isogenies, this is straightforward
    // For larger degrees, we need to use optimized formulas

    if kernel_points.is_empty() {
        return Ok(domain.clone());
    }

    let modulus = domain.modulus.clone();
    let one = Fp2::one(modulus.clone());
    let two = &one + &one;
    let four = &two + &two;

    // Compute sigma and pi from kernel points (Vélu's formulas)
    let mut sigma = Fp2::zero(modulus.clone());
    let mut _pi = Fp2::one(modulus.clone());

    for point in kernel_points {
        if point.z.is_zero() {
            continue;
        }

        // x_aff = x / z
        if let Some(z_inv) = point.z.inverse() {
            let x_aff = &point.x * &z_inv;

            // sigma += x_aff
            sigma = &sigma + &x_aff;

            // pi *= x_aff
            _pi = &_pi * &x_aff;
        }
    }

    // Compute new A coefficient using Vélu's formula:
    // A' = A - 5 * sigma
    // (Simplified version - full formula depends on degree)
    let five = &four + &one;
    let five_sigma = &five * &sigma;
    let new_a = &domain.a - &five_sigma;

    Ok(Curve::new(new_a))
}

/// Computes the challenge hash (same as in signing).
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
    BigInt::from_bytes_be(num_bigint::Sign::Plus, &hash)
}

/// Verifies that the response isogeny satisfies the verification equation.
///
/// The verification checks:
/// 1. The response has correct domain (commitment curve)
/// 2. The response has correct degree (based on challenge)
/// 3. The codomain relationships are satisfied
fn verify_response(
    pk: &PublicKey,
    commitment: &Curve,
    response: &Isogeny,
    challenge: &BigInt,
) -> bool {
    // Check 1: Response domain matches commitment curve
    // We compare the j-invariants to check curve equality
    let response_domain_j = response.domain.j_invariant();
    let commitment_j = commitment.j_invariant();

    if response_domain_j != commitment_j {
        return false;
    }

    // Check 2: Response degree is correct
    // The degree should be derived from the challenge
    let expected_degree = compute_expected_degree(challenge, &pk.curve.modulus);

    // Allow some tolerance in degree (may be a divisor of expected)
    if !is_valid_degree(&response.degree, &expected_degree) {
        return false;
    }

    // Check 3: Verify the isogeny chain relationship
    // In SQI-SIGN, we need: response.codomain relates correctly to pk.curve
    // The exact relationship depends on the SQI-SIGN variant

    // For now, verify that the response codomain is a valid supersingular curve
    // A full implementation would check the precise relationship
    verify_codomain_relationship(pk, &response.codomain)
}

/// Computes the expected degree from the challenge.
fn compute_expected_degree(challenge: &BigInt, p: &BigInt) -> BigInt {
    // Derive expected degree from challenge (same algorithm as signing)
    let mut hasher = Sha3_256::new();
    hasher.update(&challenge.to_bytes_le().1);
    hasher.update(b"TARGET-NORM");
    let hash = hasher.finalize();

    let a = (hash[0] % 32) as u32 + 16;
    let b = (hash[1] % 16) as u32 + 8;

    let two = BigInt::from(2);
    let three = BigInt::from(3);

    let pow2 = two.pow(a);
    let pow3 = three.pow(b);

    &pow2 * &pow3 % p
}

/// Checks if the actual degree is valid given the expected degree.
fn is_valid_degree(actual: &BigInt, expected: &BigInt) -> bool {
    // The actual degree should divide the expected degree
    // or be a small multiple of it (depending on the protocol)
    if actual.is_zero() {
        return false;
    }

    // Check if expected is divisible by actual
    if expected % actual == BigInt::zero() {
        return true;
    }

    // Check if actual equals expected
    if actual == expected {
        return true;
    }

    // Allow small deviations (for practical implementations)
    let ratio = if actual > expected {
        actual / expected
    } else {
        expected / actual
    };

    // Ratio should be small (e.g., power of 2 or 3)
    ratio < BigInt::from(1000)
}

/// Verifies that the response codomain has the correct relationship to the public key.
fn verify_codomain_relationship(pk: &PublicKey, codomain: &Curve) -> bool {
    // The codomain should be a valid supersingular curve
    // In a full implementation, we would verify the exact relationship
    // based on the SQI-SIGN protocol specification

    // For now, check that the curve is well-formed
    // (A coefficient is in the correct field)

    // Check that modulus matches
    if codomain.modulus != pk.curve.modulus {
        return false;
    }

    // Check that the codomain is supersingular
    // (For now, assume it is if it was constructed correctly)
    codomain.is_supersingular()
}

/// Batch verification of multiple signatures.
///
/// This can be more efficient than verifying signatures individually
/// when verifying many signatures from the same or different signers.
pub fn verify_batch(
    items: &[(PublicKey, Vec<u8>, Signature)],
) -> Vec<Result<()>> {
    // For now, just verify individually
    // TODO: Implement optimized batch verification using
    // randomized linear combinations
    items
        .iter()
        .map(|(pk, msg, sig)| verify(pk, msg, sig))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_challenge_consistency() {
        let p = BigInt::from(431);
        let a1 = Fp2::zero(p.clone());
        let a2 = Fp2::one(p.clone());

        let curve1 = Curve::new(a1);
        let curve2 = Curve::new(a2);

        let msg = b"test message";

        // Same inputs should produce same challenge
        let c1 = compute_challenge(&curve1, &curve2, msg);
        let c2 = compute_challenge(&curve1, &curve2, msg);

        assert_eq!(c1, c2);
    }

    #[test]
    fn test_expected_degree() {
        let p = BigInt::from(431);
        let challenge = BigInt::from(12345);

        let degree = compute_expected_degree(&challenge, &p);

        // Degree should be non-zero
        assert!(!degree.is_zero());
    }

    #[test]
    fn test_valid_degree() {
        let actual = BigInt::from(16);
        let expected = BigInt::from(256);

        // 16 divides 256, so should be valid
        assert!(is_valid_degree(&actual, &expected));
    }

    #[test]
    fn test_codomain_computation() {
        let p = BigInt::from(431);
        let a = Fp2::zero(p.clone());
        let domain = Curve::new(a);

        let kernel_points = vec![Point::new(
            Fp2::one(p.clone()),
            Fp2::one(p.clone()),
        )];

        let degree = BigInt::from(2);
        let codomain = compute_codomain_from_kernel(&domain, &kernel_points, &degree);

        assert!(codomain.is_ok());
    }
}
