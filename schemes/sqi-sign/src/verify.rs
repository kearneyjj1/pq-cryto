//! Verification for SQI-SIGN.
//!
//! Verification checks that the signature contains a valid response isogeny
//! that satisfies the verification equation relative to the challenge.

use crate::curve::Curve;
use crate::error::{Result, SqiSignError};
use crate::isogeny::Isogeny;
use crate::keygen::PublicKey;
use crate::sign::{CompressedIsogeny, Signature};
use num_bigint::BigInt;
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
fn decompress_isogeny(compressed: &CompressedIsogeny, commitment: &Curve) -> Result<Isogeny> {
    // TODO: Implement isogeny decompression
    unimplemented!("Isogeny decompression not yet implemented")
}

/// Computes the challenge hash (same as in signing).
fn compute_challenge(public_curve: &Curve, commitment: &Curve, message: &[u8]) -> BigInt {
    let mut hasher = Sha3_256::new();

    // TODO: Properly serialize curves
    // hasher.update(&public_curve.to_bytes());
    // hasher.update(&commitment.to_bytes());
    hasher.update(message);

    let hash = hasher.finalize();
    BigInt::from_bytes_be(num_bigint::Sign::Plus, &hash)
}

/// Verifies that the response isogeny satisfies the verification equation.
fn verify_response(
    pk: &PublicKey,
    commitment: &Curve,
    response: &Isogeny,
    challenge: &BigInt,
) -> bool {
    // TODO: Implement verification checks:
    //
    // 1. Check that response has correct domain (commitment curve)
    // 2. Check that response has correct degree (based on challenge)
    // 3. Check the codomain relationships between pk.curve, commitment, response.codomain
    // 4. Verify any auxiliary constraints required by the scheme

    unimplemented!("Response verification not yet implemented")
}

/// Batch verification of multiple signatures.
///
/// This can be more efficient than verifying signatures individually
/// when verifying many signatures from the same or different signers.
pub fn verify_batch(
    items: &[(PublicKey, Vec<u8>, Signature)],
) -> Vec<Result<()>> {
    // For now, just verify individually
    // TODO: Implement optimized batch verification
    items
        .iter()
        .map(|(pk, msg, sig)| verify(pk, msg, sig))
        .collect()
}

#[cfg(test)]
mod tests {
    // Tests will be added as functions are implemented
}
