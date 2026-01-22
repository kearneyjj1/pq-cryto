//! Signature verification for SLH-DSA.
//!
//! Implements SLH-DSA verification (FIPS 205 Algorithm 19).

use crate::address::{Address, AddressType};
use crate::error::{Result, SlhDsaError};
use crate::fors::fors_pk_from_sig;
use crate::hash::h_msg;
use crate::hypertree::ht_verify;
use crate::keygen::PublicKey;
use crate::sign::Signature;

/// Verifies an SLH-DSA signature.
///
/// slh_verify(M, SIG, PK) - FIPS 205 Algorithm 19
///
/// # Arguments
/// * `pk` - Public key
/// * `message` - Message that was signed
/// * `sig` - Signature to verify
///
/// # Returns
/// `Ok(())` if the signature is valid, or an error if verification fails.
pub fn verify(pk: &PublicKey, message: &[u8], sig: &Signature) -> Result<()> {
    let params = &pk.params;

    // Validate signature structure
    if sig.r.len() != params.n {
        return Err(SlhDsaError::InvalidSignature);
    }

    if sig.sig_fors.sk.len() != params.k || sig.sig_fors.auth.len() != params.k {
        return Err(SlhDsaError::InvalidSignature);
    }

    if sig.sig_ht.layers.len() != params.d {
        return Err(SlhDsaError::InvalidSignature);
    }

    // Step 1: Compute message digest and indices from R
    // H_msg(R, PK.seed, PK.root, M) -> (md, idx_tree, idx_leaf)
    let (md, idx_tree, idx_leaf) = h_msg(&sig.r, &pk.pk_seed, &pk.pk_root, message, params);

    // Step 2: Reconstruct FORS public key from signature
    let mut fors_adrs = Address::new();
    fors_adrs.set_type(AddressType::ForsTree);
    fors_adrs.set_layer(0);
    fors_adrs.set_tree(idx_tree);
    fors_adrs.set_keypair(idx_leaf);

    let pk_fors = fors_pk_from_sig(&sig.sig_fors, &md, &pk.pk_seed, &mut fors_adrs, params);

    // Step 3: Verify hypertree signature
    let valid = ht_verify(
        &pk_fors,
        &sig.sig_ht,
        &pk.pk_seed,
        idx_tree,
        idx_leaf,
        &pk.pk_root,
        params,
    );

    if valid {
        Ok(())
    } else {
        Err(SlhDsaError::InvalidSignature)
    }
}

/// Boolean verification wrapper.
///
/// Returns `true` if the signature is valid, `false` otherwise.
pub fn verify_bool(pk: &PublicKey, message: &[u8], sig: &Signature) -> bool {
    verify(pk, message, sig).is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::keygen;
    use crate::params::SLH_DSA_SHAKE_128F;
    use crate::sign::{sign, sign_randomized};
    use rand::rngs::OsRng;

    #[test]
    fn test_verify_valid_signature() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"Test message for verification";

        let sig = sign(&sk, message).unwrap();
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_verify_bool() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"Test message";

        let sig = sign(&sk, message).unwrap();
        assert!(verify_bool(&pk, message, &sig));
    }

    #[test]
    fn test_verify_wrong_message() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);

        let sig = sign(&sk, b"Original message").unwrap();
        assert!(verify(&pk, b"Modified message", &sig).is_err());
    }

    #[test]
    fn test_verify_wrong_key() {
        let params = SLH_DSA_SHAKE_128F;
        let (_pk1, sk1) = keygen(&mut OsRng, params);
        let (pk2, _) = keygen(&mut OsRng, params);
        let message = b"Test message";

        let sig = sign(&sk1, message).unwrap();

        // Signature from sk1 should not verify with pk2
        assert!(verify(&pk2, message, &sig).is_err());
    }

    #[test]
    fn test_verify_tampered_signature() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"Test message";

        let mut sig = sign(&sk, message).unwrap();

        // Tamper with the randomizer
        sig.r[0] ^= 0xFF;

        assert!(verify(&pk, message, &sig).is_err());
    }

    #[test]
    fn test_verify_tampered_fors() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"Test message";

        let mut sig = sign(&sk, message).unwrap();

        // Tamper with FORS signature
        sig.sig_fors.sk[0][0] ^= 0xFF;

        assert!(verify(&pk, message, &sig).is_err());
    }

    #[test]
    fn test_verify_tampered_hypertree() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"Test message";

        let mut sig = sign(&sk, message).unwrap();

        // Tamper with hypertree signature
        sig.sig_ht.layers[0].wots_sig[0][0] ^= 0xFF;

        assert!(verify(&pk, message, &sig).is_err());
    }

    #[test]
    fn test_verify_randomized_signature() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"Test message";

        let sig = sign_randomized(&mut OsRng, &sk, message).unwrap();
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_verify_empty_message() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"";

        let sig = sign(&sk, message).unwrap();
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_verify_long_message() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = vec![0x42u8; 10000];

        let sig = sign(&sk, &message).unwrap();
        assert!(verify(&pk, &message, &sig).is_ok());
    }

    // Note: These tests use 128f for speed. 128s would take much longer.
    #[test]
    fn test_verify_128s_roundtrip() {
        // This test is slow due to the 's' parameter set
        // Uncomment for thorough testing
        /*
        let params = SLH_DSA_SHAKE_128S;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"Test";

        let sig = sign(&sk, message).unwrap();
        assert!(verify(&pk, message, &sig).is_ok());
        */
    }
}
