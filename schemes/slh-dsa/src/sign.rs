//! Signature generation for SLH-DSA.
//!
//! Implements SLH-DSA signing (FIPS 205 Algorithm 18).
//!
//! # Signature Structure
//!
//! An SLH-DSA signature consists of:
//! - R: Randomizer (n bytes)
//! - SIG_FORS: FORS signature (k * (a + 1) * n bytes)
//! - SIG_HT: Hypertree signature ((h + d * wots_len) * n bytes)

use crate::address::{Address, AddressType};
use crate::error::Result;
use crate::fors::{fors_pk_from_sig, fors_sign, ForsSig};
use crate::hash::{h_msg, prf_msg};
use crate::hypertree::{ht_sign, HypertreeSig};
use crate::keygen::SecretKey;
use crate::params::Params;
use rand::{CryptoRng, RngCore};

/// SLH-DSA signature.
#[derive(Clone, Debug)]
pub struct Signature {
    /// Randomizer (n bytes)
    pub r: Vec<u8>,
    /// FORS signature
    pub sig_fors: ForsSig,
    /// Hypertree signature
    pub sig_ht: HypertreeSig,
}

impl Signature {
    /// Serializes the signature to bytes.
    ///
    /// Format: R || SIG_FORS || SIG_HT
    pub fn to_bytes(&self, params: &Params) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(params.signature_size());

        // Randomizer R
        bytes.extend(&self.r);

        // FORS signature: k secret keys + k authentication paths
        for i in 0..params.k {
            bytes.extend(&self.sig_fors.sk[i]);
            for j in 0..params.a {
                bytes.extend(&self.sig_fors.auth[i][j]);
            }
        }

        // Hypertree signature: d XMSS signatures
        for layer_sig in &self.sig_ht.layers {
            // WOTS+ signature
            for sig_elem in &layer_sig.wots_sig {
                bytes.extend(sig_elem);
            }
            // Authentication path
            for auth_elem in &layer_sig.auth {
                bytes.extend(auth_elem);
            }
        }

        bytes
    }

    /// Deserializes a signature from bytes.
    pub fn from_bytes(bytes: &[u8], params: &Params) -> Option<Self> {
        if bytes.len() != params.signature_size() {
            return None;
        }

        let n = params.n;
        let mut offset = 0;

        // Read randomizer R
        let r = bytes[offset..offset + n].to_vec();
        offset += n;

        // Read FORS signature
        let mut fors_sk = Vec::with_capacity(params.k);
        let mut fors_auth = Vec::with_capacity(params.k);

        for _ in 0..params.k {
            // Secret key value
            let sk_val = bytes[offset..offset + n].to_vec();
            offset += n;
            fors_sk.push(sk_val);

            // Authentication path
            let mut auth_path = Vec::with_capacity(params.a);
            for _ in 0..params.a {
                let auth_node = bytes[offset..offset + n].to_vec();
                offset += n;
                auth_path.push(auth_node);
            }
            fors_auth.push(auth_path);
        }

        let sig_fors = ForsSig {
            sk: fors_sk,
            auth: fors_auth,
        };

        // Read hypertree signature
        let mut ht_layers = Vec::with_capacity(params.d);

        for _ in 0..params.d {
            // WOTS+ signature
            let mut wots_sig = Vec::with_capacity(params.wots_len);
            for _ in 0..params.wots_len {
                let sig_elem = bytes[offset..offset + n].to_vec();
                offset += n;
                wots_sig.push(sig_elem);
            }

            // Authentication path
            let mut auth = Vec::with_capacity(params.hp);
            for _ in 0..params.hp {
                let auth_elem = bytes[offset..offset + n].to_vec();
                offset += n;
                auth.push(auth_elem);
            }

            ht_layers.push(crate::xmss::XmssSig { wots_sig, auth });
        }

        let sig_ht = HypertreeSig { layers: ht_layers };

        Some(Signature { r, sig_fors, sig_ht })
    }
}

/// Signs a message using SLH-DSA (deterministic).
///
/// slh_sign(M, SK) - FIPS 205 Algorithm 18 (deterministic variant)
///
/// Uses SK.prf as the randomizer source for deterministic signing.
pub fn sign(sk: &SecretKey, message: &[u8]) -> Result<Signature> {
    sign_internal(sk, message, &sk.pk_seed)
}

/// Signs a message using SLH-DSA (randomized).
///
/// slh_sign(M, SK, opt_rand) - FIPS 205 Algorithm 18 (randomized variant)
///
/// Uses externally provided randomness for the randomizer.
pub fn sign_randomized<R: RngCore + CryptoRng>(
    rng: &mut R,
    sk: &SecretKey,
    message: &[u8],
) -> Result<Signature> {
    let mut opt_rand = vec![0u8; sk.params.n];
    rng.fill_bytes(&mut opt_rand);
    sign_internal(sk, message, &opt_rand)
}

/// Internal signing implementation.
fn sign_internal(sk: &SecretKey, message: &[u8], opt_rand: &[u8]) -> Result<Signature> {
    let params = &sk.params;

    // Step 1: Compute randomizer R = PRF_msg(SK.prf, OptRand, M)
    let r = prf_msg(&sk.sk_prf, opt_rand, message, params);

    // Step 2: Compute message digest and indices
    // H_msg(R, PK.seed, PK.root, M) -> (md, idx_tree, idx_leaf)
    let (md, idx_tree, idx_leaf) = h_msg(&r, &sk.pk_seed, &sk.pk_root, message, params);

    // Step 3: Generate FORS signature
    let mut fors_adrs = Address::new();
    fors_adrs.set_type(AddressType::ForsTree);
    fors_adrs.set_layer(0);
    fors_adrs.set_tree(idx_tree);
    fors_adrs.set_keypair(idx_leaf);

    let sig_fors = fors_sign(&md, &sk.sk_seed, &sk.pk_seed, &mut fors_adrs, params);

    // Step 4: Compute FORS public key for hypertree input
    let mut fors_pk_adrs = Address::new();
    fors_pk_adrs.set_type(AddressType::ForsTree);
    fors_pk_adrs.set_layer(0);
    fors_pk_adrs.set_tree(idx_tree);
    fors_pk_adrs.set_keypair(idx_leaf);

    let pk_fors = fors_pk_from_sig(&sig_fors, &md, &sk.pk_seed, &mut fors_pk_adrs, params);

    // Step 5: Generate hypertree signature on FORS public key
    let sig_ht = ht_sign(&pk_fors, &sk.sk_seed, &sk.pk_seed, idx_tree, idx_leaf, params);

    Ok(Signature { r, sig_fors, sig_ht })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::keygen;
    use crate::params::SLH_DSA_SHAKE_128F;
    use rand::rngs::OsRng;

    #[test]
    fn test_sign_produces_signature() {
        let params = SLH_DSA_SHAKE_128F;
        let (_, sk) = keygen(&mut OsRng, params);
        let message = b"Test message for signing";

        let sig = sign(&sk, message).unwrap();

        assert_eq!(sig.r.len(), params.n);
        assert_eq!(sig.sig_fors.sk.len(), params.k);
        assert_eq!(sig.sig_ht.layers.len(), params.d);
    }

    #[test]
    fn test_sign_deterministic() {
        let params = SLH_DSA_SHAKE_128F;
        let (_, sk) = keygen(&mut OsRng, params);
        let message = b"Test message for signing";

        let sig1 = sign(&sk, message).unwrap();
        let sig2 = sign(&sk, message).unwrap();

        // Same key + same message = same signature (deterministic)
        assert_eq!(sig1.r, sig2.r);
    }

    #[test]
    fn test_sign_randomized_different() {
        let params = SLH_DSA_SHAKE_128F;
        let (_, sk) = keygen(&mut OsRng, params);
        let message = b"Test message for signing";

        let sig1 = sign_randomized(&mut OsRng, &sk, message).unwrap();
        let sig2 = sign_randomized(&mut OsRng, &sk, message).unwrap();

        // Randomized signatures should be different
        assert_ne!(sig1.r, sig2.r);
    }

    #[test]
    fn test_signature_size() {
        let params = SLH_DSA_SHAKE_128F;
        let (_, sk) = keygen(&mut OsRng, params);
        let message = b"Test message";

        let sig = sign(&sk, message).unwrap();
        let bytes = sig.to_bytes(&params);

        assert_eq!(bytes.len(), params.signature_size());
    }

    #[test]
    fn test_signature_serialization_roundtrip() {
        let params = SLH_DSA_SHAKE_128F;
        let (_, sk) = keygen(&mut OsRng, params);
        let message = b"Test message";

        let sig = sign(&sk, message).unwrap();
        let bytes = sig.to_bytes(&params);
        let sig_recovered = Signature::from_bytes(&bytes, &params).unwrap();

        assert_eq!(sig.r, sig_recovered.r);
        assert_eq!(sig.sig_fors.sk.len(), sig_recovered.sig_fors.sk.len());
        assert_eq!(sig.sig_ht.layers.len(), sig_recovered.sig_ht.layers.len());
    }

    #[test]
    fn test_sign_different_messages() {
        let params = SLH_DSA_SHAKE_128F;
        let (_, sk) = keygen(&mut OsRng, params);

        let sig1 = sign(&sk, b"Message 1").unwrap();
        let sig2 = sign(&sk, b"Message 2").unwrap();

        // Different messages should produce different signatures
        assert_ne!(sig1.r, sig2.r);
    }

    #[test]
    fn test_invalid_signature_length() {
        let params = SLH_DSA_SHAKE_128F;
        let bytes = vec![0u8; params.signature_size() - 1];

        assert!(Signature::from_bytes(&bytes, &params).is_none());
    }
}
