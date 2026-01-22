//! Key generation for SLH-DSA.
//!
//! Implements SLH-DSA key generation (FIPS 205 Algorithm 17).
//!
//! # Key Structure
//!
//! - **Public Key**: PK.seed || PK.root (2n bytes)
//!   - PK.seed: Random seed for generating public parameters
//!   - PK.root: Root of the hypertree (the actual public key value)
//!
//! - **Secret Key**: SK.seed || SK.prf || PK.seed || PK.root (4n bytes)
//!   - SK.seed: Secret seed for deriving all secret values
//!   - SK.prf: Secret key for signature randomization
//!   - PK.seed, PK.root: Copies of public key components for efficiency

use crate::address::Address;
use crate::params::Params;
use crate::xmss::xmss_node;
use rand::{CryptoRng, RngCore};
use zeroize::Zeroize;

/// SLH-DSA public key.
#[derive(Clone, Debug)]
pub struct PublicKey {
    /// Public seed (n bytes)
    pub pk_seed: Vec<u8>,
    /// Root of the hypertree (n bytes)
    pub pk_root: Vec<u8>,
    /// Parameter set
    pub params: Params,
}

impl PublicKey {
    /// Serializes the public key to bytes.
    ///
    /// Format: PK.seed || PK.root
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.params.public_key_size());
        bytes.extend(&self.pk_seed);
        bytes.extend(&self.pk_root);
        bytes
    }

    /// Deserializes a public key from bytes.
    pub fn from_bytes(bytes: &[u8], params: Params) -> Option<Self> {
        if bytes.len() != params.public_key_size() {
            return None;
        }

        let pk_seed = bytes[0..params.n].to_vec();
        let pk_root = bytes[params.n..2 * params.n].to_vec();

        Some(PublicKey {
            pk_seed,
            pk_root,
            params,
        })
    }
}

/// SLH-DSA secret key.
///
/// # Security
///
/// This struct implements `Drop` to zeroize secret key material when dropped.
#[derive(Clone, Debug)]
pub struct SecretKey {
    /// Secret seed for deriving all secret values (n bytes)
    pub sk_seed: Vec<u8>,
    /// Secret key for signature randomization (n bytes)
    pub sk_prf: Vec<u8>,
    /// Public seed (copy for efficiency, n bytes)
    pub pk_seed: Vec<u8>,
    /// Public root (copy for efficiency, n bytes)
    pub pk_root: Vec<u8>,
    /// Parameter set
    pub params: Params,
}

impl Drop for SecretKey {
    fn drop(&mut self) {
        // Zeroize all secret material
        self.sk_seed.zeroize();
        self.sk_prf.zeroize();
        // Also zeroize public copies for defense in depth
        self.pk_seed.zeroize();
        self.pk_root.zeroize();
    }
}

impl SecretKey {
    /// Serializes the secret key to bytes.
    ///
    /// Format: SK.seed || SK.prf || PK.seed || PK.root
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.params.secret_key_size());
        bytes.extend(&self.sk_seed);
        bytes.extend(&self.sk_prf);
        bytes.extend(&self.pk_seed);
        bytes.extend(&self.pk_root);
        bytes
    }

    /// Deserializes a secret key from bytes.
    pub fn from_bytes(bytes: &[u8], params: Params) -> Option<Self> {
        if bytes.len() != params.secret_key_size() {
            return None;
        }

        let n = params.n;
        let sk_seed = bytes[0..n].to_vec();
        let sk_prf = bytes[n..2 * n].to_vec();
        let pk_seed = bytes[2 * n..3 * n].to_vec();
        let pk_root = bytes[3 * n..4 * n].to_vec();

        Some(SecretKey {
            sk_seed,
            sk_prf,
            pk_seed,
            pk_root,
            params,
        })
    }

    /// Returns the public key corresponding to this secret key.
    pub fn public_key(&self) -> PublicKey {
        PublicKey {
            pk_seed: self.pk_seed.clone(),
            pk_root: self.pk_root.clone(),
            params: self.params,
        }
    }
}

/// Generates an SLH-DSA key pair.
///
/// slh_keygen() - FIPS 205 Algorithm 17
///
/// # Arguments
/// * `rng` - Cryptographically secure random number generator
/// * `params` - Parameter set to use
///
/// # Returns
/// A tuple of (public_key, secret_key).
pub fn keygen<R: RngCore + CryptoRng>(rng: &mut R, params: Params) -> (PublicKey, SecretKey) {
    let mut sk_seed = vec![0u8; params.n];
    let mut sk_prf = vec![0u8; params.n];
    let mut pk_seed = vec![0u8; params.n];

    rng.fill_bytes(&mut sk_seed);
    rng.fill_bytes(&mut sk_prf);
    rng.fill_bytes(&mut pk_seed);

    keygen_internal(&sk_seed, &sk_prf, &pk_seed, params)
}

/// Internal key generation from seeds (for deterministic testing).
///
/// # Arguments
/// * `sk_seed` - Secret seed (n bytes)
/// * `sk_prf` - PRF key (n bytes)
/// * `pk_seed` - Public seed (n bytes)
/// * `params` - Parameter set
pub fn keygen_internal(
    sk_seed: &[u8],
    sk_prf: &[u8],
    pk_seed: &[u8],
    params: Params,
) -> (PublicKey, SecretKey) {
    // Compute the root of the hypertree
    // The root is at layer d-1, tree 0, node 0 at height h'
    let mut adrs = Address::new();
    adrs.set_layer((params.d - 1) as u32);
    adrs.set_tree(0);

    let pk_root = xmss_node(sk_seed, pk_seed, 0, params.hp as u32, &mut adrs, &params);

    let pk = PublicKey {
        pk_seed: pk_seed.to_vec(),
        pk_root: pk_root.clone(),
        params,
    };

    let sk = SecretKey {
        sk_seed: sk_seed.to_vec(),
        sk_prf: sk_prf.to_vec(),
        pk_seed: pk_seed.to_vec(),
        pk_root,
        params,
    };

    (pk, sk)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{
        SLH_DSA_SHAKE_128F, SLH_DSA_SHAKE_128S, SLH_DSA_SHAKE_192F, SLH_DSA_SHAKE_192S,
        SLH_DSA_SHAKE_256F, SLH_DSA_SHAKE_256S,
    };
    use rand::rngs::OsRng;

    #[test]
    fn test_keygen_128f() {
        let (pk, sk) = keygen(&mut OsRng, SLH_DSA_SHAKE_128F);

        assert_eq!(pk.pk_seed.len(), SLH_DSA_SHAKE_128F.n);
        assert_eq!(pk.pk_root.len(), SLH_DSA_SHAKE_128F.n);
        assert_eq!(sk.sk_seed.len(), SLH_DSA_SHAKE_128F.n);
        assert_eq!(sk.sk_prf.len(), SLH_DSA_SHAKE_128F.n);
    }

    #[test]
    fn test_keygen_deterministic() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = vec![1u8; params.n];
        let sk_prf = vec![2u8; params.n];
        let pk_seed = vec![3u8; params.n];

        let (pk1, sk1) = keygen_internal(&sk_seed, &sk_prf, &pk_seed, params);
        let (pk2, sk2) = keygen_internal(&sk_seed, &sk_prf, &pk_seed, params);

        assert_eq!(pk1.pk_seed, pk2.pk_seed);
        assert_eq!(pk1.pk_root, pk2.pk_root);
        assert_eq!(sk1.sk_seed, sk2.sk_seed);
        assert_eq!(sk1.sk_prf, sk2.sk_prf);
    }

    #[test]
    fn test_keygen_different_seeds() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_prf = vec![2u8; params.n];
        let pk_seed = vec![3u8; params.n];

        let sk_seed1 = vec![1u8; params.n];
        let sk_seed2 = vec![4u8; params.n];

        let (pk1, _) = keygen_internal(&sk_seed1, &sk_prf, &pk_seed, params);
        let (pk2, _) = keygen_internal(&sk_seed2, &sk_prf, &pk_seed, params);

        // Different sk_seed should give different pk_root
        assert_ne!(pk1.pk_root, pk2.pk_root);
    }

    #[test]
    fn test_public_key_size() {
        for params in [
            SLH_DSA_SHAKE_128S,
            SLH_DSA_SHAKE_128F,
            SLH_DSA_SHAKE_192S,
            SLH_DSA_SHAKE_192F,
            SLH_DSA_SHAKE_256S,
            SLH_DSA_SHAKE_256F,
        ] {
            let (pk, _) = keygen(&mut OsRng, params);
            let bytes = pk.to_bytes();
            assert_eq!(bytes.len(), params.public_key_size());
        }
    }

    #[test]
    fn test_secret_key_size() {
        for params in [
            SLH_DSA_SHAKE_128S,
            SLH_DSA_SHAKE_128F,
            SLH_DSA_SHAKE_192S,
            SLH_DSA_SHAKE_192F,
            SLH_DSA_SHAKE_256S,
            SLH_DSA_SHAKE_256F,
        ] {
            let (_, sk) = keygen(&mut OsRng, params);
            let bytes = sk.to_bytes();
            assert_eq!(bytes.len(), params.secret_key_size());
        }
    }

    #[test]
    fn test_public_key_serialization_roundtrip() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, _) = keygen(&mut OsRng, params);

        let bytes = pk.to_bytes();
        let pk_recovered = PublicKey::from_bytes(&bytes, params).unwrap();

        assert_eq!(pk.pk_seed, pk_recovered.pk_seed);
        assert_eq!(pk.pk_root, pk_recovered.pk_root);
    }

    #[test]
    fn test_secret_key_serialization_roundtrip() {
        let params = SLH_DSA_SHAKE_128F;
        let (_, sk) = keygen(&mut OsRng, params);

        let bytes = sk.to_bytes();
        let sk_recovered = SecretKey::from_bytes(&bytes, params).unwrap();

        assert_eq!(sk.sk_seed, sk_recovered.sk_seed);
        assert_eq!(sk.sk_prf, sk_recovered.sk_prf);
        assert_eq!(sk.pk_seed, sk_recovered.pk_seed);
        assert_eq!(sk.pk_root, sk_recovered.pk_root);
    }

    #[test]
    fn test_secret_key_public_key() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);

        let pk_from_sk = sk.public_key();

        assert_eq!(pk.pk_seed, pk_from_sk.pk_seed);
        assert_eq!(pk.pk_root, pk_from_sk.pk_root);
    }

    #[test]
    fn test_invalid_public_key_length() {
        let params = SLH_DSA_SHAKE_128F;
        let bytes = vec![0u8; params.public_key_size() - 1];
        assert!(PublicKey::from_bytes(&bytes, params).is_none());
    }

    #[test]
    fn test_invalid_secret_key_length() {
        let params = SLH_DSA_SHAKE_128F;
        let bytes = vec![0u8; params.secret_key_size() + 1];
        assert!(SecretKey::from_bytes(&bytes, params).is_none());
    }
}
