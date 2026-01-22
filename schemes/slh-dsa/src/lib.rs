//! SLH-DSA (Stateless Hash-based Digital Signature Algorithm) - FIPS 205
//!
//! This crate provides an implementation of SLH-DSA, the NIST post-quantum
//! digital signature standard formerly known as SPHINCS+.
//!
//! # Security Levels
//!
//! Six parameter sets are provided (SHAKE variants):
//!
//! | Parameter Set | Security | Signature Size | Variant |
//! |---------------|----------|----------------|---------|
//! | **SLH-DSA-SHAKE-128s** | NIST Level 1 | 7,856 B | Small |
//! | **SLH-DSA-SHAKE-128f** | NIST Level 1 | 17,088 B | Fast |
//! | **SLH-DSA-SHAKE-192s** | NIST Level 3 | 16,224 B | Small |
//! | **SLH-DSA-SHAKE-192f** | NIST Level 3 | 35,664 B | Fast |
//! | **SLH-DSA-SHAKE-256s** | NIST Level 5 | 29,792 B | Small |
//! | **SLH-DSA-SHAKE-256f** | NIST Level 5 | 49,856 B | Fast |
//!
//! The 's' (small) variants have smaller signatures but slower signing.
//! The 'f' (fast) variants have faster signing but larger signatures.
//!
//! # Example Usage
//!
//! ```rust
//! use pqsigs_slh_dsa::{keygen, sign, verify, params::SLH_DSA_SHAKE_128F};
//! use rand::rngs::OsRng;
//!
//! // Generate a key pair
//! let (public_key, secret_key) = keygen(&mut OsRng, SLH_DSA_SHAKE_128F);
//!
//! // Sign a message
//! let message = b"Hello, post-quantum world!";
//! let signature = sign(&secret_key, message).expect("signing should succeed");
//!
//! // Verify the signature
//! assert!(verify(&public_key, message, &signature).is_ok());
//! ```
//!
//! # Algorithm Overview
//!
//! SLH-DSA is a stateless hash-based signature scheme that uses:
//!
//! - **WOTS+**: Winternitz One-Time Signature as the base building block
//! - **XMSS**: Extended Merkle Signature Scheme (Merkle trees of WOTS+ keys)
//! - **FORS**: Forest of Random Subsets for signing the message digest
//! - **Hypertree**: Multi-layer tree structure connecting XMSS trees
//!
//! ## Key Generation
//! 1. Generate random seeds: SK.seed, SK.prf, PK.seed
//! 2. Compute the hypertree root using XMSS trees
//! 3. Public key: (PK.seed, PK.root)
//! 4. Secret key: (SK.seed, SK.prf, PK.seed, PK.root)
//!
//! ## Signing
//! 1. Compute randomizer R = PRF(SK.prf, message)
//! 2. Hash message to get FORS indices: H_msg(R, PK, message)
//! 3. Generate FORS signature on message digest
//! 4. Generate hypertree signature on FORS public key
//! 5. Return signature: (R, FORS_sig, HT_sig)
//!
//! ## Verification
//! 1. Recompute message digest from R
//! 2. Reconstruct FORS public key from signature
//! 3. Verify hypertree signature matches PK.root
//!
//! # Security Warning
//!
//! This implementation:
//! - Is NOT constant-time and may leak information through timing
//! - Has NOT been audited by security professionals
//! - Should NOT be used in production systems
//!
//! Use only for learning, experimentation, and research.
//!
//! # References
//!
//! - FIPS 205: Stateless Hash-Based Digital Signature Standard
//! - <https://csrc.nist.gov/pubs/fips/205/final>

#![warn(missing_docs)]
#![warn(rust_2018_idioms)]

pub mod address;
pub mod error;
pub mod fors;
pub mod hash;
pub mod hypertree;
pub mod keygen;
pub mod params;
pub mod sign;
pub mod verify;
pub mod wots;
pub mod xmss;

// Re-export main types and functions for convenience
pub use error::{Result, SlhDsaError};
pub use keygen::{keygen, keygen_internal, PublicKey, SecretKey};
pub use params::{
    Params, SLH_DSA_SHAKE_128F, SLH_DSA_SHAKE_128S, SLH_DSA_SHAKE_192F, SLH_DSA_SHAKE_192S,
    SLH_DSA_SHAKE_256F, SLH_DSA_SHAKE_256S,
};
pub use sign::{sign, sign_randomized, Signature};
pub use verify::{verify, verify_bool};

#[cfg(test)]
mod integration_tests {
    use super::*;
    use rand::rngs::OsRng;

    #[test]
    fn test_full_roundtrip_128f() {
        let (pk, sk) = keygen(&mut OsRng, SLH_DSA_SHAKE_128F);
        let message = b"Integration test message for SLH-DSA-SHAKE-128f";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_key_sizes() {
        // Test all parameter sets have correct key sizes
        for params in [
            SLH_DSA_SHAKE_128S,
            SLH_DSA_SHAKE_128F,
            SLH_DSA_SHAKE_192S,
            SLH_DSA_SHAKE_192F,
            SLH_DSA_SHAKE_256S,
            SLH_DSA_SHAKE_256F,
        ] {
            let (pk, sk) = keygen(&mut OsRng, params);

            assert_eq!(pk.to_bytes().len(), params.public_key_size());
            assert_eq!(sk.to_bytes().len(), params.secret_key_size());
        }
    }

    #[test]
    fn test_signature_sizes() {
        // Test fast variants (128f, 192f, 256f) - they're faster
        for params in [
            SLH_DSA_SHAKE_128F,
            SLH_DSA_SHAKE_192F,
            SLH_DSA_SHAKE_256F,
        ] {
            let (_, sk) = keygen(&mut OsRng, params);
            let sig = sign(&sk, b"test").unwrap();

            assert_eq!(sig.to_bytes(&params).len(), params.signature_size());
        }
    }

    #[test]
    fn test_cross_key_rejection() {
        let (pk1, _) = keygen(&mut OsRng, SLH_DSA_SHAKE_128F);
        let (_, sk2) = keygen(&mut OsRng, SLH_DSA_SHAKE_128F);

        let sig = sign(&sk2, b"test").unwrap();

        // Signature from sk2 should not verify with pk1
        assert!(verify(&pk1, b"test", &sig).is_err());
    }

    #[test]
    fn test_message_modification_detected() {
        let (pk, sk) = keygen(&mut OsRng, SLH_DSA_SHAKE_128F);

        let sig = sign(&sk, b"original message").unwrap();

        // Any modification to the message should fail verification
        assert!(verify(&pk, b"original messag", &sig).is_err());
        assert!(verify(&pk, b"Original message", &sig).is_err());
        assert!(verify(&pk, b"original message!", &sig).is_err());
    }

    #[test]
    fn test_deterministic_keygen() {
        let params = SLH_DSA_SHAKE_128F;
        let sk_seed = [42u8; 16];
        let sk_prf = [43u8; 16];
        let pk_seed = [44u8; 16];

        let (pk1, sk1) = keygen_internal(&sk_seed, &sk_prf, &pk_seed, params);
        let (pk2, sk2) = keygen_internal(&sk_seed, &sk_prf, &pk_seed, params);

        assert_eq!(pk1.to_bytes(), pk2.to_bytes());
        assert_eq!(sk1.to_bytes(), sk2.to_bytes());
    }

    #[test]
    fn test_deterministic_signing() {
        let (pk, sk) = keygen(&mut OsRng, SLH_DSA_SHAKE_128F);
        let message = b"deterministic test";

        let sig1 = sign(&sk, message).unwrap();
        let sig2 = sign(&sk, message).unwrap();

        // Same key + same message = same signature (deterministic signing)
        assert_eq!(sig1.r, sig2.r);
        assert!(verify(&pk, message, &sig1).is_ok());
        assert!(verify(&pk, message, &sig2).is_ok());
    }

    #[test]
    fn test_randomized_signing() {
        let (pk, sk) = keygen(&mut OsRng, SLH_DSA_SHAKE_128F);
        let message = b"randomized test";

        let sig1 = sign_randomized(&mut OsRng, &sk, message).unwrap();
        let sig2 = sign_randomized(&mut OsRng, &sk, message).unwrap();

        // Randomized signing produces different signatures
        assert_ne!(sig1.r, sig2.r);

        // Both should still verify
        assert!(verify(&pk, message, &sig1).is_ok());
        assert!(verify(&pk, message, &sig2).is_ok());
    }

    #[test]
    fn test_signature_serialization() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"serialization test";

        let sig = sign(&sk, message).unwrap();
        let bytes = sig.to_bytes(&params);
        let sig_recovered = Signature::from_bytes(&bytes, &params).unwrap();

        // Recovered signature should still verify
        assert!(verify(&pk, message, &sig_recovered).is_ok());
    }

    #[test]
    fn test_key_serialization() {
        let params = SLH_DSA_SHAKE_128F;
        let (pk, sk) = keygen(&mut OsRng, params);
        let message = b"key serialization test";

        // Sign with original keys
        let sig = sign(&sk, message).unwrap();

        // Serialize and deserialize keys
        let pk_bytes = pk.to_bytes();
        let sk_bytes = sk.to_bytes();

        let pk_recovered = PublicKey::from_bytes(&pk_bytes, params).unwrap();
        let sk_recovered = SecretKey::from_bytes(&sk_bytes, params).unwrap();

        // Sign with recovered secret key
        let sig2 = sign(&sk_recovered, message).unwrap();

        // Verify with recovered public key
        assert!(verify(&pk_recovered, message, &sig).is_ok());
        assert!(verify(&pk_recovered, message, &sig2).is_ok());
    }

    // Slow tests - uncomment for thorough testing
    /*
    #[test]
    fn test_full_roundtrip_128s() {
        let (pk, sk) = keygen(&mut OsRng, SLH_DSA_SHAKE_128S);
        let message = b"SLH-DSA-SHAKE-128s test";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
        assert_eq!(sig.to_bytes(&SLH_DSA_SHAKE_128S).len(), SLH_DSA_SHAKE_128S.signature_size());
    }

    #[test]
    fn test_full_roundtrip_192s() {
        let (pk, sk) = keygen(&mut OsRng, SLH_DSA_SHAKE_192S);
        let message = b"SLH-DSA-SHAKE-192s test";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_full_roundtrip_256s() {
        let (pk, sk) = keygen(&mut OsRng, SLH_DSA_SHAKE_256S);
        let message = b"SLH-DSA-SHAKE-256s test";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
    }
    */
}
