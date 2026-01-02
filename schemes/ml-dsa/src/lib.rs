//! ML-DSA (Module-Lattice Digital Signature Algorithm) - FIPS 204
//!
//! This crate provides a reference implementation of ML-DSA, the NIST
//! post-quantum digital signature standard formerly known as CRYSTALS-Dilithium.
//!
//! # Security Levels
//!
//! Three parameter sets are provided:
//! - **ML-DSA-44**: NIST Level 2 (~128-bit security)
//! - **ML-DSA-65**: NIST Level 3 (~192-bit security)
//! - **ML-DSA-87**: NIST Level 5 (~256-bit security)
//!
//! # Example Usage
//!
//! ```rust
//! use pqsigs_ml_dsa::{keygen, sign, verify, params::ML_DSA_44};
//! use rand::rngs::OsRng;
//!
//! // Generate a key pair
//! let (public_key, secret_key) = keygen(&mut OsRng, ML_DSA_44);
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
//! ML-DSA is a lattice-based signature scheme based on the "Fiat-Shamir with
//! Aborts" paradigm. The security is based on the hardness of the Module
//! Learning With Errors (M-LWE) and Module Short Integer Solution (M-SIS)
//! problems.
//!
//! ## Key Generation
//! 1. Generate random seeds ρ, ρ', K
//! 2. Expand matrix A from ρ using SHAKE128
//! 3. Sample secret vectors s1, s2 with small coefficients from ρ'
//! 4. Compute t = As1 + s2
//! 5. Split t into t1 (high bits) and t0 (low bits) using Power2Round
//! 6. Public key: (ρ, t1), Secret key: (ρ, K, tr, s1, s2, t0)
//!
//! ## Signing
//! 1. Compute message hash μ = H(tr || M)
//! 2. Sample masking vector y
//! 3. Compute w = Ay and extract high bits w1
//! 4. Compute challenge c = H(μ || w1)
//! 5. Compute z = y + cs1
//! 6. If ||z||∞ or ||LowBits(w - cs2)||∞ too large, restart
//! 7. Compute hints h for recovering w1 from w - cs2 + ct0
//! 8. Return signature (c_tilde, z, h)
//!
//! ## Verification
//! 1. Recompute challenge c from c_tilde
//! 2. Compute w' = Az - ct1·2^d
//! 3. Use hints to recover w1' from w'
//! 4. Verify H(μ || w1') == c_tilde
//!
//! # References
//!
//! - FIPS 204: Module-Lattice-Based Digital Signature Standard
//! - <https://csrc.nist.gov/pubs/fips/204/final>

#![warn(missing_docs)]
#![warn(rust_2018_idioms)]

pub mod error;
pub mod field;
pub mod keygen;
pub mod ntt;
pub mod packing;
pub mod params;
pub mod poly;
pub mod polyvec;
pub mod reduce;
pub mod rounding;
pub mod sampling;
pub mod sign;
pub mod verify;

// Re-export main types and functions for convenience
pub use error::{MlDsaError, Result};
pub use keygen::{keygen, keygen_internal, PublicKey, SecretKey};
pub use params::{Params, ML_DSA_44, ML_DSA_65, ML_DSA_87};
pub use sign::{sign, Signature};
pub use verify::{verify, verify_bool};

#[cfg(test)]
mod integration_tests {
    use super::*;
    use rand::rngs::OsRng;

    #[test]
    fn test_full_roundtrip_ml_dsa_44() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
        let message = b"Integration test message for ML-DSA-44";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_full_roundtrip_ml_dsa_65() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_65);
        let message = b"Integration test message for ML-DSA-65";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_full_roundtrip_ml_dsa_87() {
        let (pk, sk) = keygen(&mut OsRng, ML_DSA_87);
        let message = b"Integration test message for ML-DSA-87";

        let sig = sign(&sk, message).expect("signing should succeed");
        assert!(verify(&pk, message, &sig).is_ok());
    }

    #[test]
    fn test_key_sizes() {
        let (pk44, sk44) = keygen(&mut OsRng, ML_DSA_44);
        let (pk65, sk65) = keygen(&mut OsRng, ML_DSA_65);
        let (pk87, sk87) = keygen(&mut OsRng, ML_DSA_87);

        // Verify key sizes match specification
        assert_eq!(pk44.to_bytes().len(), ML_DSA_44.public_key_size());
        assert_eq!(sk44.to_bytes().len(), ML_DSA_44.secret_key_size());

        assert_eq!(pk65.to_bytes().len(), ML_DSA_65.public_key_size());
        assert_eq!(sk65.to_bytes().len(), ML_DSA_65.secret_key_size());

        assert_eq!(pk87.to_bytes().len(), ML_DSA_87.public_key_size());
        assert_eq!(sk87.to_bytes().len(), ML_DSA_87.secret_key_size());
    }

    #[test]
    fn test_signature_sizes() {
        let (_, sk44) = keygen(&mut OsRng, ML_DSA_44);
        let (_, sk65) = keygen(&mut OsRng, ML_DSA_65);
        let (_, sk87) = keygen(&mut OsRng, ML_DSA_87);

        let sig44 = sign(&sk44, b"test").unwrap();
        let sig65 = sign(&sk65, b"test").unwrap();
        let sig87 = sign(&sk87, b"test").unwrap();

        assert_eq!(
            sig44
                .to_bytes(ML_DSA_44.gamma1, ML_DSA_44.omega, ML_DSA_44.k)
                .len(),
            ML_DSA_44.signature_size()
        );
        assert_eq!(
            sig65
                .to_bytes(ML_DSA_65.gamma1, ML_DSA_65.omega, ML_DSA_65.k)
                .len(),
            ML_DSA_65.signature_size()
        );
        assert_eq!(
            sig87
                .to_bytes(ML_DSA_87.gamma1, ML_DSA_87.omega, ML_DSA_87.k)
                .len(),
            ML_DSA_87.signature_size()
        );
    }

    #[test]
    fn test_cross_parameter_rejection() {
        // Keys from different parameter sets should not work together
        let (pk44, _) = keygen(&mut OsRng, ML_DSA_44);
        let (_, sk65) = keygen(&mut OsRng, ML_DSA_65);

        let sig = sign(&sk65, b"test").unwrap();

        // This will fail because dimensions don't match
        // The verify function should handle this gracefully
        assert!(verify(&pk44, b"test", &sig).is_err());
    }

    #[test]
    fn test_deterministic_signing() {
        let seed = [42u8; 32];
        let (pk, sk) = keygen_internal(&seed, ML_DSA_44);
        let message = b"deterministic test";

        let sig1 = sign(&sk, message).unwrap();
        let sig2 = sign(&sk, message).unwrap();

        // Same key + same message = same signature
        assert_eq!(sig1.c_tilde, sig2.c_tilde);
        assert!(verify(&pk, message, &sig1).is_ok());
        assert!(verify(&pk, message, &sig2).is_ok());
    }
}
