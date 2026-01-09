//! SQI-SIGN: Isogeny-Based Digital Signature Scheme
//!
//! This crate provides an educational implementation of SQI-SIGN, a post-quantum
//! digital signature scheme based on the hardness of finding isogenies between
//! supersingular elliptic curves.
//!
//! # Security Basis
//!
//! SQI-SIGN's security relies on the computational difficulty of:
//! - The Isogeny Problem: Finding an isogeny between two supersingular curves
//! - The Endomorphism Ring Problem: Computing the endomorphism ring of a curve
//!
//! # Algorithm Overview
//!
//! SQI-SIGN uses a commitment-challenge-response paradigm:
//! 1. **Key Generation**: Generate a secret isogeny φ: E₀ → E_A
//! 2. **Signing**: Commit to a random isogeny, receive challenge, respond using
//!    quaternion algebra computations to find a "short" isogeny
//! 3. **Verification**: Check that the response isogeny has the correct properties
//!
//! # Warning
//!
//! This implementation is for educational purposes only and has NOT been audited
//! for production use. See SECURITY.md for details.
//!
//! # Example
//!
//! ```ignore
//! use pqsigs_sqi_sign::{keygen, sign, verify, params::SQISIGN_NIST_I};
//! use rand::rngs::OsRng;
//!
//! let (pk, sk) = keygen(&mut OsRng, SQISIGN_NIST_I);
//! let sig = sign(&mut OsRng, &sk, b"message").unwrap();
//! assert!(verify(&pk, b"message", &sig).is_ok());
//! ```

#![warn(missing_docs)]

pub mod error;
pub mod params;
pub mod field;
pub mod fp2;
pub mod curve;
pub mod isogeny;
pub mod quaternion;
pub mod ideal;
pub mod keygen;
pub mod sign;
pub mod verify;

pub use error::{SqiSignError, Result};
pub use keygen::{keygen, PublicKey, SecretKey};
pub use sign::{sign, Signature};
pub use verify::verify;
