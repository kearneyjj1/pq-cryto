//! # pqsigs-uov
//!
//! A Rust implementation of the UOV (Unbalanced Oil and Vinegar) post-quantum
//! digital signature scheme.
//!
//! ## Overview
//!
//! UOV is a multivariate quadratic (MQ) signature scheme that is believed to be
//! secure against quantum computers. It was submitted to the NIST Post-Quantum
//! Cryptography standardization process.
//!
//! This implementation is intended for **educational and experimental purposes only**.
//! It has not been audited for production use and may contain timing side-channels.
//!
//! ## Quick Start
//!
//! ```rust
//! use rand::rngs::OsRng;
//! use pqsigs_uov::{keygen::keygen, sign::sign, verify::verify, params::PARAMS_DEMO};
//!
//! // Generate a key pair
//! let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
//!
//! // Sign a message
//! let msg = b"Hello, post-quantum world!";
//! let sig = sign(&mut OsRng, &pk, &sk, msg).expect("signing should succeed");
//!
//! // Verify the signature
//! assert!(verify(&pk, msg, &sig).is_ok());
//! ```
//!
//! ## Parameter Sets
//!
//! The library provides several parameter sets:
//!
//! - [`params::PARAMS_DEMO`]: Small parameters for fast testing (NOT secure)
//! - [`params::PARAMS_NIST_L1`]: NIST Level 1 (~128-bit security)
//! - [`params::PARAMS_NIST_L3`]: NIST Level 3 (~192-bit security)
//! - [`params::PARAMS_NIST_L5`]: NIST Level 5 (~256-bit security)
//!
//! ## Modules
//!
//! - [`error`]: Error types for UOV operations
//! - [`field`]: GF(2^8) finite field arithmetic
//! - [`matrix`]: Matrix operations over GF(256)
//! - [`params`]: Parameter sets and validation
//! - [`keygen`]: Key generation
//! - [`sign`]: Signature generation
//! - [`verify`]: Signature verification
//!
//! ## Security Warning
//!
//! This implementation:
//! - Is NOT constant-time and may leak information through timing
//! - Has NOT been audited by security professionals
//! - Should NOT be used in production systems
//!
//! Use only for learning, experimentation, and research.

#![warn(missing_docs)]

pub mod error;
pub mod field;
pub mod keygen;
pub mod matrix;
pub mod params;
pub mod sign;
pub mod verify;

// Re-export commonly used types at crate root for convenience
pub use error::{Result, UovError};
pub use field::F;
pub use keygen::{keygen, PublicKey, QuadForm, SecretKey, Signature};
pub use params::{Params, PARAMS_DEMO, PARAMS_NIST_L1, PARAMS_NIST_L3, PARAMS_NIST_L5};
pub use sign::sign;
pub use verify::{verify, verify_bool};
