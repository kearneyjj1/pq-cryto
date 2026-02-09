//! # pqsigs-fn-dsa
//!
//! A Rust implementation of FN-DSA (FALCON), a lattice-based post-quantum
//! digital signature scheme.
//!
//! ## Overview
//!
//! FALCON (Fast Fourier Lattice-based Compact Signatures over NTRU) is a
//! hash-and-sign signature scheme based on NTRU lattices. It was selected
//! by NIST for standardization as part of the post-quantum cryptography
//! standardization process.
//!
//! This implementation is intended for **educational and experimental purposes only**.
//! It has not been audited for production use and may contain timing side-channels.
//!
//! ## Parameter Sets
//!
//! - [`params::FALCON_512`]: NIST Level 1 (~128-bit security)
//! - [`params::FALCON_1024`]: NIST Level 5 (~256-bit security)
//!
//! ## Platform Requirements
//!
//! This implementation requires IEEE 754 double-precision floating-point
//! with round-to-nearest-even mode for the signing path. On x86/x86_64,
//! SSE2 is required (default on x86_64, enforced via `compile_error!`).
//! Signing results may differ across platforms with different FP behavior
//! (x87 extended precision, FMA instructions).
//!
//! ## Security Warning
//!
//! This implementation:
//! - Has partial constant-time protections (BerExp, RCDT, field arithmetic)
//! - Uses floating-point FFT which requires platform constraints (see above)
//! - Has NOT been audited by security professionals
//! - Should NOT be used in production systems
//!
//! Use only for learning, experimentation, and research.

#![warn(missing_docs)]
#![warn(rust_2018_idioms)]

pub mod error;
pub(crate) mod fft;
pub mod fft_tree;
pub(crate) mod field;
pub(crate) mod gaussian;
pub mod hash;
pub mod keygen;
pub(crate) mod ntru;
pub mod params;
pub(crate) mod poly;
pub(crate) mod sampler;
pub mod sign;
pub mod verify;

pub mod packing;

// Re-export main types for convenience
pub use error::{FnDsaError, Result};
pub use keygen::{KeyPair, PublicKey, SecretKey, keygen, keygen_512, keygen_1024, keygen_with_seed};
pub use params::{Params, FALCON_512, FALCON_1024};
pub use sign::{Signature, sign};
pub use verify::verify;

// Re-export packing functions
pub use packing::{
    encode_public_key, decode_public_key,
    encode_secret_key, decode_secret_key,
    encode_signature, decode_signature,
    encode_signature_raw, encode_signature_compressed,
    decode_nist_public_key, parse_nist_signed_message,
    to_hex, from_hex,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_params_available() {
        assert_eq!(FALCON_512.n, 512);
        assert_eq!(FALCON_1024.n, 1024);
    }
}
