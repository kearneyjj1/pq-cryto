//! Error types for SLH-DSA signature operations.

use std::fmt;

/// Errors that can occur during SLH-DSA signature operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SlhDsaError {
    /// Signature verification failed.
    InvalidSignature,

    /// The provided key is malformed or invalid.
    InvalidKey {
        /// Description of the key issue.
        reason: &'static str,
    },

    /// Invalid parameter set specified.
    InvalidParams {
        /// Description of why the parameters are invalid.
        reason: &'static str,
    },

    /// Invalid input was provided to a function.
    InvalidInput {
        /// The name of the invalid field/parameter.
        field: &'static str,
        /// Description of why the input is invalid.
        reason: &'static str,
    },

    /// Decoding/unpacking failed.
    DecodingError {
        /// What was being decoded.
        context: &'static str,
    },
}

impl fmt::Display for SlhDsaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SlhDsaError::InvalidSignature => {
                write!(f, "signature verification failed")
            }
            SlhDsaError::InvalidKey { reason } => {
                write!(f, "invalid key: {}", reason)
            }
            SlhDsaError::InvalidParams { reason } => {
                write!(f, "invalid parameters: {}", reason)
            }
            SlhDsaError::InvalidInput { field, reason } => {
                write!(f, "invalid input for '{}': {}", field, reason)
            }
            SlhDsaError::DecodingError { context } => {
                write!(f, "decoding error: {}", context)
            }
        }
    }
}

impl std::error::Error for SlhDsaError {}

/// Result type alias for SLH-DSA operations.
pub type Result<T> = std::result::Result<T, SlhDsaError>;
