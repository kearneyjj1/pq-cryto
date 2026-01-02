//! Error types for ML-DSA signature operations.

use std::fmt;

/// Errors that can occur during ML-DSA signature operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MlDsaError {
    /// Signature verification failed.
    InvalidSignature,

    /// Signing failed after exhausting retry attempts.
    /// This can happen if the norm checks consistently fail.
    SigningFailed {
        /// Number of attempts made before failure.
        attempts: u32,
    },

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

    /// A hint vector was invalid (too many hints).
    InvalidHint,
}

impl fmt::Display for MlDsaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MlDsaError::InvalidSignature => {
                write!(f, "signature verification failed")
            }
            MlDsaError::SigningFailed { attempts } => {
                write!(f, "signing failed after {} attempts", attempts)
            }
            MlDsaError::InvalidKey { reason } => {
                write!(f, "invalid key: {}", reason)
            }
            MlDsaError::InvalidParams { reason } => {
                write!(f, "invalid parameters: {}", reason)
            }
            MlDsaError::InvalidInput { field, reason } => {
                write!(f, "invalid input for '{}': {}", field, reason)
            }
            MlDsaError::DecodingError { context } => {
                write!(f, "decoding error: {}", context)
            }
            MlDsaError::InvalidHint => {
                write!(f, "invalid hint vector")
            }
        }
    }
}

impl std::error::Error for MlDsaError {}

/// Result type alias for ML-DSA operations.
pub type Result<T> = std::result::Result<T, MlDsaError>;
