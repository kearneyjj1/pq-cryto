//! Error types for UOV signature operations.

use std::fmt;

/// Specific reasons why signature verification failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerificationFailure {
    /// The public map evaluation doesn't match the expected hash.
    HashMismatch,
    /// The signature point has wrong length.
    InvalidLength,
    /// The salt has wrong length.
    InvalidSaltLength,
    /// Generic verification failure (unspecified reason).
    Generic,
}

impl fmt::Display for VerificationFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VerificationFailure::HashMismatch => write!(f, "public map hash mismatch"),
            VerificationFailure::InvalidLength => write!(f, "invalid signature length"),
            VerificationFailure::InvalidSaltLength => write!(f, "invalid salt length"),
            VerificationFailure::Generic => write!(f, "verification failed"),
        }
    }
}

/// Errors that can occur during UOV signature operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UovError {
    /// Signature verification failed with a specific reason.
    VerificationFailed(VerificationFailure),

    /// Signature verification failed (legacy variant for compatibility).
    InvalidSignature,

    /// Signing failed after exhausting all retry attempts.
    SigningFailed {
        /// Number of attempts made before failure.
        attempts: u32,
    },

    /// The provided key pair is inconsistent or malformed.
    InvalidKeyPair {
        /// Description of the inconsistency.
        reason: &'static str,
    },

    /// The provided parameters are invalid.
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
}

impl fmt::Display for UovError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UovError::VerificationFailed(reason) => {
                write!(f, "signature verification failed: {}", reason)
            }
            UovError::InvalidSignature => {
                write!(f, "signature verification failed")
            }
            UovError::SigningFailed { attempts } => {
                write!(f, "signing failed after {} attempts", attempts)
            }
            UovError::InvalidKeyPair { reason } => {
                write!(f, "invalid key pair: {}", reason)
            }
            UovError::InvalidParams { reason } => {
                write!(f, "invalid parameters: {}", reason)
            }
            UovError::InvalidInput { field, reason } => {
                write!(f, "invalid input for '{}': {}", field, reason)
            }
        }
    }
}

impl std::error::Error for UovError {}

/// Result type alias for UOV operations.
pub type Result<T> = std::result::Result<T, UovError>;
