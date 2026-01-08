//! Error types for ML-DSA signature operations.

use std::fmt;

/// Specific reasons why signature verification failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerificationFailure {
    /// The signature's z vector norm exceeds the allowed bound.
    NormBoundExceeded,
    /// The computed challenge hash doesn't match the signature's challenge.
    ChallengeMismatch,
    /// The hint count exceeds the maximum allowed (omega).
    TooManyHints,
    /// The hint vector format is invalid.
    MalformedHints,
    /// Generic verification failure (unspecified reason).
    Generic,
}

impl fmt::Display for VerificationFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VerificationFailure::NormBoundExceeded => write!(f, "z norm bound exceeded"),
            VerificationFailure::ChallengeMismatch => write!(f, "challenge hash mismatch"),
            VerificationFailure::TooManyHints => write!(f, "too many hints"),
            VerificationFailure::MalformedHints => write!(f, "malformed hint vector"),
            VerificationFailure::Generic => write!(f, "verification failed"),
        }
    }
}

/// Errors that can occur during ML-DSA signature operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MlDsaError {
    /// Signature verification failed with a specific reason.
    VerificationFailed(VerificationFailure),

    /// Signature verification failed (legacy variant for compatibility).
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
            MlDsaError::VerificationFailed(reason) => {
                write!(f, "signature verification failed: {}", reason)
            }
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
