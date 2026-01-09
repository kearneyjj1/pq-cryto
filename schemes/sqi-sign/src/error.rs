//! Error types for SQI-SIGN.

use std::fmt;

/// Result type alias for SQI-SIGN operations.
pub type Result<T> = std::result::Result<T, SqiSignError>;

/// Errors that can occur during SQI-SIGN operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SqiSignError {
    /// Key generation failed.
    KeyGenFailed {
        /// Reason for failure.
        reason: &'static str,
    },
    /// Signing failed after maximum attempts.
    SigningFailed {
        /// Number of attempts made.
        attempts: u32,
    },
    /// Signature verification failed.
    InvalidSignature,
    /// Invalid input parameters.
    InvalidParams {
        /// Description of the invalid parameter.
        reason: &'static str,
    },
    /// Isogeny computation failed.
    IsogenyFailed {
        /// Reason for failure.
        reason: &'static str,
    },
    /// Quaternion algebra computation failed.
    QuaternionError {
        /// Reason for failure.
        reason: &'static str,
    },
    /// Ideal computation failed.
    IdealError {
        /// Reason for failure.
        reason: &'static str,
    },
}

impl fmt::Display for SqiSignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::KeyGenFailed { reason } => write!(f, "Key generation failed: {}", reason),
            Self::SigningFailed { attempts } => {
                write!(f, "Signing failed after {} attempts", attempts)
            }
            Self::InvalidSignature => write!(f, "Invalid signature"),
            Self::InvalidParams { reason } => write!(f, "Invalid parameters: {}", reason),
            Self::IsogenyFailed { reason } => write!(f, "Isogeny computation failed: {}", reason),
            Self::QuaternionError { reason } => {
                write!(f, "Quaternion algebra error: {}", reason)
            }
            Self::IdealError { reason } => write!(f, "Ideal computation error: {}", reason),
        }
    }
}

impl std::error::Error for SqiSignError {}
