//! Error types for FN-DSA (FALCON).

use std::fmt;

/// Specific reasons why signature verification failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerificationFailure {
    /// The signature norm exceeds the allowed bound.
    NormBoundExceeded,
    /// The recovered hash doesn't match expected value.
    HashMismatch,
    /// The signature polynomial has invalid coefficients.
    InvalidCoefficients,
    /// Generic verification failure (unspecified reason).
    Generic,
}

impl fmt::Display for VerificationFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VerificationFailure::NormBoundExceeded => write!(f, "signature norm exceeded"),
            VerificationFailure::HashMismatch => write!(f, "hash mismatch"),
            VerificationFailure::InvalidCoefficients => write!(f, "invalid coefficients"),
            VerificationFailure::Generic => write!(f, "verification failed"),
        }
    }
}

/// Errors that can occur during FN-DSA operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FnDsaError {
    /// Signature verification failed with a specific reason.
    VerificationFailed(VerificationFailure),

    /// The signature failed verification (legacy variant).
    InvalidSignature,

    /// Signing failed after the specified number of attempts.
    SigningFailed {
        /// Number of attempts made before failure.
        attempts: u32,
    },

    /// The provided key is invalid.
    InvalidKey {
        /// Description of why the key is invalid.
        reason: &'static str,
    },

    /// The parameters are invalid.
    InvalidParams {
        /// Description of why the parameters are invalid.
        reason: &'static str,
    },

    /// Invalid input was provided.
    InvalidInput {
        /// The field that was invalid.
        field: &'static str,
        /// Description of why the input is invalid.
        reason: &'static str,
    },

    /// Error decoding data.
    DecodingError {
        /// Context about what was being decoded.
        context: &'static str,
    },

    /// Key generation failed.
    KeygenFailed {
        /// Description of why keygen failed.
        reason: &'static str,
    },

    /// Gaussian sampling failed.
    SamplingFailed,

    /// NTRU solve failed to find a solution.
    NtruSolveFailed,
}

impl fmt::Display for FnDsaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FnDsaError::VerificationFailed(reason) => {
                write!(f, "signature verification failed: {}", reason)
            }
            FnDsaError::InvalidSignature => write!(f, "invalid signature"),
            FnDsaError::SigningFailed { attempts } => {
                write!(f, "signing failed after {} attempts", attempts)
            }
            FnDsaError::InvalidKey { reason } => write!(f, "invalid key: {}", reason),
            FnDsaError::InvalidParams { reason } => write!(f, "invalid parameters: {}", reason),
            FnDsaError::InvalidInput { field, reason } => {
                write!(f, "invalid {}: {}", field, reason)
            }
            FnDsaError::DecodingError { context } => write!(f, "decoding error: {}", context),
            FnDsaError::KeygenFailed { reason } => write!(f, "key generation failed: {}", reason),
            FnDsaError::SamplingFailed => write!(f, "Gaussian sampling failed"),
            FnDsaError::NtruSolveFailed => write!(f, "NTRU solve failed to find solution"),
        }
    }
}

impl std::error::Error for FnDsaError {}

/// Result type alias for FN-DSA operations.
pub type Result<T> = std::result::Result<T, FnDsaError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        assert_eq!(
            format!("{}", FnDsaError::InvalidSignature),
            "invalid signature"
        );
        assert_eq!(
            format!("{}", FnDsaError::SigningFailed { attempts: 5 }),
            "signing failed after 5 attempts"
        );
        assert_eq!(
            format!("{}", FnDsaError::InvalidKey { reason: "corrupt" }),
            "invalid key: corrupt"
        );
    }

    #[test]
    fn test_error_is_error_trait() {
        let err: Box<dyn std::error::Error> = Box::new(FnDsaError::InvalidSignature);
        assert!(err.to_string().contains("invalid"));
    }
}
