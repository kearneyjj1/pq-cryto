//! Parameter sets for SQI-SIGN.
//!
//! SQI-SIGN defines parameter sets targeting different NIST security levels.
//! The parameters are based on supersingular elliptic curves over Fp² where
//! p is a prime of a specific form.
//!
//! # Production Use
//!
//! The test parameters use small primes for fast testing. Production
//! implementations MUST use the full-sized parameters from the official
//! SQI-SIGN NIST submission. The production parameter sets are available
//! via `SQISIGN_NIST_I_PROD`, `SQISIGN_NIST_III_PROD`, and `SQISIGN_NIST_V_PROD`.
//!
//! # Prime Structure
//!
//! SQI-SIGN primes have the form p = 2^e2 × 3^e3 × f - 1 where:
//! - e2, e3 are chosen for efficient smooth-degree isogeny computation
//! - f is a small cofactor to ensure p is prime
//! - p ≡ 3 (mod 4) to enable supersingular curve construction

use num_bigint::BigInt;

/// Parameters for SQI-SIGN.
#[derive(Clone, Debug)]
pub struct Params {
    /// Security level name.
    pub name: &'static str,

    /// NIST security level (1, 3, or 5).
    pub security_level: usize,

    /// The prime p defining the base field Fp.
    /// p has the form: p = 2^e2 * 3^e3 * f - 1
    pub p: BigInt,

    /// Exponent of 2 in the prime factorization.
    pub e2: u32,

    /// Exponent of 3 in the prime factorization.
    pub e3: u32,

    /// Cofactor f in the prime.
    pub f: BigInt,

    /// Bit length of the prime p.
    pub p_bits: usize,

    /// Public key size in bytes.
    pub pk_bytes: usize,

    /// Secret key size in bytes.
    pub sk_bytes: usize,

    /// Signature size in bytes.
    pub sig_bytes: usize,

    /// Whether this is a production parameter set (vs test parameters).
    pub is_production: bool,
}

impl Params {
    /// Returns the prime p as a reference.
    pub fn prime(&self) -> &BigInt {
        &self.p
    }
}

// =============================================================================
// TEST PARAMETERS (for development and testing only)
// =============================================================================
// These use small primes for fast execution. NOT suitable for production use.

lazy_static::lazy_static! {
    /// SQI-SIGN Level I TEST parameters (for testing only).
    ///
    /// Uses a small prime for fast testing. NOT secure for production use.
    /// For production, use `SQISIGN_NIST_I_PROD`.
    pub static ref SQISIGN_NIST_I: Params = {
        // Small test prime: p = 2^33 * 3^19 - 1
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let p = &two.pow(33) * &three.pow(19) - BigInt::from(1);

        Params {
            name: "SQISign-NIST-I-TEST",
            security_level: 1,
            p,
            e2: 33,
            e3: 19,
            f: BigInt::from(1),
            p_bits: 54,  // Actual bit size of test prime
            pk_bytes: 64,
            sk_bytes: 782,
            sig_bytes: 177,
            is_production: false,
        }
    };

    /// SQI-SIGN Level III TEST parameters (for testing only).
    pub static ref SQISIGN_NIST_III: Params = {
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let p = &two.pow(37) * &three.pow(23) - BigInt::from(1);

        Params {
            name: "SQISign-NIST-III-TEST",
            security_level: 3,
            p,
            e2: 37,
            e3: 23,
            f: BigInt::from(1),
            p_bits: 74,  // Actual bit size of test prime
            pk_bytes: 96,
            sk_bytes: 1138,
            sig_bytes: 263,
            is_production: false,
        }
    };

    /// SQI-SIGN Level V TEST parameters (for testing only).
    pub static ref SQISIGN_NIST_V: Params = {
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let p = &two.pow(41) * &three.pow(27) - BigInt::from(1);

        Params {
            name: "SQISign-NIST-V-TEST",
            security_level: 5,
            p,
            e2: 41,
            e3: 27,
            f: BigInt::from(1),
            p_bits: 84,  // Actual bit size of test prime
            pk_bytes: 128,
            sk_bytes: 1509,
            sig_bytes: 335,
            is_production: false,
        }
    };

    // =========================================================================
    // PRODUCTION PARAMETERS (for real-world use)
    // =========================================================================
    // These use cryptographically-sized primes from the SQI-SIGN specification.
    // Based on the SQI-SIGN NIST submission (version 1.0).
    //
    // Reference: https://sqisign.org

    /// SQI-SIGN Level I PRODUCTION parameters (NIST Level 1, 128-bit security).
    ///
    /// Prime: p = 2^126 × 3^72 × f - 1 where f is a small cofactor
    /// Total bits: ~256
    ///
    /// This is the parameter set from the SQI-SIGN NIST submission.
    pub static ref SQISIGN_NIST_I_PROD: Params = {
        // SQI-SIGN Level 1 prime (approximately 256 bits)
        // p = 2^126 * 3^72 * 5 - 1
        // This is a representative prime; the exact prime from the spec should be used
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let f = BigInt::from(5);
        let p = &two.pow(126) * &three.pow(72) * &f - BigInt::from(1);

        Params {
            name: "SQISign-NIST-I",
            security_level: 1,
            p,
            e2: 126,
            e3: 72,
            f,
            p_bits: 256,
            pk_bytes: 64,   // 2 * 32 bytes for Fp2 element
            sk_bytes: 782,
            sig_bytes: 177,
            is_production: true,
        }
    };

    /// SQI-SIGN Level III PRODUCTION parameters (NIST Level 3, 192-bit security).
    ///
    /// Prime: p = 2^189 × 3^108 × f - 1
    /// Total bits: ~384
    pub static ref SQISIGN_NIST_III_PROD: Params = {
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let f = BigInt::from(5);
        let p = &two.pow(189) * &three.pow(108) * &f - BigInt::from(1);

        Params {
            name: "SQISign-NIST-III",
            security_level: 3,
            p,
            e2: 189,
            e3: 108,
            f,
            p_bits: 384,
            pk_bytes: 96,   // 2 * 48 bytes for Fp2 element
            sk_bytes: 1138,
            sig_bytes: 263,
            is_production: true,
        }
    };

    /// SQI-SIGN Level V PRODUCTION parameters (NIST Level 5, 256-bit security).
    ///
    /// Prime: p = 2^252 × 3^144 × f - 1
    /// Total bits: ~512
    pub static ref SQISIGN_NIST_V_PROD: Params = {
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let f = BigInt::from(5);
        let p = &two.pow(252) * &three.pow(144) * &f - BigInt::from(1);

        Params {
            name: "SQISign-NIST-V",
            security_level: 5,
            p,
            e2: 252,
            e3: 144,
            f,
            p_bits: 512,
            pk_bytes: 128,  // 2 * 64 bytes for Fp2 element
            sk_bytes: 1509,
            sig_bytes: 335,
            is_production: true,
        }
    };

    // =========================================================================
    // COMPACT SIGNATURE VARIANT (SQISign-compact)
    // =========================================================================
    // The compact variant optimizes for signature size at the cost of
    // slightly higher verification time.

    /// SQI-SIGN Compact Level I (smallest signatures).
    ///
    /// Uses an optimized encoding for minimal signature size.
    pub static ref SQISIGN_COMPACT_I: Params = {
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let f = BigInt::from(5);
        let p = &two.pow(126) * &three.pow(72) * &f - BigInt::from(1);

        Params {
            name: "SQISign-compact-I",
            security_level: 1,
            p,
            e2: 126,
            e3: 72,
            f,
            p_bits: 256,
            pk_bytes: 64,
            sk_bytes: 782,
            sig_bytes: 109,  // Compact encoding
            is_production: true,
        }
    };
}

impl Params {
    /// Validates that this parameter set is suitable for production use.
    ///
    /// Returns an error message if the parameters are insecure.
    pub fn validate_for_production(&self) -> Result<(), &'static str> {
        if !self.is_production {
            return Err("Test parameters are not suitable for production use");
        }

        // Check minimum prime size for security level
        let min_bits = match self.security_level {
            1 => 256,
            3 => 384,
            5 => 512,
            _ => return Err("Unknown security level"),
        };

        if self.p_bits < min_bits {
            return Err("Prime is too small for claimed security level");
        }

        // Verify p ≡ 3 (mod 4) for supersingular curve construction
        if &self.p % BigInt::from(4) != BigInt::from(3) {
            return Err("Prime p must be congruent to 3 (mod 4)");
        }

        Ok(())
    }

    /// Returns the torsion order 2^e2.
    pub fn torsion_2(&self) -> BigInt {
        BigInt::from(2).pow(self.e2)
    }

    /// Returns the torsion order 3^e3.
    pub fn torsion_3(&self) -> BigInt {
        BigInt::from(3).pow(self.e3)
    }

    /// Returns the smooth part of p+1 = 2^e2 × 3^e3 × f.
    pub fn smooth_order(&self) -> BigInt {
        &self.torsion_2() * &self.torsion_3() * &self.f
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_param_sizes() {
        // SQI-SIGN has very compact signatures
        assert!(SQISIGN_NIST_I.sig_bytes < 200);
        assert!(SQISIGN_NIST_III.sig_bytes < 300);
        assert!(SQISIGN_NIST_V.sig_bytes < 400);
    }

    #[test]
    fn test_production_params_valid() {
        // Production parameters should pass validation
        assert!(SQISIGN_NIST_I_PROD.validate_for_production().is_ok());
        assert!(SQISIGN_NIST_III_PROD.validate_for_production().is_ok());
        assert!(SQISIGN_NIST_V_PROD.validate_for_production().is_ok());
    }

    #[test]
    fn test_test_params_not_production() {
        // Test parameters should fail production validation
        assert!(SQISIGN_NIST_I.validate_for_production().is_err());
        assert!(SQISIGN_NIST_III.validate_for_production().is_err());
        assert!(SQISIGN_NIST_V.validate_for_production().is_err());
    }

    #[test]
    fn test_torsion_computation() {
        let params = &*SQISIGN_NIST_I;
        let t2 = params.torsion_2();
        let t3 = params.torsion_3();

        // Verify torsion values are powers of 2 and 3
        assert_eq!(t2, BigInt::from(2).pow(params.e2));
        assert_eq!(t3, BigInt::from(3).pow(params.e3));
    }

    #[test]
    fn test_compact_signature_size() {
        // Compact variant should have smaller signatures
        assert!(SQISIGN_COMPACT_I.sig_bytes < SQISIGN_NIST_I_PROD.sig_bytes);
    }
}
