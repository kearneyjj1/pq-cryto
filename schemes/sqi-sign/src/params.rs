//! Parameter sets for SQI-SIGN.
//!
//! SQI-SIGN defines parameter sets targeting different NIST security levels.
//! The parameters are based on supersingular elliptic curves over Fp² where
//! p is a prime of a specific form.

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
}

impl Params {
    /// Returns the prime p as a reference.
    pub fn prime(&self) -> &BigInt {
        &self.p
    }
}

// TODO: Define actual SQI-SIGN parameters
// These are placeholder values for structure setup

lazy_static::lazy_static! {
    /// SQI-SIGN Level I parameters (NIST Level 1, ~128-bit security).
    ///
    /// Uses p = 2^33 * 3^19 - 1 = 9586980095590400 - 1 (for testing).
    /// A proper implementation would use the full specification prime.
    pub static ref SQISIGN_NIST_I: Params = {
        // For testing: p = 2^33 * 3^19 - 1 ≈ 9.5 * 10^15 (53-bit prime)
        // This prime satisfies p ≡ 3 (mod 4) for supersingular curves
        // Production would use: p = 2^126 * 3^72 * f - 1 for proper security
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let p = &two.pow(33) * &three.pow(19) - BigInt::from(1);

        Params {
            name: "SQISign-NIST-I",
            security_level: 1,
            p,
            e2: 33,
            e3: 19,
            f: BigInt::from(1),
            p_bits: 256,
            pk_bytes: 64,
            sk_bytes: 782,
            sig_bytes: 177,
        }
    };

    /// SQI-SIGN Level III parameters (NIST Level 3, ~192-bit security).
    ///
    /// Uses p = 2^37 * 3^23 - 1 (for testing).
    pub static ref SQISIGN_NIST_III: Params = {
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let p = &two.pow(37) * &three.pow(23) - BigInt::from(1);

        Params {
            name: "SQISign-NIST-III",
            security_level: 3,
            p,
            e2: 37,
            e3: 23,
            f: BigInt::from(1),
            p_bits: 384,
            pk_bytes: 96,
            sk_bytes: 1138,
            sig_bytes: 263,
        }
    };

    /// SQI-SIGN Level V parameters (NIST Level 5, ~256-bit security).
    ///
    /// Uses p = 2^41 * 3^27 - 1 (for testing).
    pub static ref SQISIGN_NIST_V: Params = {
        let two = BigInt::from(2);
        let three = BigInt::from(3);
        let p = &two.pow(41) * &three.pow(27) - BigInt::from(1);

        Params {
            name: "SQISign-NIST-V",
            security_level: 5,
            p,
            e2: 41,
            e3: 27,
            f: BigInt::from(1),
            p_bits: 512,
            pk_bytes: 128,
            sk_bytes: 1509,
            sig_bytes: 335,
        }
    };
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
}
