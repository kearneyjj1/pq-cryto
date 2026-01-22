//! SLH-DSA parameter sets as defined in FIPS 205.
//!
//! This module defines parameter sets for SLH-DSA-SHAKE variants:
//! - SLH-DSA-SHAKE-128s/128f: NIST Level 1 (~128-bit security)
//! - SLH-DSA-SHAKE-192s/192f: NIST Level 3 (~192-bit security)
//! - SLH-DSA-SHAKE-256s/256f: NIST Level 5 (~256-bit security)
//!
//! The 's' variants have smaller signatures but slower signing.
//! The 'f' variants have faster signing but larger signatures.

/// Winternitz parameter (always 16 for FIPS 205).
pub const W: usize = 16;

/// Log base 2 of W.
pub const LOG_W: usize = 4;

/// Parameters for a specific SLH-DSA security level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Params {
    /// Security parameter n (hash output size in bytes): 16, 24, or 32
    pub n: usize,

    /// Total tree height h
    pub h: usize,

    /// Number of layers in the hypertree d
    pub d: usize,

    /// Height of each XMSS tree (h' = h/d)
    pub hp: usize,

    /// FORS tree height (a)
    pub a: usize,

    /// Number of FORS trees (k)
    pub k: usize,

    /// Length of WOTS+ signature in n-byte strings (len = len1 + len2)
    pub wots_len: usize,

    /// NIST security level (1, 3, or 5)
    pub security_level: usize,

    /// Variant: true for small signatures ('s'), false for fast signing ('f')
    pub small: bool,
}

impl Params {
    /// Returns the public key size in bytes.
    ///
    /// Public key = PK.seed (n bytes) + PK.root (n bytes) = 2n bytes
    pub const fn public_key_size(&self) -> usize {
        2 * self.n
    }

    /// Returns the secret key size in bytes.
    ///
    /// Secret key = SK.seed (n) + SK.prf (n) + PK.seed (n) + PK.root (n) = 4n bytes
    pub const fn secret_key_size(&self) -> usize {
        4 * self.n
    }

    /// Returns the signature size in bytes.
    ///
    /// Signature = R (n) + SIG_FORS (k*(a+1)*n) + SIG_HT ((h + d*wots_len)*n)
    pub const fn signature_size(&self) -> usize {
        // R: n bytes (randomizer)
        // FORS signature: k * (a + 1) * n bytes
        //   - k secret key values (k * n)
        //   - k authentication paths of a nodes each (k * a * n)
        // HT signature: (h + d * wots_len) * n bytes
        //   - d XMSS signatures, each with:
        //     - WOTS+ signature (wots_len * n)
        //     - Authentication path (hp * n)
        //   Total: d * (wots_len + hp) * n = d * wots_len * n + d * hp * n = d * wots_len * n + h * n
        let fors_size = self.k * (self.a + 1) * self.n;
        let ht_size = (self.h + self.d * self.wots_len) * self.n;
        self.n + fors_size + ht_size
    }

    /// Returns WOTS+ len1 parameter.
    ///
    /// len1 = ceil(8n / log2(w)) = ceil(8n / 4) = 2n for w=16
    pub const fn wots_len1(&self) -> usize {
        2 * self.n
    }

    /// Returns WOTS+ len2 parameter.
    ///
    /// len2 = floor(log2(len1 * (w-1)) / log2(w)) + 1
    pub const fn wots_len2(&self) -> usize {
        self.wots_len - self.wots_len1()
    }

    /// Returns the message digest length for FORS in bytes.
    ///
    /// md_len = ceil(k * a / 8)
    pub const fn md_len(&self) -> usize {
        (self.k * self.a + 7) / 8
    }

    /// Returns the number of bits for tree index.
    ///
    /// tree_bits = h - hp
    pub const fn tree_bits(&self) -> usize {
        self.h - self.hp
    }

    /// Returns the number of bits for leaf index.
    ///
    /// leaf_bits = hp
    pub const fn leaf_bits(&self) -> usize {
        self.hp
    }

    /// Returns the total index bits needed from H_msg.
    ///
    /// idx_len = ceil((h - hp + hp) / 8) = ceil(h / 8)
    pub const fn idx_len(&self) -> usize {
        (self.h + 7) / 8
    }
}

// ============================================================================
// SLH-DSA-SHAKE-128s Parameters (NIST Level 1, Small signatures)
// ============================================================================

/// SLH-DSA-SHAKE-128s parameters (~128-bit security, small signatures).
///
/// This parameter set prioritizes smaller signature sizes at the cost of
/// slower signing operations.
pub const SLH_DSA_SHAKE_128S: Params = Params {
    n: 16,
    h: 63,
    d: 7,
    hp: 9, // h' = 63/7 = 9
    a: 12,
    k: 14,
    wots_len: 35, // len1 = 32, len2 = 3
    security_level: 1,
    small: true,
};

// ============================================================================
// SLH-DSA-SHAKE-128f Parameters (NIST Level 1, Fast signing)
// ============================================================================

/// SLH-DSA-SHAKE-128f parameters (~128-bit security, fast signing).
///
/// This parameter set prioritizes faster signing at the cost of larger
/// signature sizes.
pub const SLH_DSA_SHAKE_128F: Params = Params {
    n: 16,
    h: 66,
    d: 22,
    hp: 3, // h' = 66/22 = 3
    a: 6,
    k: 33,
    wots_len: 35, // len1 = 32, len2 = 3
    security_level: 1,
    small: false,
};

// ============================================================================
// SLH-DSA-SHAKE-192s Parameters (NIST Level 3, Small signatures)
// ============================================================================

/// SLH-DSA-SHAKE-192s parameters (~192-bit security, small signatures).
pub const SLH_DSA_SHAKE_192S: Params = Params {
    n: 24,
    h: 63,
    d: 7,
    hp: 9, // h' = 63/7 = 9
    a: 14,
    k: 17,
    wots_len: 51, // len1 = 48, len2 = 3
    security_level: 3,
    small: true,
};

// ============================================================================
// SLH-DSA-SHAKE-192f Parameters (NIST Level 3, Fast signing)
// ============================================================================

/// SLH-DSA-SHAKE-192f parameters (~192-bit security, fast signing).
pub const SLH_DSA_SHAKE_192F: Params = Params {
    n: 24,
    h: 66,
    d: 22,
    hp: 3, // h' = 66/22 = 3
    a: 8,
    k: 33,
    wots_len: 51, // len1 = 48, len2 = 3
    security_level: 3,
    small: false,
};

// ============================================================================
// SLH-DSA-SHAKE-256s Parameters (NIST Level 5, Small signatures)
// ============================================================================

/// SLH-DSA-SHAKE-256s parameters (~256-bit security, small signatures).
pub const SLH_DSA_SHAKE_256S: Params = Params {
    n: 32,
    h: 64,
    d: 8,
    hp: 8, // h' = 64/8 = 8
    a: 14,
    k: 22,
    wots_len: 67, // len1 = 64, len2 = 3
    security_level: 5,
    small: true,
};

// ============================================================================
// SLH-DSA-SHAKE-256f Parameters (NIST Level 5, Fast signing)
// ============================================================================

/// SLH-DSA-SHAKE-256f parameters (~256-bit security, fast signing).
pub const SLH_DSA_SHAKE_256F: Params = Params {
    n: 32,
    h: 68,
    d: 17,
    hp: 4, // h' = 68/17 = 4
    a: 9,
    k: 35,
    wots_len: 67, // len1 = 64, len2 = 3
    security_level: 5,
    small: false,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hp_equals_h_div_d() {
        assert_eq!(SLH_DSA_SHAKE_128S.hp, SLH_DSA_SHAKE_128S.h / SLH_DSA_SHAKE_128S.d);
        assert_eq!(SLH_DSA_SHAKE_128F.hp, SLH_DSA_SHAKE_128F.h / SLH_DSA_SHAKE_128F.d);
        assert_eq!(SLH_DSA_SHAKE_192S.hp, SLH_DSA_SHAKE_192S.h / SLH_DSA_SHAKE_192S.d);
        assert_eq!(SLH_DSA_SHAKE_192F.hp, SLH_DSA_SHAKE_192F.h / SLH_DSA_SHAKE_192F.d);
        assert_eq!(SLH_DSA_SHAKE_256S.hp, SLH_DSA_SHAKE_256S.h / SLH_DSA_SHAKE_256S.d);
        assert_eq!(SLH_DSA_SHAKE_256F.hp, SLH_DSA_SHAKE_256F.h / SLH_DSA_SHAKE_256F.d);
    }

    #[test]
    fn test_wots_len_equals_len1_plus_len2() {
        for params in [
            SLH_DSA_SHAKE_128S,
            SLH_DSA_SHAKE_128F,
            SLH_DSA_SHAKE_192S,
            SLH_DSA_SHAKE_192F,
            SLH_DSA_SHAKE_256S,
            SLH_DSA_SHAKE_256F,
        ] {
            assert_eq!(params.wots_len, params.wots_len1() + params.wots_len2());
        }
    }

    #[test]
    fn test_wots_len1_calculation() {
        // len1 = ceil(8n / log2(w)) = 2n for w=16
        assert_eq!(SLH_DSA_SHAKE_128S.wots_len1(), 32); // 2 * 16
        assert_eq!(SLH_DSA_SHAKE_192S.wots_len1(), 48); // 2 * 24
        assert_eq!(SLH_DSA_SHAKE_256S.wots_len1(), 64); // 2 * 32
    }

    #[test]
    fn test_public_key_sizes() {
        // pk = 2n bytes
        assert_eq!(SLH_DSA_SHAKE_128S.public_key_size(), 32);
        assert_eq!(SLH_DSA_SHAKE_128F.public_key_size(), 32);
        assert_eq!(SLH_DSA_SHAKE_192S.public_key_size(), 48);
        assert_eq!(SLH_DSA_SHAKE_192F.public_key_size(), 48);
        assert_eq!(SLH_DSA_SHAKE_256S.public_key_size(), 64);
        assert_eq!(SLH_DSA_SHAKE_256F.public_key_size(), 64);
    }

    #[test]
    fn test_secret_key_sizes() {
        // sk = 4n bytes
        assert_eq!(SLH_DSA_SHAKE_128S.secret_key_size(), 64);
        assert_eq!(SLH_DSA_SHAKE_128F.secret_key_size(), 64);
        assert_eq!(SLH_DSA_SHAKE_192S.secret_key_size(), 96);
        assert_eq!(SLH_DSA_SHAKE_192F.secret_key_size(), 96);
        assert_eq!(SLH_DSA_SHAKE_256S.secret_key_size(), 128);
        assert_eq!(SLH_DSA_SHAKE_256F.secret_key_size(), 128);
    }

    #[test]
    fn test_signature_sizes() {
        // Expected sizes from FIPS 205 Table 1
        assert_eq!(SLH_DSA_SHAKE_128S.signature_size(), 7856);
        assert_eq!(SLH_DSA_SHAKE_128F.signature_size(), 17088);
        assert_eq!(SLH_DSA_SHAKE_192S.signature_size(), 16224);
        assert_eq!(SLH_DSA_SHAKE_192F.signature_size(), 35664);
        assert_eq!(SLH_DSA_SHAKE_256S.signature_size(), 29792);
        assert_eq!(SLH_DSA_SHAKE_256F.signature_size(), 49856);
    }

    #[test]
    fn test_128s_parameters() {
        let p = SLH_DSA_SHAKE_128S;
        assert_eq!(p.n, 16);
        assert_eq!(p.h, 63);
        assert_eq!(p.d, 7);
        assert_eq!(p.hp, 9);
        assert_eq!(p.a, 12);
        assert_eq!(p.k, 14);
        assert_eq!(p.wots_len, 35);
        assert_eq!(p.security_level, 1);
        assert!(p.small);
    }

    #[test]
    fn test_128f_parameters() {
        let p = SLH_DSA_SHAKE_128F;
        assert_eq!(p.n, 16);
        assert_eq!(p.h, 66);
        assert_eq!(p.d, 22);
        assert_eq!(p.hp, 3);
        assert_eq!(p.a, 6);
        assert_eq!(p.k, 33);
        assert_eq!(p.wots_len, 35);
        assert_eq!(p.security_level, 1);
        assert!(!p.small);
    }

    #[test]
    fn test_192s_parameters() {
        let p = SLH_DSA_SHAKE_192S;
        assert_eq!(p.n, 24);
        assert_eq!(p.h, 63);
        assert_eq!(p.d, 7);
        assert_eq!(p.hp, 9);
        assert_eq!(p.a, 14);
        assert_eq!(p.k, 17);
        assert_eq!(p.wots_len, 51);
        assert_eq!(p.security_level, 3);
    }

    #[test]
    fn test_192f_parameters() {
        let p = SLH_DSA_SHAKE_192F;
        assert_eq!(p.n, 24);
        assert_eq!(p.h, 66);
        assert_eq!(p.d, 22);
        assert_eq!(p.hp, 3);
        assert_eq!(p.a, 8);
        assert_eq!(p.k, 33);
        assert_eq!(p.wots_len, 51);
        assert_eq!(p.security_level, 3);
    }

    #[test]
    fn test_256s_parameters() {
        let p = SLH_DSA_SHAKE_256S;
        assert_eq!(p.n, 32);
        assert_eq!(p.h, 64);
        assert_eq!(p.d, 8);
        assert_eq!(p.hp, 8);
        assert_eq!(p.a, 14);
        assert_eq!(p.k, 22);
        assert_eq!(p.wots_len, 67);
        assert_eq!(p.security_level, 5);
    }

    #[test]
    fn test_256f_parameters() {
        let p = SLH_DSA_SHAKE_256F;
        assert_eq!(p.n, 32);
        assert_eq!(p.h, 68);
        assert_eq!(p.d, 17);
        assert_eq!(p.hp, 4);
        assert_eq!(p.a, 9);
        assert_eq!(p.k, 35);
        assert_eq!(p.wots_len, 67);
        assert_eq!(p.security_level, 5);
    }
}
