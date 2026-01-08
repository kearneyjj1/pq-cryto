//! Number Theoretic Transform for ML-DSA.
//!
//! The NTT is used for efficient polynomial multiplication in the ring
//! R_q = Z_q[X]/(X^256 + 1). Using NTT, multiplication becomes pointwise
//! multiplication of the coefficient vectors, reducing complexity from
//! O(n^2) to O(n log n).
//!
//! This is a reference implementation prioritizing clarity over speed.
//! It uses simple modular arithmetic without Montgomery optimization.

use crate::params::{N, Q, ZETA};
use crate::reduce::freeze;

/// The multiplicative inverse of N (256) modulo Q (8380417).
///
/// Used to scale the result of the inverse NTT.
/// Computed as: 256^(-1) mod 8380417 = 8347681
/// Verification: 256 * 8347681 = 2136942336 = 255 * 8380417 + 1 ≡ 1 (mod Q)
const N_INV: i32 = 8347681;

/// Computes a * b mod Q with proper handling of negative values.
#[inline]
fn mulmod(a: i32, b: i32) -> i32 {
    let prod = (a as i64) * (b as i64);
    let mut r = (prod % Q as i64) as i32;
    if r < 0 {
        r += Q;
    }
    r
}

/// Computes base^exp mod Q using square-and-multiply.
fn powmod(base: i32, mut exp: u32) -> i32 {
    let mut result = 1i64;
    let mut b = ((base as i64 % Q as i64) + Q as i64) % Q as i64;
    let q = Q as i64;

    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * b) % q;
        }
        b = (b * b) % q;
        exp >>= 1;
    }

    result as i32
}

/// Precompute all 256 roots of unity needed for NTT.
/// zetas[i] = ζ^(brv(i)) where brv is bit-reversal permutation.
fn precompute_zetas() -> [i32; 256] {
    let mut zetas = [0i32; 256];

    for i in 0..256 {
        // For negacyclic NTT: ζ^(2*BitRev_7(i) + 1) for indices used in layers
        // But we'll compute ζ^(BitRev_8(i)) for full table
        let br = bitrev8(i);
        zetas[i] = powmod(ZETA, br as u32);
    }

    zetas
}

/// Bit-reverse an 8-bit integer.
fn bitrev8(x: usize) -> usize {
    let mut result = 0;
    let mut val = x;
    for _ in 0..8 {
        result = (result << 1) | (val & 1);
        val >>= 1;
    }
    result
}

/// Forward NTT (Number Theoretic Transform).
///
/// Transforms a polynomial from coefficient representation to NTT domain.
/// This implementation uses Cooley-Tukey decimation-in-time.
pub fn ntt(a: &mut [i32; N]) {
    let zetas = precompute_zetas();
    let mut k = 0usize;
    let mut len = 128;

    while len >= 1 {
        let mut start = 0;
        while start < N {
            k += 1;
            let zeta = zetas[k];
            for j in start..start + len {
                let t = mulmod(zeta, a[j + len]);
                a[j + len] = freeze(a[j] - t);
                a[j] = freeze(a[j] + t);
            }
            start += 2 * len;
        }
        len >>= 1;
    }
}

/// Inverse NTT.
///
/// Transforms a polynomial from NTT domain back to coefficient representation.
/// This implementation uses Gentleman-Sande decimation-in-frequency.
pub fn inv_ntt(a: &mut [i32; N]) {
    let zetas = precompute_zetas();
    let mut k = 256usize;
    let mut len = 1;

    while len < N {
        let mut start = 0;
        while start < N {
            k -= 1;
            // For inverse, we need -zeta
            let zeta = Q - zetas[k];
            for j in start..start + len {
                let t = a[j];
                a[j] = freeze(t + a[j + len]);
                a[j + len] = mulmod(zeta, freeze(t - a[j + len]));
            }
            start += 2 * len;
        }
        len <<= 1;
    }

    // Scale by n^(-1) mod q
    for coeff in a.iter_mut() {
        *coeff = mulmod(*coeff, N_INV);
    }
}

/// Pointwise multiplication of two polynomials in NTT domain.
pub fn pointwise_mul(c: &mut [i32; N], a: &[i32; N], b: &[i32; N]) {
    for i in 0..N {
        c[i] = mulmod(a[i], b[i]);
    }
}

/// Convert coefficients to fully reduced form in [0, Q).
pub fn reduce_coeffs(coeffs: &mut [i32; N]) {
    for coeff in coeffs.iter_mut() {
        *coeff = freeze(*coeff);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_powmod() {
        // ZETA^256 should equal -1 mod Q
        assert_eq!(powmod(ZETA, 256), Q - 1);
        // ZETA^512 should equal 1 mod Q
        assert_eq!(powmod(ZETA, 512), 1);
    }

    #[test]
    fn test_n_inv() {
        // 256 * 8347681 mod Q should equal 1
        let prod = mulmod(256, 8347681);
        assert_eq!(prod, 1);
    }

    #[test]
    fn test_ntt_inv_ntt_roundtrip() {
        let mut a = [0i32; N];
        for i in 0..N {
            a[i] = (i as i32 * 17) % 1000;
        }
        let original = a;

        ntt(&mut a);
        inv_ntt(&mut a);

        for i in 0..N {
            let recovered = freeze(a[i]);
            let expected = freeze(original[i]);
            assert_eq!(
                recovered, expected,
                "Mismatch at {}: got {}, expected {}",
                i, recovered, expected
            );
        }
    }

    #[test]
    fn test_ntt_zero() {
        let mut coeffs = [0i32; N];
        ntt(&mut coeffs);
        for &c in coeffs.iter() {
            assert_eq!(c, 0);
        }
    }

    #[test]
    fn test_pointwise_mul_identity() {
        let mut a = [0i32; N];
        for i in 0..N {
            a[i] = (i as i32 * 7) % 1000;
        }
        let original = a;

        let mut one = [0i32; N];
        one[0] = 1;

        ntt(&mut a);
        ntt(&mut one);

        let mut result = [0i32; N];
        pointwise_mul(&mut result, &a, &one);

        inv_ntt(&mut result);

        for i in 0..N {
            let r = freeze(result[i]);
            let e = freeze(original[i]);
            assert_eq!(r, e, "Mismatch at {}: got {}, expected {}", i, r, e);
        }
    }

    #[test]
    fn test_polynomial_multiplication() {
        // Test (1 + X) * (1 + X) = 1 + 2X + X^2
        let mut a = [0i32; N];
        let mut b = [0i32; N];
        a[0] = 1;
        a[1] = 1;
        b[0] = 1;
        b[1] = 1;

        ntt(&mut a);
        ntt(&mut b);

        let mut c = [0i32; N];
        pointwise_mul(&mut c, &a, &b);

        inv_ntt(&mut c);

        assert_eq!(freeze(c[0]), 1, "c[0] should be 1");
        assert_eq!(freeze(c[1]), 2, "c[1] should be 2");
        assert_eq!(freeze(c[2]), 1, "c[2] should be 1");
        for i in 3..N {
            assert_eq!(freeze(c[i]), 0, "c[{}] should be 0", i);
        }
    }
}
