//! Regression tests for audit finding CRIT-01 (2026-05-16): strict-EUF-CMA
//! malleability via non-canonical signature encoding.
//!
//! The Falcon spec (round-3 §3.11.2 `comp_decode`, mirrored by the PQClean
//! reference) requires the decoder to reject any signature with:
//!   (a) non-zero pad bits in the final byte of the compressed bitstream, or
//!   (b) extra bytes beyond the natural end of the encoded stream (in the
//!       non-padded format used by FIPS 206 and the legacy FSC1 format).
//!
//! Before the fix, `decode_golomb_rice` and `decode_falcon_comp` ignored both,
//! so every valid signature had many distinct wire encodings that all verified.
//! For any setting that relies on signature identity (blockchains, audit logs,
//! payment confirmations), that breaks strict-EUF-CMA.
//!
//! These tests assert the rejection at the decoder boundary — they should
//! fail on the pre-fix decoder and pass on the post-fix decoder.

use pqsigs_fn_dsa::packing::{
    decode_signature, encode_signature, encode_signature_compressed,
};
use pqsigs_fn_dsa::{keygen_512, sign, verify, FALCON_512};
use rand::rngs::StdRng;
use rand::SeedableRng;

/// Appending a non-zero byte to a FIPS 206 signature must be rejected at decode time.
#[test]
fn crit01_nonzero_trailing_byte_rejected_fips206() {
    let mut rng = StdRng::seed_from_u64(0xc71701_f175_2026);
    let kp = keygen_512(&mut rng).expect("keygen");
    let msg = b"crit-01 trailing-byte fips206";
    let sig = sign(&mut rng, &kp.sk, msg).expect("sign");
    verify(&kp.pk, msg, &sig).expect("baseline verify");

    let bytes = encode_signature(&sig, FALCON_512.n, FALCON_512.sig_bytes_max)
        .expect("encode");
    decode_signature(&bytes).expect("baseline decode of canonical encoding");

    for trailing in [0x01u8, 0x55, 0xAA, 0xFF] {
        let mut mutated = bytes.clone();
        mutated.push(trailing);
        assert!(
            decode_signature(&mutated).is_err(),
            "FIPS 206 decode must reject signature with appended byte {trailing:#04x}"
        );
    }
}

/// Appending a *zero* byte must also be rejected — FIPS 206 is non-padded,
/// so any extra byte (even zero) is a strict-EUF-CMA malleability vector.
#[test]
fn crit01_zero_trailing_byte_rejected_fips206() {
    let mut rng = StdRng::seed_from_u64(0xc71701_2e20_2026);
    let kp = keygen_512(&mut rng).expect("keygen");
    let msg = b"crit-01 zero-trailing-byte fips206";
    let sig = sign(&mut rng, &kp.sk, msg).expect("sign");

    let mut bytes = encode_signature(&sig, FALCON_512.n, FALCON_512.sig_bytes_max)
        .expect("encode");
    bytes.push(0x00);

    assert!(
        decode_signature(&bytes).is_err(),
        "FIPS 206 decode must reject signature with appended zero byte"
    );
}

/// Same property for the legacy FSC1 compressed format — it is also non-padded.
#[test]
fn crit01_trailing_byte_rejected_legacy_compressed() {
    let mut rng = StdRng::seed_from_u64(0xc71701_f5c1_2026);
    let kp = keygen_512(&mut rng).expect("keygen");
    let msg = b"crit-01 trailing-byte fsc1";
    let sig = sign(&mut rng, &kp.sk, msg).expect("sign");

    let bytes = encode_signature_compressed(&sig, FALCON_512.n);
    decode_signature(&bytes).expect("baseline decode of FSC1 canonical encoding");

    for trailing in [0x00u8, 0x01, 0xAA, 0xFF] {
        let mut mutated = bytes.clone();
        mutated.push(trailing);
        assert!(
            decode_signature(&mutated).is_err(),
            "FSC1 decode must reject signature with appended byte {trailing:#04x}"
        );
    }
}

/// Flipping a pad bit in the final byte of the compressed bitstream must be
/// rejected at decode time. The encoder always emits the high pad bits as zero,
/// so setting any high bit on the last byte produces a non-canonical encoding
/// whenever that bit was actually a pad bit (vs. an actual zero data bit).
///
/// Without the CRIT-01 fix, non-canonical encodings with pad-bit mutations
/// silently round-trip and verify, breaking strict-EUF-CMA. With the fix,
/// at least one of the eight bit positions in the last byte will trigger a
/// decoder rejection (and the original sig is recoverable from the mutated
/// one — proving a many-to-one map existed pre-fix).
///
/// We sweep multiple message instances to find ones whose last byte contains
/// pad bits (i.e., bit_count > 0 in the encoder's terminal flush). If every
/// sig encoded ends exactly on a byte boundary, the test logs `inconclusive`
/// rather than spuriously passing or failing.
#[test]
fn crit01_pad_bit_in_final_byte_rejected_fips206() {
    let mut rng = StdRng::seed_from_u64(0xc71701_dab1_2026);
    let kp = keygen_512(&mut rng).expect("keygen");

    let mut observed_pad_rejection = false;
    let mut probes_tried = 0usize;

    for trial in 0..40 {
        let msg = format!("crit-01 padbit trial {trial}").into_bytes();
        let sig = sign(&mut rng, &kp.sk, &msg).expect("sign");
        let bytes = encode_signature(&sig, FALCON_512.n, FALCON_512.sig_bytes_max)
            .expect("encode");

        decode_signature(&bytes).expect("baseline decode");
        verify(&kp.pk, &msg, &sig).expect("baseline verify");

        let last_idx = bytes.len() - 1;
        let last = bytes[last_idx];

        for bit_pos in 0..8u8 {
            if (last >> bit_pos) & 1 != 0 {
                continue;
            }
            let mut mutated = bytes.clone();
            mutated[last_idx] ^= 1 << bit_pos;
            probes_tried += 1;

            match decode_signature(&mutated) {
                Err(_) => {
                    observed_pad_rejection = true;
                }
                Ok((decoded, _)) => {
                    // If decode accepts the mutation, the only safety net is
                    // verify catching it. The decoded sig MUST NOT verify under
                    // the original key/message (strict-EUF-CMA).
                    assert!(
                        verify(&kp.pk, &msg, &decoded).is_err(),
                        "CRIT-01 regression: bit {bit_pos} flipped in last byte (trial {trial}) decoded AND verified — strict-EUF-CMA broken"
                    );
                }
            }
        }
        if observed_pad_rejection && probes_tried > 50 {
            break;
        }
    }

    assert!(
        observed_pad_rejection,
        "Inconclusive: no pad-bit rejection observed in {probes_tried} probes — either the encoder produced byte-aligned outputs every time (very unlikely) or the canonical check is not engaged"
    );
}

// NOTE: While developing these CRIT-01 tests we observed that the FIPS 206
// (header byte 0x39) and legacy FSC1 ("FSC1" + n) formats share an identical
// inner layout: 40-byte nonce followed by the same LSB-first Golomb-Rice
// bitstream. As a result, the body of a FIPS 206 signature re-wrapped with
// an FSC1 header decodes and verifies under the original key. That is a
// distinct strict-EUF-CMA malleability vector from CRIT-01 (encoding-internal
// pad bits and trailing bytes) — it is *cross-format* malleability and is
// best addressed by domain-separating the formats at the bitstream level
// (e.g. mixing the format tag into the SHAKE absorb in hash_to_point, or
// removing one of the format adapters from the public API). It is tracked
// separately and intentionally NOT asserted here so that the CRIT-01
// regression suite has clear semantic scope.
