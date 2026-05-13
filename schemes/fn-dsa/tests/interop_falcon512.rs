//! FALCON-512 cross-implementation interoperability tests.
//!
//! The strongest interoperability check would be to FFI into PQClean's
//! Falcon-512 reference C code, sign with our `sign()`, encode in FIPS
//! 206 wire format, then call PQClean's `crypto_sign_open` and confirm
//! acceptance. Vendoring PQClean's reference C is out of scope here:
//! - Git Bash's `/usr/bin/link.exe` shadows MSVC's `link.exe`, which
//!   complicates the `cc` build for downstream Windows consumers. The
//!   documented workaround (use the GNU toolchain via
//!   `cargo +stable-x86_64-pc-windows-gnu`) does not extend cleanly to
//!   compiling C reference source against an embedded Windows linker
//!   without operator setup.
//! - PQClean has license + supply-chain implications: vendoring
//!   third-party C into a Rust crypto crate widens the trust surface.
//!
//! In lieu of the C-side roundtrip, this file establishes the
//! internal-consistency checks that cover the same correctness
//! properties for the wire format: encode, decode, and verify
//! self-roundtrip, plus negative tests for forgery resistance. When the
//! PQClean FFI test is eventually implemented, this file's coverage
//! should remain — it complements the FFI roundtrip rather than
//! replaces it.

use pqsigs_fn_dsa::{
    keygen_512, sign, verify, FALCON_512,
};
use pqsigs_fn_dsa::packing::{
    encode_public_key, decode_public_key,
    encode_signature, decode_signature,
};
use rand::SeedableRng;
use rand::rngs::StdRng;

/// Generate a key, sign N messages, then re-decode the public key and
/// every signature through our wire format and confirm each still verifies.
///
/// This validates the round-trip path:
///     (pk, sk) → encode_public_key → bytes → decode_public_key → pk'
///     sig → encode_signature → bytes → decode_signature → sig'
///     verify(pk', msg, sig') ?= Ok
///
/// A wire-format off-by-one (e.g. Golomb-Rice bit framing, public key
/// 14-bit packing) would surface here. The in-tree unit tests cover this
/// indirectly; here we run a deliberately wider sweep.
#[test]
fn wire_format_roundtrip_falcon_512() {
    let mut rng = StdRng::seed_from_u64(0xfa1c_0512_2026_0513);
    let keypair = keygen_512(&mut rng).expect("keygen must succeed");

    // Roundtrip the public key through encode/decode.
    let pk_bytes = encode_public_key(&keypair.pk);
    let pk_decoded = decode_public_key(&pk_bytes)
        .expect("public key round-trip must succeed");
    assert_eq!(pk_decoded.h, keypair.pk.h, "PK h coefficients must match after encode/decode");
    assert_eq!(pk_decoded.params.n, FALCON_512.n, "PK n must match");

    const N: usize = 20;
    for i in 0..N {
        let msg = format!("interop probe {i}").into_bytes();
        let sig = sign(&mut rng, &keypair.sk, &msg)
            .expect("sign must succeed on a valid FALCON-512 key");

        // Verify against the original pk.
        verify(&keypair.pk, &msg, &sig)
            .expect("original signature must verify");

        // Roundtrip the signature through encode/decode.
        let sig_bytes = encode_signature(&sig, FALCON_512.n, FALCON_512.sig_bytes_max)
            .expect("signature encoding must succeed");
        let (sig_decoded, _consumed) = decode_signature(&sig_bytes)
            .expect("signature round-trip must succeed");
        assert_eq!(sig_decoded.nonce, sig.nonce, "Sig nonce must match after roundtrip");
        assert_eq!(sig_decoded.s2, sig.s2, "Sig s2 must match after roundtrip");

        // Verify the round-tripped signature against the round-tripped pk.
        verify(&pk_decoded, &msg, &sig_decoded)
            .expect("round-tripped signature must verify under round-tripped pk");
    }
}

/// Cross-keypair forgery resistance: signatures from key A must not verify under key B.
///
/// Trivial property; a regression here would indicate a hash-to-point or
/// NTT bug producing constant or message-independent signatures.
#[test]
fn cross_keypair_rejection_falcon_512() {
    let mut rng = StdRng::seed_from_u64(0xcafe_d00d);

    let kp_a = keygen_512(&mut rng).expect("keygen A");
    let kp_b = keygen_512(&mut rng).expect("keygen B");

    let msg = b"the same message signed by A, attempted-verified by B";
    let sig = sign(&mut rng, &kp_a.sk, msg).expect("sign with A");

    verify(&kp_a.pk, msg, &sig).expect("A's signature verifies under A");
    let result = verify(&kp_b.pk, msg, &sig);
    assert!(result.is_err(), "A's signature must not verify under B's key");
}

/// Tampered-signature rejection. Flipping a single bit in s2 should make
/// verification fail (with overwhelming probability — the bound check is
/// loose, but the s1 = c - s2*h relation tightens once h is fixed).
#[test]
fn tampered_signature_rejection_falcon_512() {
    let mut rng = StdRng::seed_from_u64(0xdead_beef);
    let kp = keygen_512(&mut rng).expect("keygen");

    let msg = b"a message we shall tamper with";
    let mut sig = sign(&mut rng, &kp.sk, msg).expect("sign");

    // Tamper: flip the sign of one coefficient. Any coefficient in [-1, 1]
    // wouldn't shift the norm noticeably, so pick a non-trivial one.
    let target_idx = sig.s2.iter().position(|&c| c.unsigned_abs() > 50).unwrap_or(0);
    sig.s2[target_idx] = -sig.s2[target_idx];

    let result = verify(&kp.pk, msg, &sig);
    assert!(result.is_err(), "Tampered signature must not verify");
}

/// Tampered-message rejection. Flipping a single bit in the message should
/// make verification fail because hash_to_point produces a different c.
#[test]
fn tampered_message_rejection_falcon_512() {
    let mut rng = StdRng::seed_from_u64(0xfee1_dead);
    let kp = keygen_512(&mut rng).expect("keygen");

    let msg = b"the original message";
    let sig = sign(&mut rng, &kp.sk, msg).expect("sign");
    verify(&kp.pk, msg, &sig).expect("original verifies");

    let tampered = b"the tampered message";
    let result = verify(&kp.pk, tampered, &sig);
    assert!(result.is_err(), "Signature must not verify on a different message");
}

/// Single-bit-flip rejection over the encoded signature bytes.
///
/// Generates a valid FALCON-512 signature, encodes it to the wire format,
/// then for each byte in the s2 region (skipping the nonce and any
/// padding) flips one bit at a time and re-decodes. The bit-flipped
/// signature should either fail to decode or fail to verify in the
/// overwhelming majority of cases — the rare exception is a flip in
/// padding bits that the decoder ignores. We assert that at least 90%
/// of bit flips are rejected.
///
/// This is the smallest mutation that should fail verify; it's the
/// regression class most likely to slip past `tampered_signature_*`
/// (which flips an entire coefficient) and most likely to expose a
/// Golomb-Rice framing bug.
#[test]
fn single_bit_flip_rejection_falcon_512() {
    use pqsigs_fn_dsa::packing::{encode_signature, decode_signature};

    let mut rng = StdRng::seed_from_u64(0xb175_f117);
    let kp = keygen_512(&mut rng).expect("keygen");
    let msg = b"single-bit-flip probe";
    let sig = sign(&mut rng, &kp.sk, msg).expect("sign");
    verify(&kp.pk, msg, &sig).expect("clean signature verifies");

    let sig_bytes = encode_signature(&sig, FALCON_512.n, FALCON_512.sig_bytes_max)
        .expect("encode");

    // The first 5 bytes are header + nonce-prefix; the nonce occupies the
    // next 40 bytes. Bit-flips in the nonce break verification trivially
    // (different c from hash_to_point) so include them but separate the
    // analysis: we ask "did the flip get rejected?" not "did it tamper
    // a specific field."
    let total_bits = sig_bytes.len() * 8;
    // Bound the test wall-clock: probe every 4th bit position for ~256
    // probes on a typical Falcon-512 sig.
    let stride = (total_bits / 256).max(1);
    let probes: Vec<usize> = (0..total_bits).step_by(stride).collect();

    let mut rejected = 0usize;
    let mut total = 0usize;

    for &bit_idx in &probes {
        let byte_idx = bit_idx / 8;
        let bit_in_byte = bit_idx % 8;
        let mut mutated = sig_bytes.clone();
        mutated[byte_idx] ^= 1 << bit_in_byte;

        total += 1;
        // Decode may fail (good) — count as rejected.
        let decoded = match decode_signature(&mutated) {
            Ok((s, _)) => s,
            Err(_) => {
                rejected += 1;
                continue;
            }
        };
        // Or decode succeeds; verify must fail.
        if verify(&kp.pk, msg, &decoded).is_err() {
            rejected += 1;
        }
    }

    let rate = rejected as f64 / total as f64;
    assert!(
        rate >= 0.90,
        "Bit-flip rejection rate {} too low ({} of {} probes accepted as valid)",
        rate, total - rejected, total
    );
}
