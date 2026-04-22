//! End-to-end integration tests for the SQI-SIGN implementation.
//!
//! These tests exercise the full public API (`keygen` → `sign` → `verify`)
//! plus serialization round-trips and basic forgery resistance checks.
//!
//! Uses the *test* parameter sets (`SQISIGN_NIST_I/III/V`) which are built on
//! small primes (~54/74/84 bits) so the suite runs quickly. Production
//! parameter sets are exercised in a separate, `#[ignore]`d slow test below.

use pqsigs_sqi_sign::{
    keygen, sign, verify,
    params::{SQISIGN_NIST_I, SQISIGN_NIST_III, SQISIGN_NIST_V},
    PublicKey, SecretKey, Signature,
};
use rand::{rngs::StdRng, SeedableRng};
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn seeded_rng(tag: u64) -> StdRng {
    // Deterministic RNG so test failures are reproducible.
    StdRng::seed_from_u64(tag)
}

/// Runs `f` on a worker thread and fails the test if it doesn't return within
/// `timeout`. Without this, `sign()` currently hangs on the test parameter
/// sets (see audit 2026-04-22) and the test harness would never exit.
///
/// The worker thread is leaked on timeout — acceptable for a test binary
/// that is about to exit anyway.
fn with_timeout<T, F>(timeout: Duration, label: &str, f: F) -> T
where
    F: FnOnce() -> T + Send + 'static,
    T: Send + 'static,
{
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let _ = tx.send(f());
    });
    match rx.recv_timeout(timeout) {
        Ok(v) => v,
        Err(_) => panic!(
            "`{}` did not complete within {:?} — probable non-termination \
             in sign()/verify() (see audit 2026-04-22)",
            label, timeout,
        ),
    }
}

const SIGN_TIMEOUT: Duration = Duration::from_secs(30);

/// Calls `sign` on a worker thread with a wall-clock timeout. Panics if
/// signing does not return in time.
fn timed_sign(
    rng_seed: u64,
    sk: SecretKey,
    msg: Vec<u8>,
) -> pqsigs_sqi_sign::Result<Signature> {
    with_timeout(SIGN_TIMEOUT, "sign", move || {
        let mut rng = StdRng::seed_from_u64(rng_seed);
        sign(&mut rng, &sk, &msg)
    })
}

/// Calls `verify` on a worker thread with a wall-clock timeout.
fn timed_verify(
    pk: PublicKey,
    msg: Vec<u8>,
    sig: Signature,
) -> pqsigs_sqi_sign::Result<()> {
    with_timeout(SIGN_TIMEOUT, "verify", move || {
        verify(&pk, &msg, &sig)
    })
}

// ---------------------------------------------------------------------------
// Core roundtrip: keygen -> sign -> verify
// ---------------------------------------------------------------------------

#[test]
fn roundtrip_level_i() {
    let mut rng = seeded_rng(0xA11CE);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let msg = b"hello SQI-SIGN level I".to_vec();
    let sig = timed_sign(0xA11CE, sk, msg.clone()).expect("signing must not fail");
    timed_verify(pk, msg, sig).expect("honest signature must verify");
}

#[test]
fn roundtrip_level_iii() {
    let mut rng = seeded_rng(0xB0B);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_III.clone());

    let msg = b"hello SQI-SIGN level III".to_vec();
    let sig = timed_sign(0xB0B, sk, msg.clone()).expect("signing must not fail");
    timed_verify(pk, msg, sig).expect("honest signature must verify");
}

#[test]
fn roundtrip_level_v() {
    let mut rng = seeded_rng(0xCA7);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_V.clone());

    let msg = b"hello SQI-SIGN level V".to_vec();
    let sig = timed_sign(0xCA7, sk, msg.clone()).expect("signing must not fail");
    timed_verify(pk, msg, sig).expect("honest signature must verify");
}

#[test]
fn roundtrip_empty_message() {
    let mut rng = seeded_rng(0xDEADBEEF);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let sig = timed_sign(0xDEADBEEF, sk, Vec::new()).expect("empty-message signing must not fail");
    timed_verify(pk, Vec::new(), sig).expect("empty-message signature must verify");
}

#[test]
fn roundtrip_long_message() {
    let mut rng = seeded_rng(0xFEEDFACE);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let msg = vec![0xAB_u8; 64 * 1024]; // 64 KiB
    let sig = timed_sign(0xFEEDFACE, sk, msg.clone()).expect("long-message signing must not fail");
    timed_verify(pk, msg, sig).expect("long-message signature must verify");
}

// ---------------------------------------------------------------------------
// Forgery resistance
// ---------------------------------------------------------------------------

#[test]
fn reject_tampered_message() {
    let mut rng = seeded_rng(0x1234);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let msg = b"authentic message".to_vec();
    let sig = timed_sign(0x1234, sk, msg.clone()).expect("signing must not fail");

    // Sanity: honest sig verifies first.
    timed_verify(pk.clone(), msg, sig.clone()).expect("honest signature must verify");

    let tampered = b"authentic message!".to_vec();
    assert!(
        timed_verify(pk, tampered, sig).is_err(),
        "verifier must reject signature over a different message",
    );
}

#[test]
fn reject_cross_key_verification() {
    let mut rng = seeded_rng(0x5678);
    let (pk_a, sk_a) = keygen(&mut rng, SQISIGN_NIST_I.clone());
    let (pk_b, _sk_b) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let msg = b"signed by A, verified against B's pk".to_vec();
    let sig = timed_sign(0x5678, sk_a, msg.clone()).expect("signing must not fail");

    timed_verify(pk_a, msg.clone(), sig.clone())
        .expect("honest signature must verify under signer's pk");
    assert!(
        timed_verify(pk_b, msg, sig).is_err(),
        "verifier must reject A's signature under B's public key",
    );
}

#[test]
fn reject_corrupted_signature_bytes() {
    let mut rng = seeded_rng(0x9ABC);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let msg = b"message whose signature we will mangle".to_vec();
    let sig = timed_sign(0x9ABC, sk, msg.clone()).expect("signing must not fail");

    let mut raw = sig.to_bytes();
    // Flip a byte in the middle (targets the commitment curve region).
    let mid = raw.len() / 2;
    raw[mid] ^= 0xFF;

    let modulus = pk.curve.modulus.clone();
    match Signature::from_bytes(&raw, &modulus) {
        // If deserialization rejects the malformed bytes, that's a valid outcome.
        Err(_) => {}
        // If it accepts them, verification must reject.
        Ok(bad_sig) => assert!(
            timed_verify(pk, msg, bad_sig).is_err(),
            "verifier must reject a tampered signature",
        ),
    }
}

// ---------------------------------------------------------------------------
// Serialization round-trips (orthogonal to the crypto core)
// ---------------------------------------------------------------------------

#[test]
fn public_key_serialization_roundtrip() {
    let mut rng = seeded_rng(0xDEF0);
    let (pk, _sk) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let bytes = pk.to_bytes();
    let pk2 = PublicKey::from_bytes(&bytes, &SQISIGN_NIST_I).expect("pk deserialization");

    assert_eq!(
        pk.to_bytes(),
        pk2.to_bytes(),
        "public key must round-trip through to_bytes/from_bytes",
    );
}

#[test]
fn signature_serialization_roundtrip() {
    let mut rng = seeded_rng(0x1357);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let msg = b"serialize me".to_vec();
    let sig = timed_sign(0x1357, sk, msg.clone()).expect("signing must not fail");

    let bytes = sig.to_bytes();
    let modulus = pk.curve.modulus.clone();
    let sig2 = Signature::from_bytes(&bytes, &modulus).expect("sig deserialization");

    assert_eq!(
        sig.to_bytes(),
        sig2.to_bytes(),
        "signature must round-trip through to_bytes/from_bytes",
    );

    // A deserialized signature over the same message should still verify
    // (if the honest signature did).
    if timed_verify(pk.clone(), msg.clone(), sig).is_ok() {
        timed_verify(pk, msg, sig2).expect("round-tripped signature must verify");
    }
}

/// Regression guard for CRIT-04 in the 2026-04-22 audit:
/// `SecretKey::from_bytes` currently ignores its input and returns a key
/// derived from a fixed all-zero seed, so any `to_bytes → from_bytes` round
/// trip will produce an unrelated secret key.
///
/// This test is expected to fail until `from_bytes` is implemented properly.
/// It is intentionally *not* `#[ignore]`'d so CI reflects the real state.
#[test]
fn secret_key_serialization_roundtrip() {
    let mut rng = seeded_rng(0x2468);
    let (_pk, sk) = keygen(&mut rng, SQISIGN_NIST_I.clone());

    let bytes = sk.to_bytes();
    let sk2 = SecretKey::from_bytes(&bytes, &SQISIGN_NIST_I).expect("sk deserialization");

    assert_eq!(
        sk.to_bytes(),
        sk2.to_bytes(),
        "secret key must round-trip through to_bytes/from_bytes \
         (expected failure: see audit finding CRIT-04)",
    );
}

// ---------------------------------------------------------------------------
// Determinism / keygen_internal
// ---------------------------------------------------------------------------

#[test]
fn keygen_internal_is_deterministic() {
    use pqsigs_sqi_sign::keygen::keygen_internal;

    let seed = [0x42u8; 32];
    let (pk1, sk1) = keygen_internal(&seed, SQISIGN_NIST_I.clone());
    let (pk2, sk2) = keygen_internal(&seed, SQISIGN_NIST_I.clone());

    assert_eq!(
        pk1.to_bytes(),
        pk2.to_bytes(),
        "keygen_internal must produce the same public key for the same seed",
    );
    assert_eq!(
        sk1.to_bytes(),
        sk2.to_bytes(),
        "keygen_internal must produce the same secret key for the same seed",
    );
}

#[test]
fn keygen_is_non_deterministic_across_calls() {
    // Two different RNG seeds must produce different keys (with overwhelming probability).
    let mut rng_a = seeded_rng(1);
    let mut rng_b = seeded_rng(2);

    let (pk_a, _) = keygen(&mut rng_a, SQISIGN_NIST_I.clone());
    let (pk_b, _) = keygen(&mut rng_b, SQISIGN_NIST_I.clone());

    assert_ne!(
        pk_a.to_bytes(),
        pk_b.to_bytes(),
        "distinct RNG states must yield distinct public keys",
    );
}

// ---------------------------------------------------------------------------
// Slow test against production parameter sets.
// Marked #[ignore] so `cargo test` stays fast; run with `--ignored`.
// ---------------------------------------------------------------------------

#[test]
#[ignore = "uses 256-bit production primes; opt-in via --ignored"]
fn roundtrip_production_level_i() {
    use pqsigs_sqi_sign::params::SQISIGN_NIST_I_PROD;

    let mut rng = seeded_rng(0xC0FFEE);
    let (pk, sk) = keygen(&mut rng, SQISIGN_NIST_I_PROD.clone());

    let msg = b"production-prime smoke test".to_vec();
    let sig = timed_sign(0xC0FFEE, sk, msg.clone()).expect("signing must not fail");
    timed_verify(pk, msg, sig).expect("honest signature must verify");
}
