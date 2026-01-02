//! Integration tests for the UOV signature scheme.

use pqsigs_uov::{
    error::UovError,
    field::F,
    keygen::keygen,
    params::{Params, PARAMS_DEMO, PARAMS_NIST_L1},
    sign::sign,
    verify::{verify, verify_bool},
};
use rand::rngs::OsRng;

// Backward compatibility alias
use pqsigs_uov::params::PARAMS_DEMO as PARAMS_L1;

fn flip_first_byte(v: &mut [u8]) {
    if let Some(b) = v.first_mut() {
        *b ^= 0x01;
    }
}

// ============================================================================
// Key Generation Tests
// ============================================================================

#[test]
fn keygen_shapes_are_consistent() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_L1);
    let n = pk.params.n();
    let m = pk.params.m;

    // Public key has m quadratic forms in n variables
    assert_eq!(pk.polys.len(), m);
    for q in &pk.polys {
        assert_eq!(q.coeffs.len(), n * (n + 1) / 2);
    }

    // Secret key matches params and T is n×n
    assert_eq!(sk.params.v + sk.params.m, n);
    assert_eq!(sk.t.len(), n);
    for row in &sk.t {
        assert_eq!(row.len(), n);
    }
    assert_eq!(sk.f_quads.len(), m);
    for ut in &sk.f_quads {
        assert_eq!(ut.len(), n * (n + 1) / 2);
    }
}

#[test]
fn keygen_produces_different_keys() {
    let (pk1, _) = keygen(&mut OsRng, PARAMS_DEMO);
    let (pk2, _) = keygen(&mut OsRng, PARAMS_DEMO);

    // Keys should be different (with overwhelming probability)
    assert_ne!(pk1.polys[0].coeffs, pk2.polys[0].coeffs);
}

// ============================================================================
// Sign/Verify Roundtrip Tests
// ============================================================================

#[test]
fn sign_then_verify_roundtrip() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_L1);
    let msg = b"uov level-1 demo message";

    let sig = sign(&mut OsRng, &pk, &sk, msg).expect("sign should succeed");
    assert!(verify(&pk, msg, &sig).is_ok());
}

#[test]
fn verify_rejects_when_message_changes() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_L1);

    let msg_ok = b"original message";
    let mut msg_bad = msg_ok.to_vec();
    msg_bad[0] ^= 0xFF;

    let sig = sign(&mut OsRng, &pk, &sk, msg_ok).expect("sign should succeed");

    assert!(verify(&pk, msg_ok, &sig).is_ok());
    assert!(
        verify(&pk, &msg_bad, &sig).is_err(),
        "modified message must fail"
    );
}

#[test]
fn verify_rejects_when_salt_changes() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_L1);
    let msg = b"msg";

    let mut sig = sign(&mut OsRng, &pk, &sk, msg).expect("sign should succeed");
    assert!(verify(&pk, msg, &sig).is_ok());

    // Flip one bit in salt
    flip_first_byte(&mut sig.salt);
    assert!(verify(&pk, msg, &sig).is_err(), "tampered salt must fail");
}

#[test]
fn verify_rejects_when_signature_point_changes() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_L1);
    let msg = b"msg";

    let mut sig = sign(&mut OsRng, &pk, &sk, msg).expect("sign should succeed");
    assert!(verify(&pk, msg, &sig).is_ok());

    // Flip one bit in first coordinate of x
    if let Some(first) = sig.x.first_mut() {
        first.0 ^= 0x01;
    }
    assert!(
        verify(&pk, msg, &sig).is_err(),
        "tampered signature point must fail"
    );
}

#[test]
fn multiple_independent_signatures_verify() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_L1);

    for i in 0..10 {
        let msg = format!("test-{}", i);
        let sig = sign(&mut OsRng, &pk, &sk, msg.as_bytes()).expect("sign should succeed");
        assert!(verify(&pk, msg.as_bytes(), &sig).is_ok());
    }
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn verify_returns_invalid_signature_error() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let sig = sign(&mut OsRng, &pk, &sk, b"msg1").unwrap();

    let result = verify(&pk, b"different message", &sig);
    assert!(matches!(result, Err(UovError::InvalidSignature)));
}

#[test]
fn verify_returns_invalid_input_for_wrong_length() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let mut sig = sign(&mut OsRng, &pk, &sk, b"msg").unwrap();

    // Truncate signature
    sig.x.pop();

    let result = verify(&pk, b"msg", &sig);
    assert!(matches!(
        result,
        Err(UovError::InvalidInput {
            field: "signature",
            ..
        })
    ));
}

#[test]
fn verify_bool_works() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"test";
    let sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

    assert!(verify_bool(&pk, msg, &sig));
    assert!(!verify_bool(&pk, b"wrong", &sig));
}

// ============================================================================
// Field Arithmetic Tests
// ============================================================================

#[test]
fn field_add_is_xor() {
    let a = F(0x57);
    let b = F(0x83);
    assert_eq!(a + b, F(0x57 ^ 0x83));
}

#[test]
fn field_sub_equals_add() {
    let a = F(0x57);
    let b = F(0x83);
    assert_eq!(a - b, a + b);
}

#[test]
fn field_mul_identity() {
    let a = F(0x57);
    assert_eq!(a * F::ONE, a);
    assert_eq!(a * F::ZERO, F::ZERO);
}

#[test]
fn field_inverse_roundtrip() {
    for val in 1u8..=255 {
        let a = F(val);
        let a_inv = a.inverse();
        assert_eq!(a * a_inv, F::ONE, "inverse failed for {}", val);
    }
}

// ============================================================================
// Parameter Tests
// ============================================================================

#[test]
fn params_validation() {
    assert!(Params::new(16, 8).is_ok());
    assert!(Params::new(0, 8).is_err());
    assert!(Params::new(16, 0).is_err());
    assert!(Params::new(4, 8).is_err()); // v < m is invalid
}

#[test]
fn params_demo_matches_expected() {
    assert_eq!(PARAMS_DEMO.v, 16);
    assert_eq!(PARAMS_DEMO.m, 8);
    assert_eq!(PARAMS_DEMO.n(), 24);
}

#[test]
fn params_nist_l1_matches_expected() {
    assert_eq!(PARAMS_NIST_L1.v, 68);
    assert_eq!(PARAMS_NIST_L1.m, 44);
    assert_eq!(PARAMS_NIST_L1.n(), 112);
}

// ============================================================================
// Edge Case Tests
// ============================================================================

#[test]
fn sign_empty_message() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"";

    let sig = sign(&mut OsRng, &pk, &sk, msg).expect("signing empty message should work");
    assert!(verify(&pk, msg, &sig).is_ok());
}

#[test]
fn sign_long_message() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = vec![0xAB; 10000]; // 10KB message

    let sig = sign(&mut OsRng, &pk, &sk, &msg).expect("signing long message should work");
    assert!(verify(&pk, &msg, &sig).is_ok());
}

#[test]
fn same_message_different_signatures() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"same message";

    let sig1 = sign(&mut OsRng, &pk, &sk, msg).unwrap();
    let sig2 = sign(&mut OsRng, &pk, &sk, msg).unwrap();

    // Both should verify
    assert!(verify(&pk, msg, &sig1).is_ok());
    assert!(verify(&pk, msg, &sig2).is_ok());

    // But should be different (randomized signatures)
    assert!(sig1.salt != sig2.salt || sig1.x != sig2.x);
}
