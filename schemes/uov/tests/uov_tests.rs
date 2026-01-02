//! Integration tests for the UOV signature scheme.

use pqsigs_uov::{
    error::UovError,
    field::F,
    keygen::keygen,
    params::{Params, PARAMS_DEMO, PARAMS_NIST_L1, PARAMS_NIST_L3, PARAMS_NIST_L5},
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

// ============================================================================
// NIST Parameter Set Tests
// ============================================================================

#[test]
fn roundtrip_nist_l1() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_NIST_L1);
    let msg = b"NIST Level 1 test message";

    let sig = sign(&mut OsRng, &pk, &sk, msg).expect("signing should succeed");
    assert!(verify(&pk, msg, &sig).is_ok());
}

#[test]
#[ignore = "L3 keygen is slow in debug mode, run with --release"]
fn roundtrip_nist_l3() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_NIST_L3);
    let msg = b"NIST Level 3 test message";

    let sig = sign(&mut OsRng, &pk, &sk, msg).expect("signing should succeed");
    assert!(verify(&pk, msg, &sig).is_ok());
}

#[test]
#[ignore = "L5 keygen is slow (~30s), run with --ignored"]
fn roundtrip_nist_l5() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_NIST_L5);
    let msg = b"NIST Level 5 test message";

    let sig = sign(&mut OsRng, &pk, &sk, msg).expect("signing should succeed");
    assert!(verify(&pk, msg, &sig).is_ok());
}

#[test]
fn nist_l1_multiple_messages() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_NIST_L1);

    for i in 0..5 {
        let msg = format!("NIST L1 message {}", i);
        let sig = sign(&mut OsRng, &pk, &sk, msg.as_bytes()).expect("signing should succeed");
        assert!(verify(&pk, msg.as_bytes(), &sig).is_ok());
    }
}

#[test]
#[ignore = "L3 keygen is slow in debug mode, run with --release"]
fn nist_l3_multiple_messages() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_NIST_L3);

    for i in 0..3 {
        let msg = format!("NIST L3 message {}", i);
        let sig = sign(&mut OsRng, &pk, &sk, msg.as_bytes()).expect("signing should succeed");
        assert!(verify(&pk, msg.as_bytes(), &sig).is_ok());
    }
}

// ============================================================================
// Cross-Key and Cross-Parameter Tests
// ============================================================================

#[test]
fn verify_rejects_wrong_key() {
    let (pk1, sk1) = keygen(&mut OsRng, PARAMS_DEMO);
    let (pk2, _sk2) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"test message";

    let sig = sign(&mut OsRng, &pk1, &sk1, msg).expect("signing should succeed");

    // Should verify with correct key
    assert!(verify(&pk1, msg, &sig).is_ok());

    // Should fail with wrong key
    assert!(verify(&pk2, msg, &sig).is_err());
}

#[test]
fn verify_rejects_wrong_key_nist_l1() {
    let (pk1, sk1) = keygen(&mut OsRng, PARAMS_NIST_L1);
    let (pk2, _sk2) = keygen(&mut OsRng, PARAMS_NIST_L1);
    let msg = b"test message";

    let sig = sign(&mut OsRng, &pk1, &sk1, msg).expect("signing should succeed");

    assert!(verify(&pk1, msg, &sig).is_ok());
    assert!(verify(&pk2, msg, &sig).is_err());
}

#[test]
fn cross_parameter_rejection_demo_vs_l1() {
    let (pk_demo, sk_demo) = keygen(&mut OsRng, PARAMS_DEMO);
    let (pk_l1, _sk_l1) = keygen(&mut OsRng, PARAMS_NIST_L1);
    let msg = b"cross-param test";

    let sig = sign(&mut OsRng, &pk_demo, &sk_demo, msg).expect("signing should succeed");

    // Should verify with matching params
    assert!(verify(&pk_demo, msg, &sig).is_ok());

    // Should fail with different params (different signature length)
    assert!(verify(&pk_l1, msg, &sig).is_err());
}

#[test]
#[ignore = "L3 keygen is slow in debug mode, run with --release"]
fn cross_parameter_rejection_l1_vs_l3() {
    let (pk_l1, sk_l1) = keygen(&mut OsRng, PARAMS_NIST_L1);
    let (pk_l3, _sk_l3) = keygen(&mut OsRng, PARAMS_NIST_L3);
    let msg = b"L1 vs L3 test";

    let sig = sign(&mut OsRng, &pk_l1, &sk_l1, msg).expect("signing should succeed");

    assert!(verify(&pk_l1, msg, &sig).is_ok());
    assert!(verify(&pk_l3, msg, &sig).is_err());
}

// ============================================================================
// Key and Signature Size Tests
// ============================================================================

#[test]
fn key_sizes_demo() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);

    // Public key: m quadratic forms, each with ut_size coefficients
    assert_eq!(pk.polys.len(), PARAMS_DEMO.m);
    for poly in &pk.polys {
        assert_eq!(poly.coeffs.len(), PARAMS_DEMO.ut_size());
    }

    // Secret key: n×n matrix T
    assert_eq!(sk.t.len(), PARAMS_DEMO.n());
    for row in &sk.t {
        assert_eq!(row.len(), PARAMS_DEMO.n());
    }
}

#[test]
fn key_sizes_nist_l1() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_NIST_L1);

    assert_eq!(pk.polys.len(), PARAMS_NIST_L1.m);
    for poly in &pk.polys {
        assert_eq!(poly.coeffs.len(), PARAMS_NIST_L1.ut_size());
    }

    assert_eq!(sk.t.len(), PARAMS_NIST_L1.n());
    for row in &sk.t {
        assert_eq!(row.len(), PARAMS_NIST_L1.n());
    }
}

#[test]
fn signature_sizes() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let sig = sign(&mut OsRng, &pk, &sk, b"test").unwrap();

    assert_eq!(sig.x.len(), PARAMS_DEMO.signature_size());
    assert_eq!(sig.salt.len(), PARAMS_DEMO.salt_size());
}

#[test]
fn signature_sizes_nist_l1() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_NIST_L1);
    let sig = sign(&mut OsRng, &pk, &sk, b"test").unwrap();

    assert_eq!(sig.x.len(), PARAMS_NIST_L1.signature_size());
    assert_eq!(sig.salt.len(), PARAMS_NIST_L1.salt_size());
}

// ============================================================================
// Comprehensive Field Arithmetic Tests
// ============================================================================

#[test]
fn field_add_commutative() {
    for a in 0..=255u8 {
        for b in 0..=255u8 {
            assert_eq!(F(a) + F(b), F(b) + F(a));
        }
    }
}

#[test]
fn field_add_associative() {
    let vals = [0u8, 1, 127, 128, 255];
    for &a in &vals {
        for &b in &vals {
            for &c in &vals {
                assert_eq!((F(a) + F(b)) + F(c), F(a) + (F(b) + F(c)));
            }
        }
    }
}

#[test]
fn field_mul_commutative() {
    for a in 0..=255u8 {
        for b in 0..=255u8 {
            assert_eq!(F(a) * F(b), F(b) * F(a));
        }
    }
}

#[test]
fn field_mul_associative() {
    let vals = [0u8, 1, 2, 127, 128, 255];
    for &a in &vals {
        for &b in &vals {
            for &c in &vals {
                assert_eq!((F(a) * F(b)) * F(c), F(a) * (F(b) * F(c)));
            }
        }
    }
}

#[test]
fn field_distributive() {
    let vals = [0u8, 1, 2, 127, 128, 255];
    for &a in &vals {
        for &b in &vals {
            for &c in &vals {
                // a * (b + c) = a*b + a*c
                assert_eq!(F(a) * (F(b) + F(c)), F(a) * F(b) + F(a) * F(c));
            }
        }
    }
}

#[test]
fn field_zero_absorbs() {
    for a in 0..=255u8 {
        assert_eq!(F(a) * F::ZERO, F::ZERO);
        assert_eq!(F::ZERO * F(a), F::ZERO);
    }
}

#[test]
fn field_one_identity() {
    for a in 0..=255u8 {
        assert_eq!(F(a) * F::ONE, F(a));
        assert_eq!(F::ONE * F(a), F(a));
    }
}

#[test]
fn field_add_inverse() {
    // In GF(2^8), every element is its own additive inverse
    for a in 0..=255u8 {
        assert_eq!(F(a) + F(a), F::ZERO);
    }
}

#[test]
fn field_inverse_zero() {
    // Inverse of zero is defined as zero (for safety)
    assert_eq!(F::ZERO.inverse(), F::ZERO);
}

#[test]
fn field_square_all() {
    // Verify squaring works for all elements
    for a in 0..=255u8 {
        assert_eq!(F(a).square(), F(a) * F(a));
    }
}

// ============================================================================
// Matrix Operation Tests
// ============================================================================

#[test]
fn matrix_identity_multiplication() {
    use pqsigs_uov::field::F;

    // Create 3x3 identity matrix
    let n = 3;
    let identity: Vec<Vec<F>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { F::ONE } else { F::ZERO }).collect())
        .collect();

    // Create a random-ish matrix
    let a: Vec<Vec<F>> = vec![
        vec![F(1), F(2), F(3)],
        vec![F(4), F(5), F(6)],
        vec![F(7), F(8), F(9)],
    ];

    // A * I = A
    let result = mat_mul_test(&a, &identity);
    assert_eq!(result, a);

    // I * A = A
    let result2 = mat_mul_test(&identity, &a);
    assert_eq!(result2, a);
}

// Helper function for matrix multiplication test
fn mat_mul_test(a: &[Vec<F>], b: &[Vec<F>]) -> Vec<Vec<F>> {
    let n = a.len();
    let mut c = vec![vec![F::ZERO; n]; n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i][k];
            if !aik.is_zero() {
                for j in 0..n {
                    c[i][j] += aik * b[k][j];
                }
            }
        }
    }
    c
}

#[test]
fn matrix_transpose_involutive() {
    // Transpose of transpose is original
    let a: Vec<Vec<F>> = vec![
        vec![F(1), F(2), F(3)],
        vec![F(4), F(5), F(6)],
        vec![F(7), F(8), F(9)],
    ];

    let t1 = transpose_test(&a);
    let t2 = transpose_test(&t1);
    assert_eq!(t2, a);
}

fn transpose_test(a: &[Vec<F>]) -> Vec<Vec<F>> {
    let n = a.len();
    let mut t = vec![vec![F::ZERO; n]; n];
    for i in 0..n {
        for j in 0..n {
            t[j][i] = a[i][j];
        }
    }
    t
}

// ============================================================================
// Additional Tampering Tests
// ============================================================================

#[test]
fn tampering_all_signature_bytes() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"tampering test";
    let sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

    // Tamper with each byte of signature point
    for i in 0..sig.x.len() {
        let mut tampered = sig.clone();
        tampered.x[i].0 ^= 0xFF;
        assert!(verify(&pk, msg, &tampered).is_err(), "tampering x[{}] should fail", i);
    }
}

#[test]
fn tampering_all_salt_bytes() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"salt tampering test";
    let sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

    // Tamper with each byte of salt
    for i in 0..sig.salt.len() {
        let mut tampered = sig.clone();
        tampered.salt[i] ^= 0xFF;
        assert!(verify(&pk, msg, &tampered).is_err(), "tampering salt[{}] should fail", i);
    }
}

#[test]
fn verify_consistency_multiple_runs() {
    // Verify that verification is deterministic
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"consistency test";
    let sig = sign(&mut OsRng, &pk, &sk, msg).unwrap();

    for _ in 0..100 {
        assert!(verify(&pk, msg, &sig).is_ok());
    }
}

// ============================================================================
// Stress Tests
// ============================================================================

#[test]
fn stress_many_signatures_demo() {
    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);

    for i in 0..50 {
        let msg = format!("stress test message {}", i);
        let sig = sign(&mut OsRng, &pk, &sk, msg.as_bytes()).expect("signing should succeed");
        assert!(verify(&pk, msg.as_bytes(), &sig).is_ok());
    }
}

#[test]
fn stress_many_keypairs_demo() {
    for _ in 0..10 {
        let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
        let msg = b"keypair stress test";
        let sig = sign(&mut OsRng, &pk, &sk, msg).expect("signing should succeed");
        assert!(verify(&pk, msg, &sig).is_ok());
    }
}
