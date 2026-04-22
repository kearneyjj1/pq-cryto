//! FIPS 206 / Falcon-512 Known Answer Test (KAT) verification tests.
//!
//! These tests validate our FN-DSA implementation against the official
//! NIST PQC Falcon-512 KAT vectors from the reference implementation.
//!
//! The KAT file (`falcon512-KAT.rsp`) contains 100 test vectors, each with:
//! - A deterministic seed (for the reference implementation's DRBG)
//! - A message of variable length
//! - A public key (pk) in NIST API format
//! - A secret key (sk) in NIST API format
//! - A signed message (sm = sig || msg) in NIST API format
//!
//! Since our implementation uses different internal randomness (not the NIST
//! AES-CTR-DRBG), we cannot reproduce keygen/signing byte-for-byte. Instead,
//! we perform **verification-only** testing: we parse the reference pk and
//! signature from each KAT vector and verify them using our verification code.
//! This validates that our:
//! - Public key decoding (14-bit packed coefficients) is correct
//! - Golomb-Rice signature decompression is correct
//! - Hash-to-point (SHAKE256 rejection sampling) matches the reference
//! - NTT-based polynomial multiplication is correct
//! - Norm bound checking is correct

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use pqsigs_fn_dsa::packing::{decode_nist_public_key, from_hex, parse_nist_signed_message};
use pqsigs_fn_dsa::verify::verify;

/// Parse a NIST KAT .rsp file into a vector of test vectors.
///
/// Each test vector is a HashMap with keys: count, seed, mlen, msg, pk, sk, smlen, sm
fn parse_kat_file(path: &Path) -> Vec<HashMap<String, String>> {
    let content = fs::read_to_string(path).expect("Failed to read KAT file");
    let mut vectors = Vec::new();
    let mut current: HashMap<String, String> = HashMap::new();

    for line in content.lines() {
        let line = line.trim();

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            if !current.is_empty() {
                vectors.push(current.clone());
                current.clear();
            }
            continue;
        }

        // Parse key = value
        if let Some(eq_pos) = line.find(" = ") {
            let key = line[..eq_pos].trim().to_string();
            let value = line[eq_pos + 3..].trim().to_string();
            current.insert(key, value);
        }
    }

    // Don't forget the last vector
    if !current.is_empty() {
        vectors.push(current);
    }

    vectors
}

/// Verify a single KAT vector.
///
/// Returns Ok(()) if verification succeeds, Err with details if it fails.
fn verify_kat_vector(vector: &HashMap<String, String>) -> Result<(), String> {
    let count: usize = vector
        .get("count")
        .ok_or("missing count")?
        .parse()
        .map_err(|e| format!("bad count: {}", e))?;

    let mlen: usize = vector
        .get("mlen")
        .ok_or("missing mlen")?
        .parse()
        .map_err(|e| format!("bad mlen: {}", e))?;

    let pk_hex = vector.get("pk").ok_or("missing pk")?;
    let sm_hex = vector.get("sm").ok_or("missing sm")?;
    let msg_hex = vector.get("msg").ok_or("missing msg")?;

    // Decode hex values
    let pk_bytes = from_hex(pk_hex).map_err(|e| format!("KAT {}: bad pk hex: {}", count, e))?;
    let sm_bytes = from_hex(sm_hex).map_err(|e| format!("KAT {}: bad sm hex: {}", count, e))?;
    let msg_bytes = from_hex(msg_hex).map_err(|e| format!("KAT {}: bad msg hex: {}", count, e))?;

    // Decode the NIST-format public key
    let pk = decode_nist_public_key(&pk_bytes)
        .map_err(|e| format!("KAT {}: pk decode failed: {}", count, e))?;

    // Parse the signed message to extract nonce, s2, and message
    let (extracted_msg, sig) = parse_nist_signed_message(&sm_bytes, mlen)
        .map_err(|e| format!("KAT {}: sm parse failed: {}", count, e))?;

    // Verify that the extracted message matches the expected message
    if extracted_msg != msg_bytes {
        return Err(format!(
            "KAT {}: extracted message does not match expected message",
            count
        ));
    }

    // Verify the signature using our implementation
    verify(&pk, &extracted_msg, &sig)
        .map_err(|e| format!("KAT {}: verification FAILED: {}", count, e))?;

    Ok(())
}

#[test]
fn test_kat_falcon512_first_vector() {
    let kat_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/falcon512-KAT.rsp");
    if !kat_path.exists() {
        eprintln!(
            "KAT file not found at {:?}, skipping test",
            kat_path
        );
        return;
    }

    let vectors = parse_kat_file(&kat_path);
    assert!(!vectors.is_empty(), "KAT file contains no test vectors");

    // Test just the first vector as a quick smoke test
    verify_kat_vector(&vectors[0]).expect("First KAT vector verification failed");

    println!("KAT vector 0: PASSED");
}

#[test]
fn test_kat_falcon512_first_10_vectors() {
    let kat_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/falcon512-KAT.rsp");
    if !kat_path.exists() {
        eprintln!(
            "KAT file not found at {:?}, skipping test",
            kat_path
        );
        return;
    }

    let vectors = parse_kat_file(&kat_path);
    let count = vectors.len().min(10);

    let mut passed = 0;
    let mut failed = 0;

    for i in 0..count {
        match verify_kat_vector(&vectors[i]) {
            Ok(()) => {
                passed += 1;
            }
            Err(e) => {
                eprintln!("FAIL: {}", e);
                failed += 1;
            }
        }
    }

    println!(
        "KAT Falcon-512 first 10: {}/{} passed, {} failed",
        passed, count, failed
    );
    assert_eq!(failed, 0, "{} KAT vectors failed verification", failed);
}

#[test]
fn test_kat_falcon512_all_vectors() {
    let kat_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/falcon512-KAT.rsp");
    if !kat_path.exists() {
        eprintln!(
            "KAT file not found at {:?}, skipping test",
            kat_path
        );
        return;
    }

    let vectors = parse_kat_file(&kat_path);
    assert!(
        vectors.len() >= 100,
        "Expected at least 100 KAT vectors, got {}",
        vectors.len()
    );

    let mut passed = 0;
    let mut failed = 0;
    let mut errors: Vec<String> = Vec::new();

    for (i, vector) in vectors.iter().enumerate() {
        match verify_kat_vector(vector) {
            Ok(()) => {
                passed += 1;
            }
            Err(e) => {
                errors.push(e);
                failed += 1;
            }
        }

        // Print progress every 10 vectors
        if (i + 1) % 10 == 0 {
            println!("Progress: {}/{} vectors processed...", i + 1, vectors.len());
        }
    }

    println!(
        "\nKAT Falcon-512 ALL: {}/{} passed, {} failed",
        passed,
        vectors.len(),
        failed
    );

    if !errors.is_empty() {
        for e in &errors {
            eprintln!("  FAIL: {}", e);
        }
        panic!(
            "{} of {} KAT vectors failed verification",
            failed,
            vectors.len()
        );
    }
}

/// Test that the NIST pk decoding produces the expected number of coefficients
/// and that all are in valid range [0, q).
#[test]
fn test_nist_pk_decode_format() {
    let kat_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/falcon512-KAT.rsp");
    if !kat_path.exists() {
        return;
    }

    let vectors = parse_kat_file(&kat_path);
    let pk_hex = vectors[0].get("pk").expect("missing pk");
    let pk_bytes = from_hex(pk_hex).expect("bad hex");

    assert_eq!(pk_bytes.len(), 897, "FALCON-512 pk should be 897 bytes");
    assert_eq!(pk_bytes[0], 0x09, "FALCON-512 pk header should be 0x09");

    let pk = decode_nist_public_key(&pk_bytes).expect("pk decode failed");
    assert_eq!(pk.h.len(), 512, "Should have 512 coefficients");
    assert_eq!(pk.params.n, 512, "Should be FALCON-512 params");

    // All coefficients should be in [0, q)
    for (i, &coeff) in pk.h.iter().enumerate() {
        assert!(
            coeff >= 0 && (coeff as i32) < 12289,
            "Coefficient {} = {} out of range [0, q)",
            i,
            coeff
        );
    }
}

/// Test that the NIST sm parsing correctly extracts the message.
#[test]
fn test_nist_sm_message_extraction() {
    let kat_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/falcon512-KAT.rsp");
    if !kat_path.exists() {
        return;
    }

    let vectors = parse_kat_file(&kat_path);

    for (i, vector) in vectors.iter().enumerate().take(10) {
        let mlen: usize = vector.get("mlen").unwrap().parse().unwrap();
        let msg_hex = vector.get("msg").unwrap();
        let sm_hex = vector.get("sm").unwrap();

        let msg_bytes = from_hex(msg_hex).expect("bad msg hex");
        let sm_bytes = from_hex(sm_hex).expect("bad sm hex");

        let (extracted_msg, sig) =
            parse_nist_signed_message(&sm_bytes, mlen).expect("sm parse failed");

        assert_eq!(
            extracted_msg, msg_bytes,
            "KAT {}: extracted message mismatch",
            i
        );
        assert_eq!(
            sig.s2.len(),
            512,
            "KAT {}: s2 should have 512 coefficients",
            i
        );
    }
}
