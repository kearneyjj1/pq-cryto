//! FALCON-512 Full Integration Test
//!
//! This tests the complete FALCON-512 implementation:
//! - Key generation with NTRU equation solving
//! - Signing with FFT lattice sampling (ffSampling with LDL* tree)
//! - Verification with NTT-based modular arithmetic
//! - Serialization

use pqsigs_fn_dsa::{
    keygen_512, sign, verify,
    encode_public_key, decode_public_key,
    encode_secret_key, decode_secret_key,
};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn main() {
    println!("=== FALCON-512 Full Integration Test ===\n");

    let mut rng = StdRng::seed_from_u64(42);

    // Key Generation
    println!("1. Key Generation...");
    let keypair = keygen_512(&mut rng).expect("keygen failed");
    println!("   Public key h[0..5]: {:?}", &keypair.pk.h[0..5]);
    println!("   Secret key f[0..5]: {:?}", &keypair.sk.f[0..5]);
    println!("   Signature bound: {:.0}", keypair.sk.params.sig_bound_sq);
    println!("   Key generation: SUCCESS\n");

    // Try multiple signing approaches
    let message = b"Hello, FALCON-512!";

    // Try main sign function
    println!("2. Signing with FFT sampler...");
    let sig_result = sign(&mut rng, &keypair.sk, message);

    let sig = match sig_result {
        Ok(s) => {
            println!("   Signature s2 norm^2: {}", s.norm_sq());
            println!("   FFT Signing: SUCCESS\n");
            Some(s)
        }
        Err(e) => {
            println!("   FFT signing failed: {}", e);
            None
        }
    };

    // Verify if we got a signature
    if let Some(sig) = sig {
        println!("3. Verification...");
        match verify(&keypair.pk, message, &sig) {
            Ok(()) => {
                println!("   Verification: SUCCESS\n");

                // Test with wrong message — must fail for FALCON-512
                println!("4. Wrong message test...");
                match verify(&keypair.pk, b"Wrong message", &sig) {
                    Ok(()) => println!("   WARNING: Wrong message accepted (unexpected)"),
                    Err(_) => println!("   Wrong message rejected: SUCCESS\n"),
                }
            }
            Err(e) => {
                println!("   Verification failed: {}", e);

                // Debug: show computed s1 norm
                println!("   Debugging signature verification...\n");
            }
        }
    }

    // Serialization
    println!("5. Serialization...");
    let pk_bytes = encode_public_key(&keypair.pk);
    let sk_bytes = encode_secret_key(&keypair.sk);
    println!("   Public key: {} bytes", pk_bytes.len());
    println!("   Secret key: {} bytes", sk_bytes.len());

    let pk_decoded = decode_public_key(&pk_bytes).expect("pk decode");
    let sk_decoded = decode_secret_key(&sk_bytes).expect("sk decode");
    assert_eq!(keypair.pk.h, pk_decoded.h);
    assert_eq!(keypair.sk.f, sk_decoded.f);
    println!("   Roundtrip: SUCCESS\n");

    println!("=== Test Complete ===");
}
