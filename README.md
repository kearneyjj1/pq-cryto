# pq-crypto

Rust implementations of post-quantum digital signature schemes. Built for learning and experimentation, featuring key generation, signing, verification, and comprehensive test suites.

**Warning**: These implementations are for educational and experimental purposes only and are NOT cleared for production deployment. See [SECURITY.md](SECURITY.md) for the current security posture.

## Security Status

This codebase has undergone multiple rounds of security review and hardening.

- **Secure Memory Handling**: Secret keys implement `Drop` with zeroization via the `zeroize` crate. FN-DSA additionally zeroizes recursive ffSampling intermediates and FFT-domain copies of `(f, g, F, G)` per signing call.
- **Type-System Secret Containment** (FN-DSA): All secret-material fields are `pub(crate)`; `SecretKey` does not implement `Clone` (use explicit `try_clone()`); `sign_simple` / `sign_with_nonce` (non-spec demos) are `pub(crate)`.
- **`#![forbid(unsafe_code)]`** on FN-DSA — no `unsafe` block can be added without a deliberate lint downgrade.
- **Constant-Time Operations**: FN-DSA implements FIPS 206 Algorithm 12 SamplerZ with a constant-time half-Gaussian CDT base sampler at `sigma_0 = 1.8205`, arithmetic-mask comparisons, and a constant-time FACCT polynomial BerExp acceptance test (Algorithm 14); the keygen `is_invertible` check uses a branchless accumulator; the FFT-tree `d1` clamp is branchless. UOV's GF(2^8) multiplication and matrix operations use arithmetic masking. ML-DSA's field arithmetic uses Barrett reduction with branchless correction.
- **Input Validation**: All deserialization functions validate coefficient ranges and structural integrity. FN-DSA's `decode_secret_key` checks the NTRU equation `f*G - g*F = q`; `verify` validates `pk.h` length and coefficient range.
- **Cryptographic RNG**: FN-DSA `sign` and `keygen` require `R: RngCore + CryptoRng` at the type level, rejecting non-CSPRNG sources at compile time.
- **Error Handling**: Proper `Result` types instead of panics for invalid parameters.
- **KAT Validation**: FN-DSA (FALCON-512) passes all 100 NIST PQC Known Answer Test vectors.
- **Cross-implementation testing**: Wire-format round-trip plus negative tests (cross-keypair, tampered signature, tampered message) at `tests/interop_falcon512.rs`. *Open gap:* sign-with-ours / verify-with-PQClean FFI roundtrip is deferred (license and Windows toolchain considerations).
- **Supply chain**: Pinned dependencies (`=x.y.z`), `deny.toml` at the workspace root, MSRV declared.
- **Test coverage**: FN-DSA at 126/126 passing (117 lib + 5 KAT + 4 interop) with no `#[ignore]`; a statistical norm-distribution test catches sampler regressions on the first sample.

## Implemented Schemes

### ML-DSA (FIPS 204)

Module-Lattice Digital Signature Algorithm, the NIST post-quantum signature standard (formerly CRYSTALS-Dilithium).

**Status**: FIPS 204 finalized (August 2024).

**Security basis**: Module Learning With Errors (M-LWE) and Module Short Integer Solution (M-SIS) problems.

| Parameter Set | Security Level | Public Key | Secret Key | Signature |
|---------------|----------------|------------|------------|-----------|
| ML-DSA-44     | NIST Level 2   | 1,312 B    | 2,560 B    | 2,420 B   |
| ML-DSA-65     | NIST Level 3   | 1,952 B    | 4,032 B    | 3,309 B   |
| ML-DSA-87     | NIST Level 5   | 2,592 B    | 4,896 B    | 4,627 B   |

```rust
use pqsigs_ml_dsa::{keygen, sign, verify, params::ML_DSA_44};
use rand::rngs::OsRng;

let (pk, sk) = keygen(&mut OsRng, ML_DSA_44);
let sig = sign(&sk, b"message").unwrap();
assert!(verify(&pk, b"message", &sig).is_ok());
```

### FN-DSA (FALCON)

Fast Fourier Lattice-based Compact Signatures over NTRU, a NIST post-quantum signature standard.

> **DO NOT USE IN PRODUCTION.** This implementation is functionally
> correct (all 100 NIST KAT vectors verify, our own signatures
> round-trip our own verifier, wire format round-trips encode/decode
> cleanly) and implements the FIPS 206 Algorithm 12 SamplerZ with a
> constant-time half-Gaussian base sampler and a constant-time FACCT
> polynomial BerExp acceptance test. Hash-to-point FIPS 206 domain
> separation, MXCSR FTZ/DAZ, and PQClean cross-verifier
> interoperability are deferred. See [SECURITY.md](SECURITY.md) for
> the full picture.

**Status**: FIPS 206 draft submitted to NIST for approval (August 2025). Initial Public Draft expected; final standard anticipated late 2026 / early 2027. Our implementation passes all 100 Falcon-512 KAT vectors from the NIST PQC reference implementation, validating correctness of the verification path (public key decoding, signature decompression, hash-to-point, NTT multiplication, and norm bound checking). Our own signing path (FIPS 206 Algorithm 11 ffSampling + Algorithm 12 SamplerZ) produces FALCON-distributed signatures that pass our verifier; a statistical norm-distribution canary in `src/sign.rs` guards against future sampler regressions.

**Security basis**: Hardness of the Short Integer Solution (SIS) problem over NTRU lattices.

| Parameter Set | Security Level | Public Key | Secret Key | Signature |
|---------------|----------------|------------|------------|-----------|
| FALCON-512    | NIST Level 1   | 897 B      | 1,281 B    | ~666 B    |
| FALCON-1024   | NIST Level 5   | 1,793 B    | 2,305 B    | ~1,280 B  |

```rust
use pqsigs_fn_dsa::{keygen_512, sign, verify};
use rand::rngs::OsRng;

let keypair = keygen_512(&mut OsRng).unwrap();
let sig = sign(&mut OsRng, &keypair.sk, b"message").unwrap();
assert!(verify(&keypair.pk, b"message", &sig).is_ok());
```

### UOV (Unbalanced Oil and Vinegar)

Multivariate quadratic signature scheme, a candidate in the NIST Additional Digital Signatures competition.

**Security basis**: Hardness of solving systems of multivariate quadratic equations over finite fields.

| Parameter Set | Security Level | Public Key | Signature |
|---------------|----------------|------------|-----------|
| NIST L1       | ~128-bit       | ~50 KB     | 156 B     |
| NIST L3       | ~192-bit       | ~240 KB    | 340 B     |
| NIST L5       | ~256-bit       | ~660 KB    | 564 B     |

```rust
use pqsigs_uov::{keygen, sign, verify, params::PARAMS_NIST_L1};
use rand::rngs::OsRng;

let (pk, sk) = keygen(&mut OsRng, PARAMS_NIST_L1);
let sig = sign(&mut OsRng, &pk, &sk, b"message").unwrap();
assert!(verify(&pk, b"message", &sig).is_ok());
```

## Project Structure

```
pq-crypto/
├── schemes/
│   ├── ml-dsa/          # ML-DSA (FIPS 204) implementation
│   │   └── src/
│   │       ├── ntt.rs       # Number Theoretic Transform
│   │       ├── poly.rs      # Polynomial arithmetic
│   │       ├── sampling.rs  # SHAKE-based sampling
│   │       ├── rounding.rs  # Power2Round, Decompose, Hints
│   │       ├── keygen.rs    # Key generation
│   │       ├── sign.rs      # Signing with Fiat-Shamir
│   │       └── verify.rs    # Verification
│   │
│   ├── fn-dsa/          # FN-DSA (FALCON) implementation
│   │   ├── src/
│   │   │   ├── fft.rs       # FFT for polynomial ring arithmetic
│   │   │   ├── ntru.rs      # NTRUSolve algorithm
│   │   │   ├── poly.rs      # Polynomial operations with NTT (mul_ntt) and FP-FFT (deprecated)
│   │   │   ├── hash.rs      # SHAKE256 hash-to-point (matches Falcon reference)
│   │   │   ├── gaussian.rs  # FIPS 206 Algorithm 12 SamplerZ + half-Gaussian CDT base sampler
│   │   │   ├── sampler.rs   # FIPS 206 Algorithm 11 recursive ffSampling
│   │   │   ├── packing.rs   # Key/signature serialization (FIPS 206 + NIST KAT formats)
│   │   │   ├── keygen.rs    # Key generation with NTRUGen + LDL* quality check
│   │   │   ├── sign.rs      # Signature generation
│   │   │   └── verify.rs    # Verification (with pk validation)
│   │   └── tests/
│   │       ├── kat_falcon512.rs      # NIST KAT verification tests (100 vectors)
│   │       ├── interop_falcon512.rs  # Wire-format round-trip + forgery resistance
│   │       └── falcon512-KAT.rsp     # Official Falcon-512 KAT vectors
│   │
│   └── uov/             # UOV implementation
│       └── src/
│           ├── field.rs     # GF(2^8) arithmetic
│           ├── matrix.rs    # Matrix operations
│           ├── keygen.rs    # Key generation
│           ├── sign.rs      # Signing
│           └── verify.rs    # Verification
└── README.md
```

## Building and Testing

Each scheme is a separate Cargo workspace member:

```bash
# Build all schemes (run each from the repo root)
(cd schemes/ml-dsa && cargo build)
(cd schemes/fn-dsa && cargo build)
(cd schemes/uov && cargo build)

# Run tests (run each from the repo root)
(cd schemes/ml-dsa && cargo test)
(cd schemes/fn-dsa && cargo test)
(cd schemes/uov && cargo test)

# Run FALCON-512 KAT verification (100 vectors)
(cd schemes/fn-dsa && cargo test --release --test kat_falcon512)

# Run a specific scheme's tests with optimizations
(cd schemes/fn-dsa && cargo test --release)
```

## Test Coverage

**Total: 333 tests passing**

| Scheme | Tests | Coverage |
|--------|-------|----------|
| ML-DSA | 130 | All parameter sets (44/65/87), NTT operations, polynomial arithmetic, roundtrip verification, tampering detection, packing/unpacking, edge cases |
| FN-DSA | 126 | 117 lib + 5 KAT + 4 interop. FFT operations, NTRUSolve (n=2 to 512), FIPS 206 Algorithm 12 SamplerZ distribution check, NTT polynomial arithmetic, sign/verify roundtrip (n=16, 512), statistical norm-distribution canary, wire-format round-trip, cross-keypair / tampered-signature / tampered-message rejection, NIST Falcon-512 KAT verification (100/100 vectors passing) |
| UOV | 77 | 27 unit + 46 integration + 4 doc tests covering all parameter sets, GF(2^8) field arithmetic, matrix operations, tampering detection, cross-key rejection, stress tests |

## Implementation Notes

### ML-DSA
- Reference implementation following FIPS 204 specification
- Uses plain modular arithmetic NTT (not Montgomery form)
- Deterministic signing using SHAKE256
- Implements the "Fiat-Shamir with Aborts" paradigm

### FN-DSA (FALCON)
- Implements recursive NTRUSolve algorithm for finding F, G such that `fG - gF = q`
- Uses FFT/iFFT for fast polynomial arithmetic in the negacyclic ring `Z[X]/(X^n + 1)`
- Exact modular arithmetic via NTT for both verification and the signer-side norm check (`Poly::mul_ntt`); the legacy FP-FFT `Poly::mul` is retained but `#[deprecated]`
- Modular arithmetic over `Z_q` with `q = 12289` (NTT-friendly prime)
- FIPS 206 Algorithm 11 recursive ffSampling over the LDL\* tree; FIPS 206 Algorithm 12 SamplerZ for the per-leaf integer Gaussian draw with a constant-time half-Gaussian CDT base sampler at `sigma_0 = 1.8205` and BerExp acceptance
- SHAKE256 hash-to-point matching the Falcon reference implementation
- Supports both FIPS 206 wire format and original Falcon NIST API format
- Verified against all 100 official Falcon-512 NIST PQC KAT vectors
- `#![forbid(unsafe_code)]`, pinned dependencies, `deny.toml`, declared MSRV
- **Open gaps** (see [SECURITY.md](SECURITY.md)):
   - Hash-to-point domain-separation prefix deferred until FIPS 206 final standard locks the format (we match PQClean today)
   - MXCSR FTZ/DAZ left at platform default to keep `forbid(unsafe_code)` in force
   - PQClean cross-verifier FFI roundtrip substituted with in-tree wire-format roundtrip
- **Note**: Keygen and signing use our own randomness (not the NIST AES-CTR-DRBG), so byte-for-byte reproducibility of reference keygen/signatures is not possible. KAT testing is verification-only.
- **Note**: FIPS 206 is not yet finalized. Once official FIPS 206 KAT vectors are published, they will be integrated separately. The current KAT vectors are from the original Falcon NIST PQC submission, which uses a slightly different wire format (header bytes, compression scheme).

### UOV
- Implements GF(2^8) field arithmetic with irreducible polynomial x^8 + x^4 + x^3 + x + 1
- Uses SHAKE256 for deterministic key expansion
- Includes demo parameters for fast testing

## References

- [FIPS 204: ML-DSA Standard](https://csrc.nist.gov/pubs/fips/204/final)
- [FIPS 205: SLH-DSA Standard](https://csrc.nist.gov/pubs/fips/205/final)
- [FIPS 206: FN-DSA (FALCON)](https://csrc.nist.gov/presentations/2025/fips-206-fn-dsa-falcon) (draft, not yet finalized)
- [NIST Post-Quantum Cryptography](https://csrc.nist.gov/projects/post-quantum-cryptography)
- [NIST Additional Digital Signatures](https://csrc.nist.gov/projects/pqc-dig-sig)
- [FALCON Specification](https://falcon-sign.info/)

## License

See [LICENSE](LICENSE) for details.
