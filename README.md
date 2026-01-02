# pq-crypto

Rust implementations of post-quantum digital signature schemes. Built for learning and experimentation, featuring key generation, signing, verification, and comprehensive test suites.

**Warning**: These implementations are for educational purposes only and have not been audited for production use.

## Implemented Schemes

### ML-DSA (FIPS 204)

Module-Lattice Digital Signature Algorithm, the NIST post-quantum signature standard (formerly CRYSTALS-Dilithium).

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
# Build all schemes
cd schemes/ml-dsa && cargo build
cd schemes/uov && cargo build

# Run tests
cd schemes/ml-dsa && cargo test
cd schemes/uov && cargo test

# Run with optimizations (recommended for benchmarking)
cargo test --release
```

## Test Coverage

- **ML-DSA**: 126 tests covering all parameter sets, roundtrip verification, tampering detection, and edge cases
- **UOV**: 20 tests plus 4 doc tests covering key generation, signing, and verification

## Implementation Notes

### ML-DSA
- Reference implementation following FIPS 204 specification
- Uses plain modular arithmetic NTT (not Montgomery form)
- Deterministic signing using SHAKE256
- Implements the "Fiat-Shamir with Aborts" paradigm

### UOV
- Implements GF(2^8) field arithmetic with irreducible polynomial x^8 + x^4 + x^3 + x + 1
- Uses SHAKE256 for deterministic key expansion
- Includes demo parameters for fast testing

## References

- [FIPS 204: ML-DSA Standard](https://csrc.nist.gov/pubs/fips/204/final)
- [NIST Post-Quantum Cryptography](https://csrc.nist.gov/projects/post-quantum-cryptography)
- [NIST Additional Digital Signatures](https://csrc.nist.gov/projects/pqc-dig-sig)

## License

See [LICENSE](LICENSE) for details.
