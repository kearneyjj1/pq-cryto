# pq-crypto

Rust implementations of post-quantum digital signature schemes. Built for learning and experimentation, featuring key generation, signing, verification, and comprehensive test suites.

**Warning**: These implementations are for educational purposes only and are NOT intended for production use. See [SECURITY.md](SECURITY.md) for details.

## Security Status

This codebase has undergone security review with the following measures implemented:

- **Secure Memory Handling**: Secret keys implement `Drop` with zeroization via the `zeroize` crate
- **Constant-Time Operations**: GF(2^8) multiplication and matrix operations use constant-time techniques
- **Input Validation**: All deserialization functions validate coefficient ranges and structural integrity
- **Error Handling**: Proper `Result` types instead of panics for invalid parameters
- **Documentation**: Comprehensive security limitations documented in [SECURITY.md](SECURITY.md)

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

### FN-DSA (FALCON)

Fast Fourier Lattice-based Compact Signatures over NTRU, a NIST post-quantum signature standard.

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

**Note**: This is an educational implementation with relaxed signature bounds. Production FALCON requires the full `ffSampling` algorithm for cryptographic security.

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
│   │   └── src/
│   │       ├── fft.rs       # FFT for polynomial ring arithmetic
│   │       ├── ntru.rs      # NTRUSolve algorithm
│   │       ├── poly.rs      # Polynomial operations with NTT
│   │       ├── gaussian.rs  # Discrete Gaussian sampling
│   │       ├── sampler.rs   # FFT sampler for signing
│   │       ├── keygen.rs    # Key generation with NTRU equation
│   │       ├── sign.rs      # Signature generation
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
cd schemes/fn-dsa && cargo build
cd schemes/uov && cargo build

# Run tests
cd schemes/ml-dsa && cargo test
cd schemes/fn-dsa && cargo test
cd schemes/uov && cargo test

# Run FALCON-512 integration test
cd schemes/fn-dsa && cargo run --bin falcon512_test

# Run with optimizations (recommended for benchmarking)
cargo test --release
```

## Test Coverage

**Total: 307 tests passing**

| Scheme | Tests | Coverage |
|--------|-------|----------|
| ML-DSA | 130 | All parameter sets (44/65/87), NTT operations, polynomial arithmetic, roundtrip verification, tampering detection, packing/unpacking, edge cases |
| FN-DSA | 100 | FFT operations, NTRUSolve (n=2 to 512), Gaussian sampling, polynomial arithmetic, sign/verify integration |
| UOV | 77 | 27 unit + 46 integration + 4 doc tests covering all parameter sets, GF(2^8) field arithmetic, matrix operations, tampering detection, cross-key rejection, stress tests |

## Implementation Notes

### ML-DSA
- Reference implementation following FIPS 204 specification
- Uses plain modular arithmetic NTT (not Montgomery form)
- Deterministic signing using SHAKE256
- Implements the "Fiat-Shamir with Aborts" paradigm

### FN-DSA (FALCON)
- Implements recursive NTRUSolve algorithm for finding F, G such that fG - gF = q
- Uses FFT/iFFT for fast polynomial arithmetic in the negacyclic ring Z[X]/(X^n + 1)
- Modular arithmetic over Z_q with q = 12289 (NTT-friendly prime)
- Discrete Gaussian sampling using the ziggurat method
- **Educational note**: Uses relaxed signature bounds; production use requires full `ffSampling`

### UOV
- Implements GF(2^8) field arithmetic with irreducible polynomial x^8 + x^4 + x^3 + x + 1
- Uses SHAKE256 for deterministic key expansion
- Includes demo parameters for fast testing

## References

- [FIPS 204: ML-DSA Standard](https://csrc.nist.gov/pubs/fips/204/final)
- [FIPS 206: FN-DSA (FALCON) Standard](https://csrc.nist.gov/pubs/fips/206/final)
- [NIST Post-Quantum Cryptography](https://csrc.nist.gov/projects/post-quantum-cryptography)
- [NIST Additional Digital Signatures](https://csrc.nist.gov/projects/pqc-dig-sig)
- [FALCON Specification](https://falcon-sign.info/)

## License

See [LICENSE](LICENSE) for details.
