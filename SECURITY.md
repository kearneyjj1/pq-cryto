# Security Policy

## Disclaimer

**This repository contains EDUCATIONAL implementations of post-quantum cryptographic signature schemes. These implementations are NOT intended for production use and have NOT been audited by professional cryptographers.**

## Known Limitations

### Side-Channel Vulnerabilities

All implementations in this repository may be vulnerable to timing side-channel attacks. Specifically:

1. **Variable-Time Operations**: Many operations are not constant-time, including:
   - Field arithmetic operations
   - Polynomial operations
   - Matrix operations (Gaussian elimination)
   - Gaussian sampling

2. **Rejection Sampling Timing**: ML-DSA signing uses rejection sampling which inherently has variable execution time based on the message and key. The number of loop iterations reveals information about the signing process. This is an accepted trade-off in the standard, but callers should be aware that signing time is not constant.

3. **Floating-Point Arithmetic**: The FN-DSA (FALCON) implementation uses floating-point arithmetic which has platform-dependent timing characteristics. The implementation assumes IEEE 754 double-precision floating-point support and may produce incorrect results on platforms with non-standard floating-point behavior (e.g., flush-to-zero mode, non-IEEE rounding).

4. **Memory Access Patterns**: Array indexing patterns may leak information through cache timing.

### Implementation Gaps

1. **FN-DSA (FALCON)**: Uses simplified sampling instead of the full `ffSampling` algorithm. Signature bounds are relaxed for educational purposes.

2. **Secret Key Handling**: While zeroization is implemented, it may not be effective against all memory disclosure attacks (e.g., speculative execution, cold boot attacks).

3. **Fault Attacks**: No countermeasures against fault injection attacks are implemented.

### Platform Requirements

The FN-DSA (FALCON) implementation has specific platform requirements:

1. **IEEE 754 Floating-Point**: Requires IEEE 754 double-precision (64-bit) floating-point arithmetic with proper rounding. May produce incorrect results on platforms with:
   - Flush-to-zero (FTZ) mode enabled
   - Denormals-are-zero (DAZ) mode enabled
   - Non-default rounding modes
   - Extended precision floating-point (x87 80-bit mode)

2. **Platform Verification**: If deploying on embedded systems or unusual architectures, verify floating-point behavior first. The `f64` operations must match IEEE 754 semantics.

3. **x86/x64**: Works correctly on modern x86/x64 with SSE2 (default for 64-bit Rust).

4. **ARM**: Works correctly on ARMv7+ and AArch64 with standard FP configuration.

### What This Means

- **DO NOT** use these implementations to protect real data
- **DO NOT** deploy these implementations in any security-sensitive context
- **DO NOT** assume the cryptographic security proofs apply to these implementations

## Security Measures Implemented

The following security best practices have been implemented:

### Memory Safety
- **Secret Key Zeroization**: All secret key types implement `Drop` trait with secure zeroization using the `zeroize` crate
- **ML-DSA**: `SecretKey` zeroizes `rho`, `k`, `tr`, `s1`, `s2`, and `t0` polynomials
- **FN-DSA**: `SecretKey` zeroizes `f`, `g`, `F`, `G`, `h`, and all FFT data
- **UOV**: `SecretKey` zeroizes transformation matrix `t` and quadratic forms `f_quads`

### Constant-Time Operations
- **GF(2^8) Multiplication**: Uses fixed 8-iteration loop with arithmetic masking (no data-dependent branches)
- **Matrix Operations**: `mat_mul`, `rank`, and `solve` use constant-time pivot selection and elimination
- **CDT Sampling**: Base Gaussian sampler uses linear scan with arithmetic masking

### Input Validation
- **Coefficient Range Checks**: All unpacking functions validate coefficient bounds
- **Structural Validation**: Deserialization validates lengths, indices, and duplicate detection
- **Error Types**: Granular `VerificationFailure` enums provide debugging info without creating oracles

### Code Quality
- **No Panics**: Pack functions return `Result` types instead of panicking on invalid parameters
- **Bounds Checking**: Security-critical functions use `assert!` (not `debug_assert!`) for release-mode safety
- **Named Constants**: Magic numbers replaced with documented constants

## Security Guarantees

This implementation provides:

- Functionally correct cryptographic operations (signatures verify correctly)
- Educational demonstration of post-quantum signature algorithms
- A starting point for learning about lattice-based and multivariate cryptography
- Secure memory handling for secret key material
- Constant-time field and matrix operations where implemented

This implementation does NOT provide:

- Full protection against all timing side-channel attacks (rejection sampling timing varies)
- Protection against power analysis or electromagnetic emanation attacks
- Protection against fault injection attacks
- Guaranteed constant-time execution for all operations
- Memory-safe handling of secrets against all attack vectors (e.g., speculative execution)

## Threat Model

These implementations assume:

- The attacker cannot measure execution time
- The attacker cannot observe power consumption
- The attacker cannot inject faults
- The attacker cannot access process memory
- The system is trusted and not compromised

If any of these assumptions are violated, secret key material may be recoverable.

## Reporting Security Issues

Since this is an educational project, there is no formal security vulnerability process. However, if you discover issues that could improve the educational value of this repository:

1. Open a GitHub issue describing the concern
2. If the issue involves a novel attack, please provide a clear explanation suitable for educational purposes

## Recommendations for Production Use

If you need production-quality post-quantum cryptography, consider:

1. **liboqs** - Open Quantum Safe project (https://openquantumsafe.org/)
2. **pqcrypto** - Rust bindings for reference implementations
3. **AWS-LC** - Amazon's cryptographic library with PQ support
4. **BoringSSL** - Google's cryptographic library with experimental PQ support

These implementations have undergone more rigorous security review and testing.

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 0.1.0 | 2026-01-08 | Initial security documentation |
| 0.2.0 | 2026-01-08 | Security audit remediation: constant-time GF(2^8) multiplication, zeroization for all secret keys, input validation, error handling improvements |
| 0.2.1 | 2026-01-08 | Follow-up audit fixes: Result-based error handling, release-mode bounds checking, named constants, platform requirements documentation, FALCON-512 runtime warnings |
| 0.3.0 | 2026-01-29 | FN-DSA ffSampling implementation: FIPS 206 compliant Gaussian sampler, LDL* tree construction, recursive ffSampling algorithm. Signatures now mathematically correct; working toward standard bounds. |

## References

- [NIST Post-Quantum Cryptography](https://csrc.nist.gov/projects/post-quantum-cryptography)
- [FIPS 204: ML-DSA Standard](https://csrc.nist.gov/pubs/fips/204/final)
- [FIPS 206: FN-DSA Standard](https://csrc.nist.gov/pubs/fips/206/final)
- [Side-Channel Attacks on Implementations of Curve25519](https://cr.yp.to/papers.html)
- [BearSSL Constant-Time Guide](https://bearssl.org/constanttime.html)
