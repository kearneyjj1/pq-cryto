# Security Policy

## Disclaimer

**This repository contains EDUCATIONAL implementations of post-quantum cryptographic signature schemes. These implementations are NOT intended for production use and have NOT been audited by professional cryptographers.**

## Known Limitations

### Side-Channel Vulnerabilities

All implementations in this repository may be vulnerable to timing side-channel attacks. Specifically:

1. **Variable-Time Operations**: Many operations are not constant-time, including:
   - Field arithmetic operations (ML-DSA, UOV)
   - Polynomial operations
   - Matrix operations (Gaussian elimination, UOV-only)
   - Some Gaussian sampling paths (see FN-DSA notes below)

2. **Rejection Sampling Timing**: ML-DSA signing uses rejection sampling which inherently has variable execution time based on the message and key. The number of loop iterations reveals information about the signing process. This is an accepted trade-off in the standard, but callers should be aware that signing time is not constant. FN-DSA signing also has an outer rejection loop per FALCON Algorithm 17; this is inherent to the scheme and treated as public in the security proof.

3. **Floating-Point Arithmetic**: The FN-DSA (FALCON) implementation uses floating-point arithmetic which has platform-dependent timing characteristics. The implementation assumes IEEE 754 double-precision floating-point support and may produce incorrect results on platforms with non-standard floating-point behavior (e.g., flush-to-zero mode, non-IEEE rounding). MXCSR FTZ/DAZ are intentionally *not* enabled at crate init because the required intrinsic is `unsafe`; see `#![forbid(unsafe_code)]` discussion in `schemes/fn-dsa/src/fft.rs`.

4. **Memory Access Patterns**: Array indexing patterns may leak information through cache timing.

5. **FN-DSA BerExp**: The integer-Gaussian acceptance test (FIPS 206 Algorithm 14, BerExp) is implemented via a 13-coefficient fixed-point Horner polynomial (FACCT, matching PQClean `fpr_expm_p63`), so the per-leaf BerExp acceptance is constant-time over its inner loop. The outer rejection-loop iteration count remains a function of the per-leaf target; this leakage is inherent to FALCON's design and is treated as public in the security proof.

### Implementation Gaps

1. **FN-DSA (FALCON)**: Signing produces FALCON-distributed signatures via FIPS 206 Algorithm 11 (recursive ffSampling) with the per-leaf integer Gaussian draw from Algorithm 12 (SamplerZ), including the constant-time FACCT polynomial BerExp acceptance test (Algorithm 14). All 100 official Falcon-512 KAT vectors verify, and a self-roundtrip distribution test confirms the empirical norm distribution is consistent with the spec. Outstanding deviations from the FIPS 206 final standard (still in draft):
   - Hash-to-point input mix matches PQClean (`SHAKE256(salt || msg)`) for KAT compatibility; the FIPS 206 final header-byte domain separator is not yet locked.
   - Keygen `generate_small_poly` uses a variable-time rejection sampler (keygen-time leak only; the signing path is independent and constant-time).
   - Custom RNG path (not NIST AES-CTR-DRBG); byte-level KAT reproduction is not in scope.

2. **Secret Key Handling**: While zeroization is implemented, it may not be effective against all memory disclosure attacks (e.g., speculative execution, cold boot attacks). Heap pages are not `mlock`ed.

3. **Fault Attacks**: No countermeasures against fault injection attacks are implemented.

4. **No external interoperability test**: We round-trip FN-DSA signatures through our own encode/decode/verify path (see `tests/interop_falcon512.rs`) but do not yet sign with our impl and verify with PQClean's reference C verifier. Vendoring PQClean for FFI testing is deferred (license + Windows toolchain considerations).

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
- **FN-DSA**: `SecretKey` zeroizes `f`, `g`, `F`, `G`, `h`, all FFT-domain copies, the LDL\* tree, and recursive ffSampling intermediates (per-call, ~10k secret-bearing vectors otherwise lingering on the heap)
- **UOV**: `SecretKey` zeroizes transformation matrix `t` and quadratic forms `f_quads`
- **`#![forbid(unsafe_code)]`** on the FN-DSA crate prevents any future `unsafe` block from being added without a deliberate lint downgrade

### Type-System Secret Containment (FN-DSA)
- All secret-material fields on `SecretKey`, `GramSchmidt`, `FftTree`, and `FftNode` are `pub(crate)` so external consumers cannot read raw `(f, g, F, G)` or reconstruct them from FFT-domain copies
- `SecretKey` does not implement `Clone`; duplication requires the explicit `try_clone()` method to keep the secret-key lifecycle visible at the call site
- `sign_simple` and `sign_with_nonce` (which bypass FALCON's proper ffSampling and would produce non-spec signatures) are `pub(crate)` and not part of the published API

### Constant-Time Operations
- **GF(2^8) Multiplication** (UOV): Uses fixed 8-iteration loop with arithmetic masking (no data-dependent branches)
- **Matrix Operations** (UOV): `mat_mul`, `rank`, and `solve` use constant-time pivot selection and elimination
- **FN-DSA Half-Gaussian Base Sampler**: CDT-based linear scan with arithmetic-mask comparisons at `sigma_0 = 1.8205` per FIPS 206 Algorithm 12. The full Algorithm 12 outer rejection loop is implemented; the remaining non-CT operation is `libm::exp` inside BerExp (see Side-Channel #5 above)
- **FN-DSA `is_invertible`**: Replaced early-exit pattern with a branchless accumulator that scans every NTT coefficient regardless of input
- **FN-DSA FFT-tree `d1` clamp**: Branchless `f64::max(0.0)` instead of `if d1 < 0.0` on secret-derived value
- **FN-DSA field arithmetic** (`Zq`): Barrett reduction, arithmetic-mask add/sub, branchless negation, `pow_ct` for exponent-sensitive operations

### Input Validation
- **Coefficient Range Checks**: All unpacking functions validate coefficient bounds
- **Structural Validation**: Deserialization validates lengths, indices, and duplicate detection
- **NTRU Equation Validation** (FN-DSA): `decode_secret_key` verifies `f*G - g*F = q` on the decoded basis polynomials before constructing a `SecretKey`. Prevents acceptance of crafted-blob secret keys that would otherwise sign garbage.
- **Public Key Validation** (FN-DSA): `verify` checks `pk.h.len() == params.n` and `pk.h[i] in [0, q)` before any FFT operation, preventing panics on a directly-constructed malformed `PublicKey`
- **Cryptographic RNG**: `sign`, `keygen`, `keygen_512`, and `keygen_1024` require `R: RngCore + CryptoRng`. The type system rejects weak PRNGs (XorShift, PCG, etc.) at compile time.
- **Error Types**: Granular `VerificationFailure` enums provide debugging info without creating oracles

### Wire-Format Correctness (FN-DSA)
- Sign-side polynomial multiplication uses NTT (`Poly::mul_ntt`) so the signer-side norm check is bitwise-identical to the verifier. The FP-FFT `Poly::mul` is retained but `#[deprecated]` to deter future misuse.
- New `tests/interop_falcon512.rs`: encode → decode round-trip for `(PublicKey, Signature)`, plus negative tests for cross-keypair forgery, tampered signatures, and tampered messages.

### Test Coverage (FN-DSA)
- 117 lib unit tests + 5 NIST KAT integration tests + 4 wire-format / forgery-resistance integration tests = 126/126 passing
- `test_sign_falcon_512_norm_distribution`: signs 50 messages on a fresh FALCON-512 keypair and asserts the empirical norm distribution is within a sensible band of the spec-predicted distribution. This is the regression-canary that would have caught the 2026-05-12 sampler bug on the first sample.

### Supply Chain
- **Dependency pinning**: `Cargo.toml` uses `=x.y.z` exact-version specifiers (committed `Cargo.lock` further fixes the resolution)
- **`deny.toml`**: cargo-deny configuration at the workspace root with allow-listed licenses, denial of unknown registries/git sources, and a yanked-crate gate
- **MSRV declared**: `rust-version = "1.76"` in `Cargo.toml`

### Code Quality
- **No Panics**: Pack functions return `Result` types instead of panicking on invalid parameters
- **Bounds Checking**: Security-critical functions use `assert!` (not `debug_assert!`) for release-mode safety. FN-DSA's leaf Hermitian-symmetry check (tightened to 1e-9) runs in release.
- **Named Constants**: Magic numbers replaced with documented constants (e.g., `TOY_KEYGEN_SIGMA`, `MAX_SIGN_ATTEMPTS`)

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
| 0.3.0 | 2026-02-10 | Updated FN-DSA documentation to reflect full ffSampling implementation and 100/100 KAT vector validation. Fixed broken FIPS 206 reference link. |
| 0.4.0 | 2026-05-13 | **FN-DSA hardening pass.** FIPS 206 Algorithm 12 SamplerZ implemented with a constant-time half-Gaussian CDT base sampler at `sigma_0 = 1.8205` and a constant-time FACCT polynomial BerExp acceptance test (Algorithm 14, matching PQClean's `fpr_expm_p63`). Visibility tightened: secret-material fields are `pub(crate)`, `Clone` removed from `SecretKey` (replaced with explicit `try_clone()`), `sign_simple` / `sign_with_nonce` gated to `pub(crate)`. NTRU equation validated in `decode_secret_key`. `verify` validates public-key length and coefficient range. `CryptoRng` bound added to `sign` and `keygen`. `is_invertible` rewritten as a branchless accumulator. `d1[i] < 0` branch in the FFT tree replaced with branchless `max`. Hermitian-symmetry check at the leaf tightened to `assert!(... < 1e-9)` and promoted to release. `Poly::mul` (FP-FFT) marked `#[deprecated]`. Runtime canary on the `sigma_min` leaf clamp. `#![forbid(unsafe_code)]`, `#![deny(noop_method_call)]`, MSRV declared, dependencies pinned, `deny.toml` added. New tests: signature norm-distribution canary, wire-format round-trip, cross-keypair / tampered-signature / tampered-message / single-bit-flip rejection. Dead-code cleanup in `gaussian.rs` (legacy variable-time samplers and the orphaned `CDT_BITS` constant gated `#[cfg(test)]` or removed). Remaining open: hash-to-point domain separation (deferred until FIPS 206 final standard locks the format), MXCSR FTZ/DAZ (intentional defer under `forbid(unsafe_code)`), PQClean FFI roundtrip (substituted with in-tree wire-format roundtrip), keygen-time `generate_small_poly` sampler (CDT refactor pending). |

## References

- [NIST Post-Quantum Cryptography](https://csrc.nist.gov/projects/post-quantum-cryptography)
- [FIPS 204: ML-DSA Standard](https://csrc.nist.gov/pubs/fips/204/final)
- [FIPS 206: FN-DSA (FALCON)](https://csrc.nist.gov/presentations/2025/fips-206-fn-dsa-falcon) (draft, not yet finalized)
- [Side-Channel Attacks on Implementations of Curve25519](https://cr.yp.to/papers.html)
- [BearSSL Constant-Time Guide](https://bearssl.org/constanttime.html)
- [PQClean Falcon-512 reference C](https://github.com/PQClean/PQClean/tree/master/crypto_sign/falcon-512) — the implementation our wire format is compared against.
