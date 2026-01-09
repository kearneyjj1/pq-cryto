//! Benchmarks for SQI-SIGN.

use criterion::{criterion_group, criterion_main, Criterion};

fn bench_keygen(_c: &mut Criterion) {
    // TODO: Add keygen benchmark
}

fn bench_sign(_c: &mut Criterion) {
    // TODO: Add signing benchmark
}

fn bench_verify(_c: &mut Criterion) {
    // TODO: Add verification benchmark
}

criterion_group!(benches, bench_keygen, bench_sign, bench_verify);
criterion_main!(benches);
