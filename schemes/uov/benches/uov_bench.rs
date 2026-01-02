//! Benchmarks for the UOV signature scheme.
//!
//! Run with: cargo bench

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::rngs::OsRng;

use pqsigs_uov::{
    field::F,
    keygen::keygen,
    params::{PARAMS_DEMO, PARAMS_NIST_L1},
    sign::sign,
    verify::verify,
};

/// Benchmark field operations.
fn bench_field_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("field");

    let a = F(0x57);
    let b = F(0x83);

    group.bench_function("add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
    });

    group.bench_function("mul", |bencher| {
        bencher.iter(|| black_box(a) * black_box(b))
    });

    group.bench_function("inverse", |bencher| bencher.iter(|| black_box(a).inverse()));

    group.finish();
}

/// Benchmark key generation.
fn bench_keygen(c: &mut Criterion) {
    let mut group = c.benchmark_group("keygen");

    group.bench_function("demo", |bencher| {
        bencher.iter(|| keygen(&mut OsRng, black_box(PARAMS_DEMO)))
    });

    // NIST L1 is slower, so we use fewer samples
    group.sample_size(10);
    group.bench_function("nist_l1", |bencher| {
        bencher.iter(|| keygen(&mut OsRng, black_box(PARAMS_NIST_L1)))
    });

    group.finish();
}

/// Benchmark signing.
fn bench_sign(c: &mut Criterion) {
    let mut group = c.benchmark_group("sign");

    // Demo params
    let (pk_demo, sk_demo) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"benchmark message for signing";

    group.bench_function("demo", |bencher| {
        bencher.iter(|| sign(&mut OsRng, &pk_demo, &sk_demo, black_box(msg)))
    });

    // NIST L1 params
    group.sample_size(10);
    let (pk_l1, sk_l1) = keygen(&mut OsRng, PARAMS_NIST_L1);

    group.bench_function("nist_l1", |bencher| {
        bencher.iter(|| sign(&mut OsRng, &pk_l1, &sk_l1, black_box(msg)))
    });

    group.finish();
}

/// Benchmark verification.
fn bench_verify(c: &mut Criterion) {
    let mut group = c.benchmark_group("verify");

    // Demo params
    let (pk_demo, sk_demo) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"benchmark message for verification";
    let sig_demo = sign(&mut OsRng, &pk_demo, &sk_demo, msg).unwrap();

    group.bench_function("demo", |bencher| {
        bencher.iter(|| verify(&pk_demo, black_box(msg), &sig_demo))
    });

    // NIST L1 params
    group.sample_size(10);
    let (pk_l1, sk_l1) = keygen(&mut OsRng, PARAMS_NIST_L1);
    let sig_l1 = sign(&mut OsRng, &pk_l1, &sk_l1, msg).unwrap();

    group.bench_function("nist_l1", |bencher| {
        bencher.iter(|| verify(&pk_l1, black_box(msg), &sig_l1))
    });

    group.finish();
}

/// Benchmark full sign+verify cycle.
fn bench_full_cycle(c: &mut Criterion) {
    let mut group = c.benchmark_group("full_cycle");

    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);
    let msg = b"benchmark message";

    group.bench_function("demo", |bencher| {
        bencher.iter(|| {
            let sig = sign(&mut OsRng, &pk, &sk, black_box(msg)).unwrap();
            verify(&pk, msg, &sig)
        })
    });

    group.finish();
}

/// Benchmark with varying message sizes.
fn bench_message_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("message_size");

    let (pk, sk) = keygen(&mut OsRng, PARAMS_DEMO);

    for size in [64, 256, 1024, 4096].iter() {
        let msg = vec![0xABu8; *size];

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &msg, |bencher, msg| {
            bencher.iter(|| sign(&mut OsRng, &pk, &sk, black_box(msg)))
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_field_ops,
    bench_keygen,
    bench_sign,
    bench_verify,
    bench_full_cycle,
    bench_message_sizes,
);
criterion_main!(benches);
