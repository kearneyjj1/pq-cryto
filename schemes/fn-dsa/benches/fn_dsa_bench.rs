//! Benchmarks for FN-DSA (FALCON).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pqsigs_fn_dsa::fft::{fft, ifft, Complex};
use pqsigs_fn_dsa::poly::Poly;

fn bench_fft(c: &mut Criterion) {
    let mut group = c.benchmark_group("fft");

    // FFT-512
    let data_512: Vec<Complex> = (0..512)
        .map(|i| Complex::new(i as f64, 0.0))
        .collect();

    group.bench_function("fft_512", |b| {
        b.iter(|| {
            let mut d = data_512.clone();
            fft(black_box(&mut d));
            d
        })
    });

    group.bench_function("ifft_512", |b| {
        b.iter(|| {
            let mut d = data_512.clone();
            ifft(black_box(&mut d));
            d
        })
    });

    // FFT-1024
    let data_1024: Vec<Complex> = (0..1024)
        .map(|i| Complex::new(i as f64, 0.0))
        .collect();

    group.bench_function("fft_1024", |b| {
        b.iter(|| {
            let mut d = data_1024.clone();
            fft(black_box(&mut d));
            d
        })
    });

    group.finish();
}

fn bench_poly_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_mul");

    // Polynomial multiplication for n=512
    let a_512: Vec<i16> = (0..512).map(|i| (i % 100) as i16).collect();
    let b_512: Vec<i16> = (0..512).map(|i| ((i + 50) % 100) as i16).collect();

    let a_poly_512 = Poly::from_i16(&a_512);
    let b_poly_512 = Poly::from_i16(&b_512);

    group.bench_function("poly_mul_512", |b| {
        b.iter(|| {
            let a = a_poly_512.clone();
            let b = b_poly_512.clone();
            black_box(a.mul(&b))
        })
    });

    // Polynomial multiplication for n=1024
    let a_1024: Vec<i16> = (0..1024).map(|i| (i % 100) as i16).collect();
    let b_1024: Vec<i16> = (0..1024).map(|i| ((i + 50) % 100) as i16).collect();

    let a_poly_1024 = Poly::from_i16(&a_1024);
    let b_poly_1024 = Poly::from_i16(&b_1024);

    group.bench_function("poly_mul_1024", |b| {
        b.iter(|| {
            let a = a_poly_1024.clone();
            let b = b_poly_1024.clone();
            black_box(a.mul(&b))
        })
    });

    group.finish();
}

fn bench_hash(c: &mut Criterion) {
    use pqsigs_fn_dsa::hash::hash_to_point;
    use pqsigs_fn_dsa::params::{FALCON_512, NONCE_SIZE};

    let mut group = c.benchmark_group("hash");

    let message = b"This is a test message for FALCON signing.";
    let nonce = [0u8; NONCE_SIZE];

    group.bench_function("hash_to_point_512", |b| {
        b.iter(|| {
            black_box(hash_to_point(message, &nonce, &FALCON_512))
        })
    });

    group.finish();
}

criterion_group!(benches, bench_fft, bench_poly_mul, bench_hash);
criterion_main!(benches);
