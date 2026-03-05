use binary_fields::GF2_128;
use criterion::{Criterion, black_box, criterion_group, criterion_main};

fn bench_mul(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);
    let b = GF2_128::new(0x1234567890abcdef, 0xfedcba9876543210);

    c.bench_function("mul", |bencher| {
        bencher.iter(|| black_box(a).mul(black_box(b)))
    });
}

fn bench_square(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);

    c.bench_function("square", |bencher| bencher.iter(|| black_box(a).square()));
}

fn bench_invert(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);

    c.bench_function("invert", |bencher| bencher.iter(|| black_box(a).invert()));
}

fn bench_pow(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);

    c.bench_function("pow_large", |bencher| {
        bencher.iter(|| black_box(a).pow(black_box(u128::MAX - 1)))
    });
}

fn bench_mul_comb(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);
    let b = GF2_128::new(0x1234567890abcdef, 0xfedcba9876543210);

    c.bench_function("mul_comb", |bencher| {
        bencher.iter(|| black_box(a).mul_comb(black_box(b)))
    });
}

fn bench_invert_binary(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);

    c.bench_function("invert_binary", |bencher| {
        bencher.iter(|| black_box(a).invert_binary())
    });
}

criterion_group!(
    benches,
    bench_mul,
    bench_square,
    bench_invert,
    bench_pow,
    bench_mul_comb,
    bench_invert_binary
);
criterion_main!(benches);
