use binary_fields::GF2_128;
use binary_fields::reduce::{reduce_2_40, reduce_2_41};
use criterion::{Criterion, black_box, criterion_group, criterion_main};

fn bench_add_2_32(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);
    let b = GF2_128::new(0x1234567890abcdef, 0xfedcba9876543210);

    c.bench_function("add_2_32", |bencher| {
        bencher.iter(|| black_box(a).add_2_32(black_box(b)))
    });
}

fn bench_mul_2_33(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);
    let b = GF2_128::new(0x1234567890abcdef, 0xfedcba9876543210);

    c.bench_function("mul_2_33", |bencher| {
        bencher.iter(|| black_box(a).mul_2_33(black_box(b)))
    });
}

fn bench_mul_2_34(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);
    let b = GF2_128::new(0x1234567890abcdef, 0xfedcba9876543210);

    c.bench_function("mul_2_34", |bencher| {
        bencher.iter(|| black_box(a).mul_2_34(black_box(b)))
    });
}

fn bench_square_2_39(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);

    c.bench_function("square_2_39", |bencher| {
        bencher.iter(|| black_box(a).square_2_39())
    });
}

fn bench_invert_2_48(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);

    c.bench_function("invert_2_48", |bencher| {
        bencher.iter(|| black_box(a).invert_2_48())
    });
}

fn bench_invert_2_49(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);

    c.bench_function("invert_2_49", |bencher| {
        bencher.iter(|| black_box(a).invert_2_49())
    });
}

fn bench_reduce_2_40(c: &mut Criterion) {
    let input = [0xdeadbeefcafe1234u64, 0x0102030405060708, 0xfedcba9876543210, 0x1111222233334444];

    c.bench_function("reduce_2_40", |bencher| {
        bencher.iter(|| reduce_2_40(black_box(input)))
    });
}

fn bench_reduce_2_41(c: &mut Criterion) {
    let input = [0xdeadbeefcafe1234u64, 0x0102030405060708, 0xfedcba9876543210, 0x1111222233334444];

    c.bench_function("reduce_2_41", |bencher| {
        bencher.iter(|| reduce_2_41(black_box(input)))
    });
}

fn bench_pow_rtl_2_34_2_39(c: &mut Criterion) {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);

    c.bench_function("pow_rtl_2_34_2_39_large", |bencher| {
        bencher.iter(|| black_box(a).pow_rtl_2_34_2_39(black_box(u128::MAX - 1)))
    });
}

criterion_group!(
    benches,
    bench_add_2_32,
    bench_reduce_2_40,
    bench_reduce_2_41,
    bench_mul_2_33,
    bench_mul_2_34,
    bench_square_2_39,
    bench_invert_2_48,
    bench_invert_2_49,
    bench_pow_rtl_2_34_2_39,
);
criterion_main!(benches);
