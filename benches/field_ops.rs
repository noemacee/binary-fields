use binary_fields::generic::arithmetic::{add_2_32, mul_2_33};
use binary_fields::generic::invert::{invert_2_48, invert_2_49};
use binary_fields::generic::reduce::reduce_2_40;
use binary_fields::fields::z128_z7_z2_z1::{mul_2_34, square_2_39, reduce_2_41, pow_rtl_2_34_2_39};
use binary_fields::ark::{Gf128, Gf128Generic};
use ark_ff::Field;
use criterion::{Criterion, black_box, criterion_group, criterion_main};

const GCM_POLY: [u64; 2] = [0x87, 0];
const DEGREE: usize = 128;

fn bench_add_2_32(c: &mut Criterion) {
    let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
    let b = [0x1234567890abcdef_u64, 0xfedcba9876543210];
    c.bench_function("add_2_32", |bencher| {
        bencher.iter(|| add_2_32(black_box(a), black_box(b)))
    });
}

fn bench_mul_2_33(c: &mut Criterion) {
    let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
    let b = [0x1234567890abcdef_u64, 0xfedcba9876543210];
    c.bench_function("mul_2_33", |bencher| {
        bencher.iter(|| mul_2_33(black_box(&a), black_box(&b), &GCM_POLY, DEGREE))
    });
}

fn bench_mul_2_34(c: &mut Criterion) {
    let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
    let b = [0x1234567890abcdef_u64, 0xfedcba9876543210];
    c.bench_function("mul_2_34", |bencher| {
        bencher.iter(|| mul_2_34(black_box(a), black_box(b)))
    });
}

fn bench_square_via_mul_2_33(c: &mut Criterion) {
    let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
    c.bench_function("square_via_mul_2_33", |bencher| {
        bencher.iter(|| mul_2_33(black_box(&a), black_box(&a), &GCM_POLY, DEGREE))
    });
}

fn bench_square_2_39(c: &mut Criterion) {
    let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
    c.bench_function("square_2_39", |bencher| {
        bencher.iter(|| square_2_39(black_box(a)))
    });
}

fn bench_invert_2_48(c: &mut Criterion) {
    let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
    c.bench_function("invert_2_48", |bencher| {
        bencher.iter(|| invert_2_48(black_box(&a), &GCM_POLY, DEGREE))
    });
}

fn bench_invert_2_49(c: &mut Criterion) {
    let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
    c.bench_function("invert_2_49", |bencher| {
        bencher.iter(|| invert_2_49(black_box(&a), &GCM_POLY, DEGREE))
    });
}

fn bench_reduce_2_40(c: &mut Criterion) {
    let input = [0xdeadbeefcafe1234_u64, 0x0102030405060708, 0xfedcba9876543210, 0x1111222233334444];
    c.bench_function("reduce_2_40", |bencher| {
        bencher.iter(|| {
            let mut v = black_box(input).to_vec();
            v.resize(8, 0);
            reduce_2_40(&mut v, &GCM_POLY, DEGREE);
            v[0]
        })
    });
}

fn bench_reduce_2_41(c: &mut Criterion) {
    let input = [0xdeadbeefcafe1234_u64, 0x0102030405060708, 0xfedcba9876543210, 0x1111222233334444];
    c.bench_function("reduce_2_41", |bencher| {
        bencher.iter(|| reduce_2_41(black_box(input)))
    });
}

fn bench_pow_rtl_2_34_2_39(c: &mut Criterion) {
    let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
    c.bench_function("pow_rtl_2_34_2_39_large", |bencher| {
        bencher.iter(|| pow_rtl_2_34_2_39(black_box(a), black_box(u128::MAX - 1)))
    });
}

fn bench_gf128_field(c: &mut Criterion) {
    let a_limbs = [0xdeadbeefcafe1234_u64, 0xabcd123400000000];
    let b_limbs = [0x1234567890abcdef_u64, 0xfedcba9876543210];

    let a_gen = Gf128Generic::new(a_limbs);
    let b_gen = Gf128Generic::new(b_limbs);
    let a_opt = Gf128::new(a_limbs);
    let b_opt = Gf128::new(b_limbs);

    let mut group = c.benchmark_group("gf128");

    group.bench_function("mul/generic",   |bencher| bencher.iter(|| black_box(a_gen) * black_box(b_gen)));
    group.bench_function("mul/optimized", |bencher| bencher.iter(|| black_box(a_opt) * black_box(b_opt)));

    group.bench_function("square/generic",   |bencher| bencher.iter(|| black_box(a_gen).square()));
    group.bench_function("square/optimized", |bencher| bencher.iter(|| black_box(a_opt).square()));

    group.bench_function("inverse/generic",   |bencher| bencher.iter(|| black_box(a_gen).inverse()));
    group.bench_function("inverse/optimized", |bencher| bencher.iter(|| black_box(a_opt).inverse()));

    group.finish();
}

criterion_group!(
    benches,
    bench_add_2_32,
    bench_reduce_2_40,
    bench_reduce_2_41,
    bench_mul_2_33,
    bench_mul_2_34,
    bench_square_via_mul_2_33,
    bench_square_2_39,
    bench_invert_2_48,
    bench_invert_2_49,
    bench_pow_rtl_2_34_2_39,
    bench_gf128_field,
);
criterion_main!(benches);
