use binary_fields::generic::arithmetic::mul_2_33;
use binary_fields::generic::invert::invert_2_48;
use binary_fields::fields::z128_z7_z2_z1::{square_2_39, pow_rtl_2_34_2_39};

const GCM_POLY: [u64; 2] = [0x87, 0];
const DEGREE: usize = 128;

#[test]
fn mul_known_vector_1() {
    assert_eq!(mul_2_33(&[0x1234, 0], &[0x5678, 0], &GCM_POLY, DEGREE), [0x5c58160, 0]);
}

#[test]
fn mul_known_vector_2() {
    assert_eq!(
        mul_2_33(&[0xdeadbeefcafe1234, 0], &[0xabcd1234, 0], &GCM_POLY, DEGREE),
        [0xe958a400e6980510, 0x7348531a],
    );
}

#[test]
fn mul_known_vector_3() {
    assert_eq!(mul_2_33(&[1, 0], &[2, 0], &GCM_POLY, DEGREE), [2, 0]);
}

#[test]
fn mul_known_vector_4() {
    assert_eq!(
        mul_2_33(&[0xffffffffffffffff, 0xffffffffffffffff], &[1, 0], &GCM_POLY, DEGREE),
        [0xffffffffffffffff, 0xffffffffffffffff],
    );
}

#[test]
fn inv_known_vector_1() {
    assert_eq!(
        invert_2_48(&[0x1234, 0], &GCM_POLY, DEGREE),
        Some([0xed291313806d0de5, 0x1a1bd097720be28e]),
    );
}

#[test]
fn inv_known_vector_2() {
    assert_eq!(
        invert_2_48(&[0xdeadbeefcafe1234, 0], &GCM_POLY, DEGREE),
        Some([0xc996631b509a63ad, 0xedc7fb439d22b363]),
    );
}

#[test]
fn inv_known_vector_3() {
    assert_eq!(invert_2_48(&[1, 0], &GCM_POLY, DEGREE), Some([1, 0]));
}

#[test]
fn inv_known_vector_4() {
    assert_eq!(
        invert_2_48(&[0xffffffffffffffff, 0xffffffffffffffff], &GCM_POLY, DEGREE),
        Some([0xfc10c53d1c96eca8, 0xfe08629e8e4b766a]),
    );
}

#[test]
fn pow_known_vector_0() {
    assert_eq!(pow_rtl_2_34_2_39([0xdeadbeefcafe1234, 0], 0), [1, 0]);
}

#[test]
fn pow_known_vector_1() {
    let a = [0xdeadbeefcafe1234_u64, 0];
    assert_eq!(pow_rtl_2_34_2_39(a, 1), a);
}

#[test]
fn pow_known_vector_3() {
    assert_eq!(
        pow_rtl_2_34_2_39([0xdeadbeefcafe1234, 0], 3),
        [0xce2b4a77e44ec116, 0xe555280bb487b18a],
    );
}

#[test]
fn pow_known_vector_17() {
    assert_eq!(
        pow_rtl_2_34_2_39([0xdeadbeefcafe1234, 0], 17),
        [0xb8b59dbbe13f7c2a, 0xd4b6b856a3200f2c],
    );
}

#[test]
fn pow_known_vector_5() {
    assert_eq!(
        pow_rtl_2_34_2_39([0x1234567890abcdef, 0x1234567890abcdef], 5),
        [0x21d292a413a4f9ef, 0x242c45156eccd574],
    );
}

#[test]
fn square_known_vector_1() {
    assert_eq!(
        square_2_39([0xdeadbeefcafe1234, 0]),
        [0x5044555401040510, 0x5154445145545455],
    );
}

#[test]
fn square_known_vector_2() {
    assert_eq!(
        square_2_39([0x1234567890abcdef, 0x1234567890abcdef]),
        [0x0623bb37c94dd37e, 0x841a9668ec72dfa1],
    );
}
