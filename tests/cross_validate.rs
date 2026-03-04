use binary_fields::GF2_128;

#[test]
fn mul_known_vector_1() {
    let a = GF2_128::new(0x1234, 0);
    let b = GF2_128::new(0x5678, 0);
    let expected = GF2_128::new(0x5c58160, 0);
    assert_eq!(a.mul(b), expected);
}

#[test]
fn mul_known_vector_2() {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0);
    let b = GF2_128::new(0xabcd1234, 0);
    let expected = GF2_128::new(0xe958a400e6980510, 0x7348531a);
    assert_eq!(a.mul(b), expected);
}

#[test]
fn mul_known_vector_3() {
    let a = GF2_128::new(0x1, 0);
    let b = GF2_128::new(0x2, 0);
    let expected = GF2_128::new(0x2, 0);
    assert_eq!(a.mul(b), expected);
}

#[test]
fn mul_known_vector_4() {
    let a = GF2_128::new(0xffffffffffffffff, 0xffffffffffffffff);
    let b = GF2_128::new(0x1, 0);
    let expected = GF2_128::new(0xffffffffffffffff, 0xffffffffffffffff);
    assert_eq!(a.mul(b), expected);
}

#[test]
fn inv_known_vector_1() {
    let a = GF2_128::new(0x1234, 0);
    let expected = GF2_128::new(0xed291313806d0de5, 0x1a1bd097720be28e);
    assert_eq!(a.invert().unwrap(), expected);
}

#[test]
fn inv_known_vector_2() {
    let a = GF2_128::new(0xdeadbeefcafe1234, 0);
    let expected = GF2_128::new(0xc996631b509a63ad, 0xedc7fb439d22b363);
    assert_eq!(a.invert().unwrap(), expected);
}

#[test]
fn inv_known_vector_3() {
    let a = GF2_128::new(0x1, 0);
    assert_eq!(a.invert().unwrap(), GF2_128::one());
}

#[test]
fn inv_known_vector_4() {
    let a = GF2_128::new(0xffffffffffffffff, 0xffffffffffffffff);
    let expected = GF2_128::new(0xfc10c53d1c96eca8, 0xfe08629e8e4b766a);
    assert_eq!(a.invert().unwrap(), expected);
}

#[test]
fn pow_known_vector_0() {
    let a = GF2_128::from_u64(0xdeadbeefcafe1234);
    assert_eq!(a.pow(0), GF2_128::one());
}

#[test]
fn pow_known_vector_1() {
    let a = GF2_128::from_u64(0xdeadbeefcafe1234);
    assert_eq!(a.pow(1), a);
}

#[test]
fn pow_known_vector_3() {
    let a = GF2_128::from_u64(0xdeadbeefcafe1234);
    let expected = GF2_128::new(0xce2b4a77e44ec116, 0xe555280bb487b18a);
    assert_eq!(a.pow(3), expected);
}

#[test]
fn pow_known_vector_17() {
    let a = GF2_128::from_u64(0xdeadbeefcafe1234);
    let expected = GF2_128::new(0xb8b59dbbe13f7c2a, 0xd4b6b856a3200f2c);
    assert_eq!(a.pow(17), expected);
}

#[test]
fn pow_known_vector_5() {
    let a = GF2_128::new(0x1234567890abcdef, 0x1234567890abcdef);
    let expected = GF2_128::new(0x21d292a413a4f9ef, 0x242c45156eccd574);
    assert_eq!(a.pow(5), expected);
}

#[test]
fn square_known_vector_1() {
    let a = GF2_128::from_u64(0xdeadbeefcafe1234);
    let expected = GF2_128::new(0x5044555401040510, 0x5154445145545455);
    assert_eq!(a.square(), expected);
}

#[test]
fn square_known_vector_2() {
    let a = GF2_128::new(0x1234567890abcdef, 0x1234567890abcdef);
    let expected = GF2_128::new(0x0623bb37c94dd37e, 0x841a9668ec72dfa1);
    assert_eq!(a.square(), expected);
}
