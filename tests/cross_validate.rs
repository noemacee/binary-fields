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
