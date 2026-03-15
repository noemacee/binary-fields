//! `Gf128Config` — the concrete config for GF(2^128) with irreducible polynomial
//! f(z) = z^128 + z^7 + z^2 + z + 1  (the GCM polynomial, also NIST GF(2^128)).
//!
//! All arithmetic delegates to the existing algorithms in `crate::arithmetic`,
//! `crate::reduce`, and `crate::invert`.

use crate::ark_compat::{BinaryField, BinaryFieldConfig};

/// Config for GF(2^128).
pub struct Gf128Config;

impl BinaryFieldConfig<2> for Gf128Config {
    // DEGREE defaults to N * 64 = 128 — no override needed.

    const ZERO: BinaryField<Self, 2> = BinaryField([0u64; 2], core::marker::PhantomData);

    const ONE: BinaryField<Self, 2> = BinaryField([1u64, 0u64], core::marker::PhantomData);

    // -1 = 1 in characteristic 2
    const NEG_ONE: BinaryField<Self, 2> = BinaryField([1u64, 0u64], core::marker::PhantomData);

    // A standard primitive element for GF(2^128) with this polynomial.
    // z (the polynomial variable) has multiplicative order 2^128 - 1.
    const GENERATOR: BinaryField<Self, 2> = BinaryField([2u64, 0u64], core::marker::PhantomData);

    /// Carry-less multiplication modulo f(z) = z^128 + z^7 + z^2 + z + 1.
    /// Uses Algorithm 2.34 (comb method) + Algorithm 2.41 (fast reduction).
    fn mul_assign(a: &mut [u64; 2], b: &[u64; 2]) {
        use crate::GF2_128;
        let av = GF2_128::new(a[0], a[1]);
        let bv = GF2_128::new(b[0], b[1]);
        let r = av.mul_2_34(bv);
        a[0] = r.0[0];
        a[1] = r.0[1];
    }

    /// Squaring via bit-expansion table (Algorithm 2.39) + fast reduction (Algorithm 2.41).
    fn square_in_place(a: &mut [u64; 2]) {
        use crate::GF2_128;
        let av = GF2_128::new(a[0], a[1]);
        let r = av.square_2_39();
        a[0] = r.0[0];
        a[1] = r.0[1];
    }

    /// Inversion via binary GCD algorithm (Algorithm 2.49).
    fn inverse(a: &[u64; 2]) -> Option<[u64; 2]> {
        use crate::GF2_128;
        let av = GF2_128::new(a[0], a[1]);
        av.invert_2_49().map(|r| [r.0[0], r.0[1]])
    }
}

/// GF(2^128) as an arkworks `Field`.
pub type Gf128 = BinaryField<Gf128Config, 2>;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{AdditiveGroup, Field, Zero, One};
    use ark_std::UniformRand;

    // Convenience: construct a GF(2^128) element from low/high 64-bit limbs.
    fn gf(lo: u64, hi: u64) -> Gf128 { Gf128::new([lo, hi]) }

    #[test]
    fn zero_is_zero() {
        assert!(Gf128::zero().is_zero());
    }

    #[test]
    fn one_is_one() {
        assert!(Gf128::one().is_one());
    }

    // From<integer> is the ring homomorphism Z → GF(2^128): n maps to n mod 2.
    #[test]
    fn from_integer_is_mod2() {
        assert!(Gf128::from(0u64).is_zero());
        assert!(Gf128::from(1u64).is_one());
        assert!(Gf128::from(2u64).is_zero());
        assert!(Gf128::from(3u64).is_one());
        assert!(Gf128::from(0xdeadbeefcafe1234u64).is_zero()); // even
        assert!(Gf128::from(0xdeadbeefcafe1235u64).is_one());  // odd
        assert!(Gf128::from(-1i64).is_one());
        assert!(Gf128::from(-2i64).is_zero());
    }

    #[test]
    fn add_self_is_zero() {
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        assert!((a + a).is_zero());
    }

    #[test]
    fn mul_by_one_is_identity() {
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        assert_eq!(a * Gf128::one(), a);
    }

    #[test]
    fn mul_by_zero_is_zero() {
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        assert!((a * Gf128::zero()).is_zero());
    }

    #[test]
    fn inverse_roundtrip() {
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let a_inv = a.inverse().unwrap();
        assert_eq!(a * a_inv, Gf128::one());
    }

    #[test]
    fn inverse_of_zero_is_none() {
        assert!(Gf128::zero().inverse().is_none());
    }

    #[test]
    fn square_consistent_with_mul() {
        let a = gf(0x1234abcd5678ef01, 0);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn square_consistent_with_mul_both_limbs() {
        let a = gf(0x1234abcd5678ef01, 0xfedcba9876543210);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn frobenius_is_squaring() {
        let a = gf(0xdeadbeef_cafebabe, 0x0102030405060708);
        let mut b = a;
        b.frobenius_map_in_place(1);
        assert_eq!(b, a.square());
    }

    #[test]
    fn frobenius_power_3() {
        let a = gf(0xdeadbeef_cafebabe, 0x0102030405060708);
        let mut b = a;
        b.frobenius_map_in_place(3); // x^(2^3) = x^8
        let expected = a.pow([8u64]);
        assert_eq!(b, expected);
    }

    #[test]
    fn sqrt_roundtrip() {
        let a = gf(0xabcd_1234_5678_ef01, 0);
        let s = a.sqrt().unwrap();
        assert_eq!(s.square(), a);
    }

    #[test]
    fn from_gf2_element() {
        use crate::ark_compat::Gf2;
        let zero = Gf128::from_base_prime_field(Gf2::zero());
        let one  = Gf128::from_base_prime_field(Gf2::one());
        assert!(zero.is_zero());
        assert!(one.is_one());
    }

    #[test]
    fn to_base_prime_field_roundtrip() {
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let bits: ark_std::vec::Vec<_> = a.to_base_prime_field_elements().collect();
        assert_eq!(bits.len(), 128);
        let b = Gf128::from_base_prime_field_elems(bits).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn serialization_roundtrip() {
        use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let mut bytes = ark_std::vec::Vec::new();
        a.serialize_uncompressed(&mut bytes).unwrap();
        assert_eq!(bytes.len(), 16); // 128 bits = 16 bytes
        let b = Gf128::deserialize_uncompressed(bytes.as_slice()).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn mul_z127_times_z_is_reduction() {
        // z^127 * z = z^128 = z^7 + z^2 + z + 1 = 0x87
        let z127 = gf(0, 1u64 << 63);
        let z1   = gf(0b10, 0);
        let expected = gf(0x87, 0);
        assert_eq!(z127 * z1, expected);
    }

    #[test]
    fn characteristic_is_2() {
        let one = Gf128::one();
        assert!((one + one).is_zero());
        assert_eq!(Gf128::characteristic(), &[2u64]);
    }

    #[test]
    fn extension_degree_is_128() {
        assert_eq!(Gf128::extension_degree(), 128);
    }

    #[test]
    fn uniform_rand_not_always_zero() {
        let mut rng = ark_std::test_rng();
        let samples: ark_std::vec::Vec<Gf128> = (0..10).map(|_| Gf128::rand(&mut rng)).collect();
        assert!(samples.iter().any(|x| !x.is_zero()));
    }

    #[test]
    fn pow_is_consistent() {
        let a = gf(0x1234567890abcdef, 0xfedcba0987654321);
        // a^4 = (a^2)^2
        assert_eq!(a.pow([4u64]), a.square().square());
        // a^0 = 1
        assert!(a.pow([0u64]).is_one());
        // a^1 = a
        assert_eq!(a.pow([1u64]), a);
    }

    #[test]
    fn mul_associative() {
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let b = gf(0x0102030405060708, 0xfedcba9876543210);
        let c = gf(0x1111222233334444, 0x5555666677778888);
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn mul_by_base_prime_field_zero_gives_zero() {
        use crate::ark_compat::Gf2;
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        assert!(a.mul_by_base_prime_field(&Gf2::zero()).is_zero());
    }

    #[test]
    fn mul_by_base_prime_field_one_gives_self() {
        use crate::ark_compat::Gf2;
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        assert_eq!(a.mul_by_base_prime_field(&Gf2::one()), a);
    }

    #[test]
    fn legendre_zero_is_zero_symbol() {
        use ark_ff::LegendreSymbol;
        assert_eq!(Gf128::zero().legendre(), LegendreSymbol::Zero);
    }

    #[test]
    fn legendre_nonzero_is_quadratic_residue() {
        use ark_ff::LegendreSymbol;
        // In GF(2^m), every nonzero element is a square (Frobenius is bijective).
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        assert_eq!(a.legendre(), LegendreSymbol::QuadraticResidue);
        assert_eq!(Gf128::one().legendre(), LegendreSymbol::QuadraticResidue);
    }

    #[test]
    fn double_in_place_is_zero() {
        let mut a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        a.double_in_place();
        assert!(a.is_zero());
    }

    #[test]
    fn neg_in_place_is_identity() {
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let mut b = a;
        b.neg_in_place();
        assert_eq!(b, a);
    }

    #[test]
    fn inverse_in_place_roundtrip() {
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let mut b = a;
        b.inverse_in_place().unwrap();
        assert_eq!(a * b, Gf128::one());
    }

    #[test]
    fn inverse_in_place_zero_returns_none() {
        let mut z = Gf128::zero();
        assert!(z.inverse_in_place().is_none());
    }

    #[test]
    fn from_base_prime_field_elems_wrong_count_is_none() {
        use crate::ark_compat::Gf2;
        // Too few
        let bits: ark_std::vec::Vec<Gf2> = (0..127).map(|_| Gf2::zero()).collect();
        assert!(Gf128::from_base_prime_field_elems(bits).is_none());
        // Too many
        let bits: ark_std::vec::Vec<Gf2> = (0..129).map(|_| Gf2::zero()).collect();
        assert!(Gf128::from_base_prime_field_elems(bits).is_none());
    }

    #[test]
    fn sum_of_elements() {
        let a = gf(0x1234, 0);
        let b = gf(0x5678, 0);
        let c = gf(0xabcd, 0);
        let expected = a + b + c;
        let got: Gf128 = [a, b, c].iter().copied().sum();
        assert_eq!(got, expected);
    }

    #[test]
    fn product_of_elements() {
        let a = gf(0x1234, 0);
        let b = gf(0x5678, 0);
        let c = gf(0xabcd, 0);
        let expected = a * b * c;
        let got: Gf128 = [a, b, c].iter().copied().product();
        assert_eq!(got, expected);
    }

    #[test]
    fn from_random_bytes_roundtrip() {
        let input: [u8; 16] = [
            0xde, 0xad, 0xbe, 0xef, 0xca, 0xfe, 0x12, 0x34,
            0xab, 0xad, 0x1d, 0xea, 0x12, 0x34, 0x56, 0x78,
        ];
        use ark_ff::Field;
        let (elem, _flags) = Gf128::from_random_bytes_with_flags::<ark_serialize::EmptyFlags>(&input).unwrap();
        // Serialise back and compare bytes.
        use ark_serialize::CanonicalSerialize;
        let mut bytes = ark_std::vec::Vec::new();
        elem.serialize_uncompressed(&mut bytes).unwrap();
        assert_eq!(&bytes, &input);
    }
}
