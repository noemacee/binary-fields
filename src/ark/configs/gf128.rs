//! Config for GF(2^128) with the GCM polynomial — optimized algorithms.

use crate::ark::{BinaryField, BinaryFieldConfig};
use crate::polynomials::{self, IrreduciblePoly};

/// Config for GF(2^128) with optimized multiplication and squaring.
///
/// Uses Algorithm 2.34 (comb) for multiplication and Algorithm 2.39
/// (bit-expansion table) for squaring, both with Algorithm 2.41 fast
/// reduction. Inverse falls back to the generic Algorithm 2.48.
pub struct Gf128Config;

impl Gf128Config {
    /// The irreducible polynomial that defines this field.
    pub const POLY: IrreduciblePoly = polynomials::GF128_GCM;
}

impl BinaryFieldConfig<2> for Gf128Config {
    // f(z) = z^128 + z^7 + z^2 + z + 1  (GCM polynomial, NIST SP 800-38D)
    const ALPHA_POW_M: [u64; 2] = [
        polynomials::GF128_GCM.alpha_pow_m[0],
        polynomials::GF128_GCM.alpha_pow_m[1],
    ];

    const ZERO: BinaryField<Self, 2> = BinaryField([0u64; 2], core::marker::PhantomData);
    const ONE: BinaryField<Self, 2> = BinaryField([1u64, 0u64], core::marker::PhantomData);

    fn mul_assign(a: &mut [u64; 2], b: &[u64; 2]) {
        *a = crate::fields::z128_z7_z2_z1::mul_2_34(*a, *b);
    }

    fn square_in_place(a: &mut [u64; 2]) {
        *a = crate::fields::z128_z7_z2_z1::square_2_39(*a);
    }
}

/// GF(2^128) as an arkworks `Field`, with optimized arithmetic.
pub type Gf128 = BinaryField<Gf128Config, 2>;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{AdditiveGroup, Field, One, Zero};
    use ark_std::UniformRand;

    fn gf(lo: u64, hi: u64) -> Gf128 {
        Gf128::new([lo, hi])
    }

    #[test]
    fn zero_is_zero() {
        assert!(Gf128::zero().is_zero());
    }

    #[test]
    fn one_is_one() {
        assert!(Gf128::one().is_one());
    }

    #[test]
    fn from_integer_is_mod2() {
        assert!(Gf128::from(0u64).is_zero());
        assert!(Gf128::from(1u64).is_one());
        assert!(Gf128::from(2u64).is_zero());
        assert!(Gf128::from(3u64).is_one());
        assert!(Gf128::from(0xdeadbeefcafe1234u64).is_zero()); // even
        assert!(Gf128::from(0xdeadbeefcafe1235u64).is_one()); // odd
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
        b.frobenius_map_in_place(3);
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
        use crate::ark::Gf2;
        let zero = Gf128::from_base_prime_field(Gf2::zero());
        let one = Gf128::from_base_prime_field(Gf2::one());
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
        use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let mut bytes = ark_std::vec::Vec::new();
        a.serialize_uncompressed(&mut bytes).unwrap();
        assert_eq!(bytes.len(), 16);
        let b = Gf128::deserialize_uncompressed(bytes.as_slice()).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn mul_z127_times_z_is_reduction() {
        let z127 = gf(0, 1u64 << 63);
        let z1 = gf(0b10, 0);
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
        assert_eq!(a.pow([4u64]), a.square().square());
        assert!(a.pow([0u64]).is_one());
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
        use crate::ark::Gf2;
        let a = gf(0xdeadbeefcafe1234, 0xabad1dea12345678);
        assert!(a.mul_by_base_prime_field(&Gf2::zero()).is_zero());
    }

    #[test]
    fn mul_by_base_prime_field_one_gives_self() {
        use crate::ark::Gf2;
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
        use crate::ark::Gf2;
        let bits: ark_std::vec::Vec<Gf2> = (0..127).map(|_| Gf2::zero()).collect();
        assert!(Gf128::from_base_prime_field_elems(bits).is_none());
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
    fn poly_constant_matches_config() {
        assert_eq!(Gf128Config::POLY.degree, 128);
        assert_eq!(Gf128Config::POLY.terms, &[7, 2, 1, 0]);
        assert_eq!(Gf128Config::POLY.alpha_pow_m, &[0x87u64, 0]);
    }

    // ── Known-vector tests (Python `galois` library) ─────────────────────────

    #[test]
    fn mul_known_vectors() {
        let cases: &[([u64; 2], [u64; 2], [u64; 2])] = &[
            (
                [0x6d2e9cb704e5f153, 0xe65ea5827839d4f0],
                [0xee408ef1996ab63a, 0xe4dd43466d146df7],
                [0x455a77734b91cf8e, 0x670ac21e287ccb0b],
            ),
            (
                [0x60ac95819cb87ff0, 0xc5bb01ffba43a322],
                [0xfd023762a0c39500, 0x77109b38f72041f0],
                [0xccd1226c225306a5, 0xbe424c494361e678],
            ),
            (
                [0x65138aaf54ded307, 0xb1dfd59a5d604585],
                [0x686c514eab1f1cba, 0xec6fd3b6b9101400],
                [0x6ca978fa2b10ab4e, 0x5a6e88cc1f8b54b2],
            ),
            (
                [0xd8790dfed23e6153, 0xb68403d3619ed109],
                [0xa9309903820d4638, 0x46b94a5d4486fe1b],
                [0x0e5856d8dca261e9, 0x62aded75cc9cb38e],
            ),
            (
                [0x404b5b85fe11370d, 0x4043e4a9f9072dd5],
                [0x852e48672d75dadc, 0x63fd0c3450ae1132],
                [0x499b3c175bd18a0b, 0xb9033cc9a69c6fb4],
            ),
            (
                [0x7cd2999122d31814, 0x513ce8596968b954],
                [0xef91d8d492cf8b1c, 0x792918de94ca7438],
                [0x227fdc899f525e12, 0x3b258f950e9f18ab],
            ),
        ];
        for &(a, b, exp) in cases {
            assert_eq!(gf(a[0], a[1]) * gf(b[0], b[1]), gf(exp[0], exp[1]));
        }
    }

    #[test]
    fn square_known_vectors() {
        let cases: &[([u64; 2], [u64; 2])] = &[
            (
                [0xf328bc0e21a4f22b, 0xed80b1e4fffd2abf],
                [0x040146256afb613f, 0xd013c4631efca40f],
            ),
            (
                [0x85ec45f9adc68838, 0x27c1246601f2a83a],
                [0x44fd792bbd82bc92, 0x56d3e4d5046b332d],
            ),
            (
                [0xa969180d6454ce16, 0x8c5e7590b0289ba9],
                [0x4f122cf1342cec34, 0xadf9c8e7e2964772],
            ),
            (
                [0x9a9d0b9a2d62473e, 0x4c72dfcbf051186c],
                [0x2f59e383b6df44dc, 0x18fea8641c47d234],
            ),
            (
                [0x5bc843e7294f608f, 0x7ee729d6b5e27424],
                [0x579dbe43f7725113, 0xd1c2f6a92d6b695a],
            ),
            (
                [0x739c0b24557cea39, 0xc9ae1a7f63380165],
                [0x7f88af9054c975e9, 0x8461b6fca493c431],
            ),
        ];
        for &(a, exp) in cases {
            assert_eq!(gf(a[0], a[1]).square(), gf(exp[0], exp[1]));
        }
    }

    #[test]
    fn inverse_known_vectors() {
        let cases: &[([u64; 2], [u64; 2])] = &[
            (
                [0xf9f88e42711b4eb1, 0x425306909080d43e],
                [0x56901c591c2d1fd1, 0x01f71baa4c8b6ca1],
            ),
            (
                [0x8a6c15d3d2b5d226, 0x93dfca97387ba07a],
                [0x653505fc7b509e39, 0x9e23f9f61b5f5615],
            ),
            (
                [0x758dce0815ab5499, 0x557be408ecf4842a],
                [0x53f0552aa5605461, 0x3052fd630fb353da],
            ),
            (
                [0xe94435dc7fccbbe8, 0x02ba68a0310fbbd5],
                [0xaad21a58c89e4d5b, 0x652db4eca00e266b],
            ),
            (
                [0x273804ca0cdfe470, 0xceb5e77d473a6c8a],
                [0xe34835721f20df0d, 0xff244027e1b23829],
            ),
            (
                [0x88b065c1977be8cb, 0xbe7a72ffb139c4f1],
                [0x4f12e624e1a0fe4c, 0xa1eaf28f16e329f4],
            ),
        ];
        for &(a, exp) in cases {
            assert_eq!(gf(a[0], a[1]).inverse(), Some(gf(exp[0], exp[1])));
        }
    }

    #[test]
    fn from_random_bytes_roundtrip() {
        let input: [u8; 16] = [
            0xde, 0xad, 0xbe, 0xef, 0xca, 0xfe, 0x12, 0x34, 0xab, 0xad, 0x1d, 0xea, 0x12, 0x34,
            0x56, 0x78,
        ];
        use ark_ff::Field;
        let (elem, _flags) =
            Gf128::from_random_bytes_with_flags::<ark_serialize::EmptyFlags>(&input).unwrap();
        use ark_serialize::CanonicalSerialize;
        let mut bytes = ark_std::vec::Vec::new();
        elem.serialize_uncompressed(&mut bytes).unwrap();
        assert_eq!(&bytes, &input);
    }
}
