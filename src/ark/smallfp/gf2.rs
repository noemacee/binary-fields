use ark_ff::{BigInt, LegendreSymbol, SmallFp, SmallFpConfig, SqrtPrecomputation};

pub struct Gf2SmallFpConfig;

impl SmallFpConfig for Gf2SmallFpConfig {
    type T = u8;

    const MODULUS: u8 = 2;
    const MODULUS_U128: u128 = 2;

    const GENERATOR: SmallFp<Self> = SmallFp::from_raw(1);
    const ZERO: SmallFp<Self> = SmallFp::from_raw(0);
    const ONE: SmallFp<Self> = SmallFp::from_raw(1);
    const NEG_ONE: SmallFp<Self> = SmallFp::from_raw(1); // -1 = 1 in char 2

    const TWO_ADICITY: u32 = 0; // |GF(2)*| = 1 = 2^0 * 1
    const TWO_ADIC_ROOT_OF_UNITY: SmallFp<Self> = SmallFp::from_raw(1);

    const SQRT_PRECOMP: Option<SqrtPrecomputation<SmallFp<Self>>> =
        Some(SqrtPrecomputation::CharTwo); // sqrt(x) = x in GF(2)

    #[inline(always)]
    fn add_assign(a: &mut SmallFp<Self>, b: &SmallFp<Self>) {
        a.value ^= b.value;
    }

    #[inline(always)]
    fn sub_assign(a: &mut SmallFp<Self>, b: &SmallFp<Self>) {
        a.value ^= b.value;
    }

    #[inline(always)]
    fn double_in_place(a: &mut SmallFp<Self>) {
        a.value = 0;
    }

    #[inline(always)]
    fn neg_in_place(_a: &mut SmallFp<Self>) {}

    #[inline(always)]
    fn mul_assign(a: &mut SmallFp<Self>, b: &SmallFp<Self>) {
        a.value &= b.value;
    }

    #[inline(always)]
    fn square_in_place(_a: &mut SmallFp<Self>) {}

    fn sum_of_products<const T: usize>(
        a: &[SmallFp<Self>; T],
        b: &[SmallFp<Self>; T],
    ) -> SmallFp<Self> {
        let mut acc = 0u8;
        for (x, y) in a.iter().zip(b.iter()) {
            acc ^= x.value & y.value;
        }
        SmallFp::from_raw(acc & 1)
    }

    #[inline(always)]
    fn inverse(a: &SmallFp<Self>) -> Option<SmallFp<Self>> {
        if a.value == 0 { None } else { Some(*a) }
    }

    #[inline(always)]
    fn new(value: u8) -> SmallFp<Self> {
        SmallFp::from_raw(value % 2)
    }

    fn from_bigint(r: BigInt<1>) -> Option<SmallFp<Self>> {
        if r.0[0] < 2 { Some(SmallFp::from_raw(r.0[0] as u8)) } else { None }
    }

    fn into_bigint(a: SmallFp<Self>) -> BigInt<1> {
        BigInt([a.value as u64])
    }

    fn legendre(_: &SmallFp<Self>) -> LegendreSymbol {
        unimplemented!("Legendre symbol is undefined in characteristic 2")
    }
}

pub type Gf2SmallFp = SmallFp<Gf2SmallFpConfig>;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{AdditiveGroup, Field, One, PrimeField, Zero};
    use ark_std::UniformRand;

    fn g0() -> Gf2SmallFp { Gf2SmallFp::zero() }
    fn g1() -> Gf2SmallFp { Gf2SmallFp::one() }

    #[test]
    fn add_is_xor() {
        assert_eq!(g0() + g0(), g0());
        assert_eq!(g0() + g1(), g1());
        assert_eq!(g1() + g0(), g1());
        assert_eq!(g1() + g1(), g0());
    }

    #[test]
    fn mul_is_and() {
        assert_eq!(g0() * g0(), g0());
        assert_eq!(g0() * g1(), g0());
        assert_eq!(g1() * g1(), g1());
    }

    #[test]
    fn neg_is_identity() {
        assert_eq!(-g0(), g0());
        assert_eq!(-g1(), g1());
    }

    #[test]
    fn square_is_identity() {
        assert_eq!(g0().square(), g0());
        assert_eq!(g1().square(), g1());
    }

    #[test]
    fn inverse() {
        assert!(g0().inverse().is_none());
        assert_eq!(g1().inverse(), Some(g1()));
    }

    #[test]
    fn double_is_zero() {
        assert_eq!(g0().double(), g0());
        assert_eq!(g1().double(), g0());
    }

    #[test]
    fn from_int_mod2() {
        assert_eq!(Gf2SmallFp::from(0u8), g0());
        assert_eq!(Gf2SmallFp::from(1u8), g1());
        assert_eq!(Gf2SmallFp::from(2u8), g0());
        assert_eq!(Gf2SmallFp::from(3u8), g1());
    }

    #[test]
    fn characteristic_is_2() {
        assert_eq!(Gf2SmallFp::characteristic(), &[2u64]);
    }

    #[test]
    fn prime_field_constants() {
        assert_eq!(Gf2SmallFp::MODULUS, ark_ff::BigInt([2u64]));
        assert_eq!(Gf2SmallFp::MODULUS_MINUS_ONE_DIV_TWO, ark_ff::BigInt([0u64]));
        assert_eq!(Gf2SmallFp::TRACE, ark_ff::BigInt([1u64]));
        assert_eq!(Gf2SmallFp::TRACE_MINUS_ONE_DIV_TWO, ark_ff::BigInt([0u64]));
    }

    #[test]
    fn uniform_rand_samples_zero() {
        let mut rng = ark_std::test_rng();
        let samples: Vec<_> = (0..1000).map(|_| Gf2SmallFp::rand(&mut rng)).collect();
        assert!(samples.iter().any(|x| x.is_zero()), "0 should be reachable");
        assert!(samples.iter().any(|x| x.is_one()), "1 should be reachable");
    }

    #[test]
    fn serialization_roundtrip() {
        use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
        for val in [g0(), g1()] {
            let mut bytes = ark_std::vec::Vec::new();
            val.serialize_uncompressed(&mut bytes).unwrap();
            assert_eq!(bytes.len(), 1);
            let back = Gf2SmallFp::deserialize_uncompressed(bytes.as_slice()).unwrap();
            assert_eq!(back, val);
        }
    }
}
