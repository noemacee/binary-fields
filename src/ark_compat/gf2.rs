//! GF(2) = Z/2Z — the base prime field for all binary extension fields.
//!
//! Implemented via `Fp<Gf2Config, 1>` following the same pattern as all
//! arkworks prime fields.  `Gf2Config` implements `FpConfig<1>` with
//! characteristic-2 arithmetic (XOR for addition, AND for multiplication).
//! All operator impls, serialization, `From<integer>`, `PrimeField`,
//! `FftField`, `UniformRand`, etc. are provided for free by `Fp<P, N>`.

use ark_ff::{BigInt, Fp, FpConfig, SqrtPrecomputation};
use core::marker::PhantomData;

/// Config for GF(2).  The modulus is 2; elements are stored as 0 or 1.
pub struct Gf2Config;

impl FpConfig<1> for Gf2Config {
    /// p = 2.
    const MODULUS: BigInt<1> = BigInt([2u64]);

    const ZERO: Fp<Self, 1> = Fp(BigInt([0u64]), PhantomData);
    const ONE:  Fp<Self, 1> = Fp(BigInt([1u64]), PhantomData);

    /// The multiplicative group GF(2)* = {1} is trivial; its only generator is 1.
    const GENERATOR: Fp<Self, 1> = Fp(BigInt([1u64]), PhantomData);

    /// |GF(2)*| = 1 = 2^0 * 1  →  TWO_ADICITY = 0.
    const TWO_ADICITY: u32 = 0;

    /// The unique 2^0 = 1st root of unity is 1.
    const TWO_ADIC_ROOT_OF_UNITY: Fp<Self, 1> = Fp(BigInt([1u64]), PhantomData);

    /// Use Tonelli-Shanks with two_adicity = 0.
    /// For GF(2) the outer loop never executes (both elements are squares),
    /// so `quadratic_nonresidue_to_trace` is never read.
    const SQRT_PRECOMP: Option<SqrtPrecomputation<Fp<Self, 1>>> =
        Some(SqrtPrecomputation::TonelliShanks {
            two_adicity: 0,
            quadratic_nonresidue_to_trace: Fp(BigInt([1u64]), PhantomData),
            trace_of_modulus_minus_one_div_two: &[],
        });

    /// Addition = XOR in characteristic 2.
    fn add_assign(a: &mut Fp<Self, 1>, b: &Fp<Self, 1>) {
        a.0.0[0] ^= b.0.0[0];
    }

    /// Subtraction = addition in characteristic 2.
    fn sub_assign(a: &mut Fp<Self, 1>, b: &Fp<Self, 1>) {
        a.0.0[0] ^= b.0.0[0];
    }

    /// x + x = 0 in characteristic 2.
    fn double_in_place(a: &mut Fp<Self, 1>) {
        a.0.0[0] = 0;
    }

    /// -x = x in characteristic 2.
    fn neg_in_place(_a: &mut Fp<Self, 1>) {}

    /// Multiplication = AND in GF(2).
    fn mul_assign(a: &mut Fp<Self, 1>, b: &Fp<Self, 1>) {
        a.0.0[0] &= b.0.0[0];
    }

    /// x^2 = x in GF(2) (Frobenius is the identity).
    fn square_in_place(_a: &mut Fp<Self, 1>) {}

    /// Only 1 has an inverse; 0 does not.
    fn inverse(a: &Fp<Self, 1>) -> Option<Fp<Self, 1>> {
        if a.0.0[0] == 0 { None } else { Some(*a) }
    }

    /// Accept 0 or 1; reject anything ≥ 2.
    fn from_bigint(r: BigInt<1>) -> Option<Fp<Self, 1>> {
        if r.0[0] < 2 {
            Some(Fp(BigInt([r.0[0]]), PhantomData))
        } else {
            None
        }
    }

    fn into_bigint(a: Fp<Self, 1>) -> BigInt<1> {
        BigInt([a.0.0[0] & 1])
    }

    /// Inner product over GF(2): XOR of pairwise ANDs.
    fn sum_of_products<const T: usize>(
        a: &[Fp<Self, 1>; T],
        b: &[Fp<Self, 1>; T],
    ) -> Fp<Self, 1> {
        let mut acc = 0u64;
        for (x, y) in a.iter().zip(b.iter()) {
            acc ^= x.0.0[0] & y.0.0[0];
        }
        Fp(BigInt([acc & 1]), PhantomData)
    }
}

/// GF(2) as an arkworks `PrimeField`.
pub type Gf2 = Fp<Gf2Config, 1>;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{AdditiveGroup, BigInt, Field, One, PrimeField, Zero};

    fn g0() -> Gf2 { Gf2::zero() }
    fn g1() -> Gf2 { Gf2::one() }

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
    fn from_int_mod2() {
        assert_eq!(Gf2::from(0u8), g0());
        assert_eq!(Gf2::from(1u8), g1());
        assert_eq!(Gf2::from(2u8), g0());
        assert_eq!(Gf2::from(3u8), g1());
        assert_eq!(Gf2::from(-1i64), g1());
        assert_eq!(Gf2::from(-2i64), g0());
    }

    #[test]
    fn double_is_zero() {
        let mut x = g1();
        x.double_in_place();
        assert_eq!(x, g0());
    }

    #[test]
    fn serialization_roundtrip() {
        use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
        for val in [g0(), g1()] {
            let mut bytes = ark_std::vec::Vec::new();
            val.serialize_uncompressed(&mut bytes).unwrap();
            assert_eq!(bytes.len(), 1);
            let back = Gf2::deserialize_uncompressed(bytes.as_slice()).unwrap();
            assert_eq!(back, val);
        }
    }

    #[test]
    fn prime_field_roundtrip() {
        let bi: BigInt<1> = g1().into_bigint();
        assert_eq!(bi.0[0], 1);
        assert_eq!(Gf2::from_bigint(bi), Some(g1()));
        assert_eq!(Gf2::from_bigint(BigInt([2u64])), None);
    }
}
