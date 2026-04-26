use ark_ff::{define_field, Fp2, Fp2Config, SmallFp};

define_field!(
    modulus = "18446744069414584321", // 2^64 - 2^32 + 1
    generator = "7",
    name = Goldilocks,
);

pub struct GoldilocksExt2Config;

impl Fp2Config for GoldilocksExt2Config {
    type Fp = Goldilocks;

    // 7 is the generator of GF(p)* and satisfies 7^((p-1)/2) = -1 (Euler's criterion) → QNR.
    // Montgomery form: 7 * R mod p where R = 2^64, R mod p = 2^32 - 1 = 4294967295
    // 7 * 4294967295 mod p = 30064771065
    const NONRESIDUE: Goldilocks = SmallFp::from_raw(30064771065u64);

    // c1[k] = NONRESIDUE^((p^k - 1) / 2)
    //   k=0: NONRESIDUE^0 = 1       → mont: 1 * (2^32-1) mod p = 4294967295
    //   k=1: NONRESIDUE^((p-1)/2) = -1  → mont: (p-1) * (2^32-1) mod p = 18446744065119617026
    const FROBENIUS_COEFF_FP2_C1: &'static [Goldilocks] = &[
        SmallFp::from_raw(4294967295u64),
        SmallFp::from_raw(18446744065119617026u64),
    ];
}

pub type GoldilocksExt2 = Fp2<GoldilocksExt2Config>;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, One, Zero};

    #[test]
    fn nonresidue_is_qnr() {
        // 7^((p-1)/2) must equal -1 (Euler's criterion)
        let seven = GoldilocksExt2Config::NONRESIDUE;
        let p_minus_1_over_2 = (Goldilocks::characteristic()[0] - 1) / 2;
        assert_eq!(seven.pow([p_minus_1_over_2]), -Goldilocks::one());
    }

    #[test]
    fn ext2_mul_inverse() {
        use ark_std::UniformRand;
        let mut rng = ark_std::test_rng();
        let a = GoldilocksExt2::rand(&mut rng);
        let a_inv = a.inverse().unwrap();
        assert_eq!(a * a_inv, GoldilocksExt2::one());
    }

    #[test]
    fn ext2_sub_self_is_zero() {
        use ark_std::UniformRand;
        let mut rng = ark_std::test_rng();
        let a = GoldilocksExt2::rand(&mut rng);
        assert!((a - a).is_zero());
    }

    #[test]
    fn ext2_frobenius() {
        use ark_std::UniformRand;
        let mut rng = ark_std::test_rng();
        let a = GoldilocksExt2::rand(&mut rng);
        let mut b = a;
        b.frobenius_map_in_place(1);
        // Frobenius twice should be identity (degree-2 extension)
        b.frobenius_map_in_place(1);
        assert_eq!(a, b);
    }
}
