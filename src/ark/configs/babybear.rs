use ark_ff::{define_field, Fp2, Fp2Config, Fp4, Fp4Config, SmallFp};

define_field!(
    modulus = "2013265921", // 2^31 - 2^27 + 1
    generator = "31",
    name = BabyBear,
);

// ── Fp2 = BabyBear[x] / (x^2 - 11) ──────────────────────────────────────────
// 11 is a QNR in BabyBear: 11^((p-1)/2) = -1 mod p (Euler's criterion).
// Montgomery constants use R = 2^31, R mod p = 134217727.

pub struct BabyBearExt2Config;

impl Fp2Config for BabyBearExt2Config {
    type Fp = BabyBear;

    // 11 in Montgomery form: 11 * 134217727 mod p = 1476394997
    const NONRESIDUE: BabyBear = SmallFp::from_raw(1476394997u32);

    // c1[i] = 11^((p^i - 1) / 2) mod p, in Montgomery form
    //   i=0: 1          → 1 * R mod p = 134217727
    //   i=1: p-1 = -1   → (p-1) * R mod p = 1879048194
    const FROBENIUS_COEFF_FP2_C1: &'static [BabyBear] = &[
        SmallFp::from_raw(134217727u32),
        SmallFp::from_raw(1879048194u32),
    ];
}

pub type BabyBearExt2 = Fp2<BabyBearExt2Config>;

// ── Fp4 = BabyBearExt2[y] / (y^2 - (0,1)) ───────────────────────────────────
// Tower: Fp4 is a degree-2 extension of Fp2, giving 4*31 = 124-bit field.
// NONRESIDUE = (0, 1) in Fp2, i.e. the element x (root of x^2 - 11 = 0).
// This satisfies the Fp4Config constraint that NONRESIDUE must equal (0, 1).

pub struct BabyBearExt4Config;

impl Fp4Config for BabyBearExt4Config {
    type Fp2Config = BabyBearExt2Config;

    // NONRESIDUE = (0, 1) in Fp2 — required by the Fp4Config tower construction.
    const NONRESIDUE: BabyBearExt2 = BabyBearExt2::new(
        SmallFp::from_raw(0u32),
        SmallFp::from_raw(134217727u32), // 1 in Montgomery form
    );

    // c1[i] = 11^((p^i - 1) / 4) mod p, in Montgomery form
    //   i=0: 1           → from_raw(134217727)
    //   i=1: 1728404513  → from_raw(1243376265)
    //   i=2: p-1 = -1    → from_raw(1879048194)
    //   i=3: 284861408   → from_raw(769889656)
    const FROBENIUS_COEFF_FP4_C1: &'static [BabyBear] = &[
        SmallFp::from_raw(134217727u32),
        SmallFp::from_raw(1243376265u32),
        SmallFp::from_raw(1879048194u32),
        SmallFp::from_raw(769889656u32),
    ];
}

pub type BabyBearExt4 = Fp4<BabyBearExt4Config>;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, One, Zero};
    use ark_std::UniformRand;

    #[test]
    fn nonresidue_is_qnr() {
        let eleven = BabyBearExt2Config::NONRESIDUE;
        let p_minus_1_over_2 = (BabyBear::characteristic()[0] - 1) / 2;
        assert_eq!(eleven.pow([p_minus_1_over_2]), -BabyBear::one());
    }

    #[test]
    fn ext2_mul_inverse() {
        let mut rng = ark_std::test_rng();
        let a = BabyBearExt2::rand(&mut rng);
        assert_eq!(a * a.inverse().unwrap(), BabyBearExt2::one());
    }

    #[test]
    fn ext2_sub_self_is_zero() {
        let mut rng = ark_std::test_rng();
        let a = BabyBearExt2::rand(&mut rng);
        assert!((a - a).is_zero());
    }

    #[test]
    fn ext2_frobenius() {
        let mut rng = ark_std::test_rng();
        let a = BabyBearExt2::rand(&mut rng);
        let mut b = a;
        b.frobenius_map_in_place(1);
        b.frobenius_map_in_place(1);
        assert_eq!(a, b);
    }

    #[test]
    fn ext4_mul_inverse() {
        let mut rng = ark_std::test_rng();
        let a = BabyBearExt4::rand(&mut rng);
        assert_eq!(a * a.inverse().unwrap(), BabyBearExt4::one());
    }

    #[test]
    fn ext4_sub_self_is_zero() {
        let mut rng = ark_std::test_rng();
        let a = BabyBearExt4::rand(&mut rng);
        assert!((a - a).is_zero());
    }

    #[test]
    fn ext4_frobenius() {
        let mut rng = ark_std::test_rng();
        let a = BabyBearExt4::rand(&mut rng);
        let mut b = a;
        // Frobenius 4 times = identity on degree-4 extension
        for _ in 0..4 {
            b.frobenius_map_in_place(1);
        }
        assert_eq!(a, b);
    }
}
