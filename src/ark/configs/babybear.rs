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

    // SmallFp::new converts actual field values to Montgomery form internally.
    fn fp(x: u32) -> BabyBear { SmallFp::new(x) }
    fn ext2(c0: u32, c1: u32) -> BabyBearExt2 { BabyBearExt2::new(fp(c0), fp(c1)) }
    fn ext4(a0: u32, a1: u32, b0: u32, b1: u32) -> BabyBearExt4 {
        BabyBearExt4::new(ext2(a0, a1), ext2(b0, b1))
    }

    // Reference values computed independently via Python galois library (actual field values).

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
    fn ext2_known_vectors() {
        // Python galois: GF(2013265921^2, x^2 - 11), element = c0 + c1*x
        let cases: &[((u32,u32),(u32,u32),(u32,u32),(u32,u32),(u32,u32))] = &[
            // (a, b, a*b, a^2, a^-1)
            ((1804778276,882461865),(1823380466,1303846403),(888764645,283232932),(1651346491,1502787040),(1918449551,93233662)),
            ((652660206,1977080709),(872477412,1357532441),(1254130409,343743121),(330801865,360943128),(752219132,883467310)),
            ((1795463164,694750683),(211159283,585134930),(195343840,306542283),(769337788,1277640436),(152414028,1365467642)),
        ];
        for &(a,b,ab,a2,ainv) in cases {
            let fa = ext2(a.0, a.1);
            let fb = ext2(b.0, b.1);
            assert_eq!(fa * fb, ext2(ab.0, ab.1), "mul mismatch");
            assert_eq!(fa.square(), ext2(a2.0, a2.1), "square mismatch");
            assert_eq!(fa.inverse().unwrap(), ext2(ainv.0, ainv.1), "inverse mismatch");
        }
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
        for _ in 0..4 {
            b.frobenius_map_in_place(1);
        }
        assert_eq!(a, b);
    }

    #[test]
    fn ext4_known_vectors() {
        // Python galois: GF(p^4, y^4-11), tower repr (fp2_c0.c0, fp2_c0.c1, fp2_c1.c0, fp2_c1.c1)
        // Basis: element = (c0+c2*x) + (c1+c3*x)*y  where x^2=11, y^2=x
        let cases: &[((u32,u32,u32,u32),(u32,u32,u32,u32),(u32,u32,u32,u32),(u32,u32,u32,u32),(u32,u32,u32,u32))] = &[
            // (a, b, a*b, a^2, a^-1)
            (
                (1351175553,608336906,1652885869,1257676251),
                (563203088,1494985531,196806800,1177587890),
                (1515297101,1580164789,674001329,1505776063),
                (303661535,1539484767,772876596,326653598),
                (1935570845,1256699769,59571902,980104888),
            ),
            (
                (1641869219,231830436,1872526691,1379046964),
                (1857382466,241661782,854703159,1736245966),
                (1368841940,392004790,904677987,1721034113),
                (578001524,264368266,279074372,1565007363),
                (95934421,1932578482,223618845,1229294099),
            ),
            (
                (538690774,1777457902,66736977,1968253536),
                (568509383,183667247,155776591,1358682039),
                (842834796,1534896902,1229077902,1042495550),
                (535042214,940428058,1452624508,1256366215),
                (872111318,1750449300,297932768,1255367450),
            ),
        ];
        for &(a,b,ab,a2,ainv) in cases {
            let fa = ext4(a.0,a.1,a.2,a.3);
            let fb = ext4(b.0,b.1,b.2,b.3);
            assert_eq!(fa * fb, ext4(ab.0,ab.1,ab.2,ab.3), "mul mismatch");
            assert_eq!(fa.square(), ext4(a2.0,a2.1,a2.2,a2.3), "square mismatch");
            assert_eq!(fa.inverse().unwrap(), ext4(ainv.0,ainv.1,ainv.2,ainv.3), "inverse mismatch");
        }
    }
}
