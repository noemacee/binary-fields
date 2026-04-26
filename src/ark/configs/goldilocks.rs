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
    use ark_std::UniformRand;

    // SmallFp::new converts actual field values to Montgomery form internally.
    fn fp(x: u64) -> Goldilocks { SmallFp::new(x) }
    fn ext2(c0: u64, c1: u64) -> GoldilocksExt2 { GoldilocksExt2::new(fp(c0), fp(c1)) }

    #[test]
    fn nonresidue_is_qnr() {
        let seven = GoldilocksExt2Config::NONRESIDUE;
        let p_minus_1_over_2 = (Goldilocks::characteristic()[0] - 1) / 2;
        assert_eq!(seven.pow([p_minus_1_over_2]), -Goldilocks::one());
    }

    #[test]
    fn ext2_mul_inverse() {
        let mut rng = ark_std::test_rng();
        let a = GoldilocksExt2::rand(&mut rng);
        assert_eq!(a * a.inverse().unwrap(), GoldilocksExt2::one());
    }

    #[test]
    fn ext2_sub_self_is_zero() {
        let mut rng = ark_std::test_rng();
        let a = GoldilocksExt2::rand(&mut rng);
        assert!((a - a).is_zero());
    }

    #[test]
    fn ext2_frobenius() {
        let mut rng = ark_std::test_rng();
        let a = GoldilocksExt2::rand(&mut rng);
        let mut b = a;
        b.frobenius_map_in_place(1);
        b.frobenius_map_in_place(1);
        assert_eq!(a, b);
    }

    #[test]
    fn ext2_known_vectors() {
        // Python galois: GF((2^64-2^32+1)^2, x^2 - 7), element = c0 + c1*x
        let cases: &[((u64,u64),(u64,u64),(u64,u64),(u64,u64),(u64,u64))] = &[
            // (a, b, a*b, a^2, a^-1)
            (
                (7545715477402998348, 14511670991623694876),
                (1987831781048552284, 7711293904151680427),
                (14683984057213515593, 13206543657924453654),
                (13045554628676694267, 10216530311954636359),
                (11415525327949330282, 17307565716904923058),
            ),
            (
                (880500214052543661, 4475494819024296210),
                (13429016922422520271, 5348018159025638973),
                (6884390467873344479, 8637939372288633548),
                (5996747058241258858, 6312094989139069224),
                (8649671574881506069, 12749859047259116419),
            ),
            (
                (4107359397723636147, 10750370325581404789),
                (12825561557292058881, 15560515093349535989),
                (5739795020184660833, 10814813959178477247),
                (3378973819513744472, 14743809710019143122),
                (16614819526829983410, 3310698401362174152),
            ),
        ];
        for &(a,b,ab,a2,ainv) in cases {
            let fa = ext2(a.0, a.1);
            let fb = ext2(b.0, b.1);
            assert_eq!(fa * fb, ext2(ab.0, ab.1), "mul mismatch");
            assert_eq!(fa.square(), ext2(a2.0, a2.1), "square mismatch");
            assert_eq!(fa.inverse().unwrap(), ext2(ainv.0, ainv.1), "inverse mismatch");
        }
    }
}
