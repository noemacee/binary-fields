//! `Gf233Config` — config for GF(2^233) with irreducible polynomial
//! f(z) = z^233 + z^74 + 1  (NIST B-233 / K-233 field polynomial).
//!
//! This config provides only the three required constants (`ALPHA_POW_M`, `ZERO`, `ONE`)
//! and relies entirely on the generic default implementations of `mul_assign`,
//! `square_in_place`, and `inverse` from `BinaryFieldConfig`.  No
//! polynomial-specific optimizations are needed — this demonstrates that
//! the generic layer is sufficient to instantiate a fully working field.

use crate::ark::{BinaryField, BinaryFieldConfig};

/// Config for GF(2^233).
pub struct Gf233Config;

impl BinaryFieldConfig<4> for Gf233Config {
    /// 233 bits fit in 4 × 64-bit limbs (256 bits), with 23 bits unused in limb 3.
    const DEGREE: usize = 233;

    /// α^233 in the polynomial basis, derived from f(α) = 0:
    ///   f(z) = z^233 + z^74 + 1  →  α^233 = α^74 + 1
    ///   α^0  → bit  0 of limb 0  →  ALPHA_POW_M[0] = 1
    ///   α^74 → bit 10 of limb 1 (74 = 64 + 10)  →  ALPHA_POW_M[1] = 1 << 10
    const ALPHA_POW_M: [u64; 4] = [1, 1 << 10, 0, 0];

    const ZERO: BinaryField<Self, 4> = BinaryField([0u64; 4], core::marker::PhantomData);

    const ONE: BinaryField<Self, 4> = BinaryField([1, 0, 0, 0], core::marker::PhantomData);

    // No method overrides: mul_assign, square_in_place, and inverse all use
    // the generic implementations driven by ALPHA_POW_M and DEGREE above.
}

/// GF(2^233) as an arkworks `Field`.
pub type Gf233 = BinaryField<Gf233Config, 4>;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, One, Zero};

    fn g(lo: u64) -> Gf233 {
        Gf233::new([lo, 0, 0, 0])
    }

    #[test]
    fn zero_is_zero() {
        assert!(Gf233::zero().is_zero());
    }

    #[test]
    fn one_is_one() {
        assert!(Gf233::one().is_one());
    }

    #[test]
    fn add_self_is_zero() {
        let a = Gf233::new([0x1234, 0x5678, 0, 0]);
        assert!((a + a).is_zero());
    }

    #[test]
    fn mul_by_one_is_identity() {
        let a = Gf233::new([0xdeadbeef, 0xabcd, 0, 0]);
        assert_eq!(a * Gf233::one(), a);
    }

    #[test]
    fn mul_by_zero_is_zero() {
        let a = Gf233::new([0xdeadbeef, 0xabcd, 0, 0]);
        assert!((a * Gf233::zero()).is_zero());
    }

    #[test]
    fn mul_commutative() {
        let a = Gf233::new([0x1234, 0x5678, 0, 0]);
        let b = Gf233::new([0xabcd, 0xef01, 0, 0]);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn mul_associative() {
        let a = g(0x1234);
        let b = g(0x5678);
        let c = g(0xabcd);
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn mul_distributive() {
        let a = g(0x1234);
        let b = g(0x5678);
        let c = g(0xabcd);
        assert_eq!(a * (b + c), a * b + a * c);
    }

    #[test]
    fn square_consistent_with_mul() {
        let a = Gf233::new([0x1234, 0x5678, 0, 0]);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn inverse_of_zero_is_none() {
        assert!(Gf233::zero().inverse().is_none());
    }

    #[test]
    fn inverse_of_one_is_one() {
        assert_eq!(Gf233::one().inverse(), Some(Gf233::one()));
    }

    #[test]
    fn inverse_roundtrip() {
        let a = Gf233::new([0xdeadbeef, 0xabcd1234, 0, 0]);
        let a_inv = a.inverse().unwrap();
        assert_eq!(a * a_inv, Gf233::one());
    }

    #[test]
    fn inverse_roundtrip_high_limbs() {
        // Element with bits set in limbs 1 and 2 (but within DEGREE=233 bits).
        let a = Gf233::new([0xdeadbeef, 0x12345678, 0xabcd, 0]);
        let a_inv = a.inverse().unwrap();
        assert_eq!(a * a_inv, Gf233::one());
    }

    #[test]
    fn characteristic_is_2() {
        assert!((Gf233::one() + Gf233::one()).is_zero());
    }

    #[test]
    fn extension_degree_is_233() {
        assert_eq!(Gf233::extension_degree(), 233);
    }

    // ── Known-vector tests generated with the Python `galois` library ────────
    // seed 0xdeadbeef; galois.GF(2**233, irreducible_poly="x^233+x^74+1")

    #[test]
    fn mul_known_vectors() {
        // Generated with Python galois: GF(2**233, irreducible_poly='x^233+x^74+1')
        type L = [u64; 4];
        let cases: &[(L,L,L)] = &[
            ([0x5f6a04567a20f9e7, 0x77dcf44bb14a97a9, 0x22f738320b02af21, 0x000000e2aa5a1584],
             [0x23500bd63d1b066d, 0xec5a0cca4e86c0d8, 0xef402a6920736baa, 0x000001ed706336c9],
             [0x4f3b4cc4507c62d5, 0xec684e0a4d8d64ce, 0x5bff407911e27ea1, 0x0000004a1c33218a]),
            ([0xb264beb461e648c3, 0x1ca0a5672911970c, 0xac09b9d860192c7e, 0x0000001294e06c6c],
             [0x2d90d6b154797aa2, 0x6addf44af94f7d8c, 0xa6f41b0146d0fd30, 0x0000005d88706d36],
             [0x6fafc9912d5db8bb, 0x1ced7526661a67e8, 0xdf02efc9dae31581, 0x0000018313759f36]),
            ([0x1bf40996f816cdbe, 0xc2e3ace50da6cbc1, 0xdc6f2321c2ec5567, 0x000001afb591a180],
             [0x5f30251c87006405, 0x2650262a5a8551a2, 0x7d6de5a83a6c692e, 0x000000db1bff01e3],
             [0xfd60d3eb796f474e, 0x9c723636a327d8bd, 0xb8f591c78b47c0c7, 0x00000120b29245ff]),
            ([0x6b69dc01fc911746, 0x7065b82fd8a0a075, 0xe443e17a74a3f087, 0x000000f29dd55b57],
             [0xbf65c1b1ecef6800, 0xd0d3ecd3f44beeee, 0x3a22b66bbe1f3bef, 0x0000011b76aaeb9f],
             [0x7af2f89c176f8c03, 0xe0639d72dc69ea6e, 0xac8f235b732a08c2, 0x0000009997fcfa4e]),
            ([0x019a5cc5fb91e5c3, 0x27a7e320625b2746, 0xd7b41b8327b2cb2b, 0x0000010d5d25fe0d],
             [0x5684aa7e389e6e0b, 0xd68a9fa7e87c8c1b, 0xf4f323422fc18905, 0x000001bcbb953e21],
             [0x5965972f4f9aecd8, 0xd04a2600929122b8, 0x8fb8c6cd5dcd18e8, 0x000000351f5e4664]),
            ([0xa551a9c7d62e5e0a, 0xaec87f8c61aaf5b5, 0x16b47d383e7cc2ba, 0x000000a3fca6be28],
             [0xd2dcf8dff8cbabb2, 0xf3ab3e045c03065e, 0xc89945a2bbf2814c, 0x00000177b672c4a5],
             [0x9308177c30c9235e, 0xc45f5d416a858f19, 0x7c4410428bd49e13, 0x0000008e34211fb0]),
        ];
        for &(a,b,exp) in cases {
            assert_eq!(Gf233::new(a) * Gf233::new(b), Gf233::new(exp));
        }
    }

    #[test]
    fn square_known_vectors() {
        // Generated with Python galois.
        type L = [u64; 4];
        let cases: &[(L,L)] = &[
            ([0x94932d53a025641f, 0xb24a9e8c0eaf0361, 0x89f3c2c597c907b8, 0x000001e992d647e1],
             [0x25c10ae2a0229609, 0xc7838dfdc6aceaf9, 0x2a54eedd828ac52c, 0x000000ece1766872]),
            ([0x8bbbb12ada2fb268, 0x1bc37cf287e76c27, 0x5fce7c1592648f06, 0x0000010d1fbe92b7],
             [0x08657b901b45b2ad, 0xefb0f83ccb3a334e, 0xc095d4151cd487d9, 0x000000af1f7a5f0e]),
            ([0x4b619e14bb709bf1, 0x206655417dd39182, 0x7acc0c80d7724c31, 0x00000123097803fa],
             [0xc24d7c1781c7505e, 0x19e163234969f69e, 0xb5f1d3af4125d830, 0x0000001e3919902b]),
            ([0x53b054d467067786, 0x18f282c1e9d157fa, 0x3667888bbdf007dd, 0x0000006ca551b974],
             [0xc4511badbcd0a12b, 0x0bcb8387243605ae, 0x5463710113b5240f, 0x0000012e620e70a3]),
            ([0x50e3ef62d9b413cf, 0x9176a891e0e492db, 0x65456a2b62ecbaf9, 0x000001756c8cc1b1],
             [0x6877b77af1f5d19a, 0xc243ddc4b4d82b3b, 0x54a05c38e982a529, 0x0000019ce4e84903]),
            ([0x9354a310681b3626, 0xfdf0ad97df723c6f, 0x8136d48ec79af1da, 0x00000081cbd1e071],
             [0xb73ac4bca66aacab, 0x203a7fbd948ff7e0, 0x73df95a625db965a, 0x00000122c4fb6995]),
        ];
        for &(a,exp) in cases {
            assert_eq!(Gf233::new(a).square(), Gf233::new(exp));
        }
    }

    #[test]
    fn inverse_known_vectors() {
        // Generated with Python galois.
        type L = [u64; 4];
        let cases: &[(L,L)] = &[
            ([0x79e54327663bf248, 0xa72247efcb3a7e61, 0x5dd42821a04f30eb, 0x000001c645923bd5],
             [0x4037e829a230bae5, 0xb57a29903a25f6d8, 0x8d6b106563a2120b, 0x00000042fb626806]),
            ([0xe8e6699c89531818, 0x2c174dfde122d54d, 0x0f7804ff5aa75144, 0x0000017de9d0b305],
             [0xa99e37079dc2775e, 0xd1fe05beddf1fdb4, 0x3b198e30b0c467cf, 0x00000172eaa32878]),
            ([0x13d3b1540200d51a, 0xdad475376b2bd6c8, 0x7c9591857e4b0f46, 0x00000118a16db5cf],
             [0x06c57a47607d229e, 0xe3d6688afa3fcbf9, 0x0110355b5e52515a, 0x000000b6e884f0d4]),
            ([0xde7daba03f3d433a, 0x0012f351dcbe8cbb, 0x7e8831d7ca84f6e6, 0x000000b56d7b6548],
             [0x2696ba66d5aacd24, 0xbbb1f115f802de0d, 0x9e411d4bfef4b78a, 0x000000a3f34f51dd]),
            ([0x7ee2c6a6b9e3b308, 0x2b7e8d7036041bbd, 0x45c2b200eb15b7ef, 0x000001c1fb7f8970],
             [0xa0a1c5678f99ffe3, 0x7984dbadd9d44cf1, 0xdf621dc8aba97bd8, 0x00000031f247f4b9]),
            ([0xd132072ef366ff8a, 0x2cf090dcedabd8f7, 0x0eecba504ebe58c7, 0x00000122d74e4c0c],
             [0x990e209c2389de8d, 0xab81bef66722c75c, 0x1e8d3ca737491db7, 0x000001dce12e1e97]),
        ];
        for &(a,exp) in cases {
            assert_eq!(Gf233::new(a).inverse(), Some(Gf233::new(exp)));
        }
    }

    #[test]
    fn mul_reduction_basic() {
        // For any nonzero a: (a * b) * a^{-1} = b — exercises reduction.
        let a = Gf233::new([0xdeadbeef, 0, 0, 0]);
        let b = Gf233::new([0x12345678, 0, 0, 0]);
        let product = a * b;
        let recovered = product * a.inverse().unwrap();
        assert_eq!(recovered, b);
    }

    #[test]
    fn invert_2_49_roundtrip() {
        use crate::generic::invert::invert_2_49;
        use crate::generic::arithmetic::mul_2_33;

        const GF233_POLY: [u64; 4] = [1, 1 << 10, 0, 0];
        const DEGREE_233: usize = 233;

        let a = [0xdeadbeef_abcd1234u64, 0x12345678u64, 0xef01u64, 0u64];
        let inv = invert_2_49(&a, &GF233_POLY, DEGREE_233).unwrap();
        assert_eq!(mul_2_33(&a, &inv, &GF233_POLY, DEGREE_233), [1, 0, 0, 0]);
    }

    #[test]
    fn invert_2_49_agrees_with_2_48() {
        use crate::generic::invert::{invert_2_48, invert_2_49};

        const GF233_POLY: [u64; 4] = [1, 1 << 10, 0, 0];
        const DEGREE_233: usize = 233;

        let cases: &[[u64; 4]] = &[
            [0xdeadbeef, 0, 0, 0],
            [0xdeadbeef_abcd1234, 0x12345678, 0xef01, 0],
            [1, 0, 0, 0],
        ];
        for a in cases {
            assert_eq!(
                invert_2_49(a, &GF233_POLY, DEGREE_233),
                invert_2_48(a, &GF233_POLY, DEGREE_233),
            );
        }
    }
}
