//! Optimized algorithms for f(z) = z^128 + z^7 + z^2 + z + 1 (GCM polynomial).
//!
//! Contains only algorithms specific to this polynomial.
//! Generic algorithms (2.33, 2.40, 2.48, 2.49) live in `crate::generic`.

/// Algorithm 2.41 — Fast word-level reduction modulo f(z) = z^128 + z^7 + z^2 + z + 1.
///
/// Reduces 14 XOR/shift operations versus 127 loop iterations in Algorithm 2.40.
/// INPUT:  c(z) of degree at most 254, stored as [u64; 4].
/// OUTPUT: c(z) mod f(z), stored as [u64; 2].
pub fn reduce_2_41(c: [u64; 4]) -> [u64; 2] {
    let mut c = c;

    let t = c[3];
    c[1] ^= t;
    c[1] ^= t << 1;
    c[1] ^= t << 2;
    c[1] ^= t << 7;
    c[2] ^= t >> 63;
    c[2] ^= t >> 62;
    c[2] ^= t >> 57;

    let t = c[2];
    c[0] ^= t;
    c[0] ^= t << 1;
    c[0] ^= t << 2;
    c[0] ^= t << 7;
    c[1] ^= t >> 63;
    c[1] ^= t >> 62;
    c[1] ^= t >> 57;

    [c[0], c[1]]
}

/// Algorithm 2.34 — Right-to-left comb multiplication + Algorithm 2.41 fast reduction.
pub fn mul_2_34(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    let mut b_ext = [b[0], b[1], 0u64];
    let mut c = [0u64; 4];
    for k in 0..64usize {
        for j in 0..2usize {
            if (a[j] >> k) & 1 == 1 {
                c[j] ^= b_ext[0];
                c[j + 1] ^= b_ext[1];
                c[j + 2] ^= b_ext[2];
            }
        }
        if k != 63 {
            let carry_01 = b_ext[0] >> 63;
            let carry_12 = b_ext[1] >> 63;
            b_ext[0] <<= 1;
            b_ext[1] = (b_ext[1] << 1) | carry_01;
            b_ext[2] = (b_ext[2] << 1) | carry_12;
        }
    }
    reduce_2_41(c)
}

/// Algorithm 2.39 — Squaring via bit-expansion table + Algorithm 2.41 fast reduction.
///
/// T[d] is the 16-bit expansion with a 0 bit inserted between each bit of d.
pub fn square_2_39(a: [u64; 2]) -> [u64; 2] {
    const SQUARE_TABLE: [u16; 256] = {
        let mut table = [0u16; 256];
        let mut d = 0usize;
        while d < 256 {
            let mut result = 0u16;
            let mut i = 0;
            while i < 8 {
                if (d >> i) & 1 == 1 { result |= 1 << (2 * i); }
                i += 1;
            }
            table[d] = result;
            d += 1;
        }
        table
    };

    let mut expanded = [0u64; 4];
    for (word_idx, &w) in a.iter().enumerate() {
        let mut lo = 0u128;
        for i in 0..8 {
            lo |= (SQUARE_TABLE[((w >> (i * 8)) & 0xff) as usize] as u128) << (i * 16);
        }
        expanded[word_idx * 2] = lo as u64;
        expanded[word_idx * 2 + 1] = (lo >> 64) as u64;
    }
    reduce_2_41(expanded)
}

/// Exponentiation by square-and-multiply (right-to-left), using Algorithm 2.34 + 2.39.
pub fn pow_rtl_2_34_2_39(a: [u64; 2], n: u128) -> [u64; 2] {
    if n == 0 { return [1, 0]; }
    let mut result = [1u64, 0u64];
    let mut base = a;
    let mut exp = n;
    while exp > 0 {
        if exp & 1 == 1 { result = mul_2_34(result, base); }
        base = square_2_39(base);
        exp >>= 1;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generic::arithmetic::mul_2_33;
    use crate::generic::invert::invert_2_48;

    const GCM_POLY: [u64; 2] = [0x87, 0];

    #[test]
    fn reduce_zero() { assert_eq!(reduce_2_41([0, 0, 0, 0]), [0, 0]); }

    #[test]
    fn reduce_z128() { assert_eq!(reduce_2_41([0, 0, 1, 0]), [0x87, 0]); }

    #[test]
    fn reduce_z192() { assert_eq!(reduce_2_41([0, 0, 0, 1]), [0, 0x87]); }

    #[test]
    fn reduce_z254() {
        assert_eq!(reduce_2_41([0, 0, 0, 1 << 62]), [0x1067, 0xc000000000000000]);
    }

    #[test]
    fn mul_comb_agrees_with_shift_add() {
        let pairs: &[([u64; 2], [u64; 2])] = &[
            ([0x1234, 0], [0x5678, 0]),
            ([0xdeadbeefcafe1234, 0], [0xabcd1234, 0]),
            ([0, 1u64 << 63], [0b10, 0]),
            ([0xffffffffffffffff, 0xffffffffffffffff], [1, 0]),
        ];
        for &(a, b) in pairs {
            assert_eq!(mul_2_34(a, b), mul_2_33(&a, &b, &GCM_POLY, 128));
        }
    }

    #[test]
    fn square_consistent_with_mul() {
        let a = [0x1234_u64, 0x5678];
        assert_eq!(square_2_39(a), mul_2_34(a, a));
    }

    #[test]
    fn pow_zero_is_one() {
        assert_eq!(pow_rtl_2_34_2_39([0xdeadbeef, 0xcafe], 0), [1, 0]);
    }

    #[test]
    fn pow_one_is_self() {
        let a = [0xdeadbeef_u64, 0xcafe];
        assert_eq!(pow_rtl_2_34_2_39(a, 1), a);
    }

    #[test]
    fn pow_two_is_square() {
        let a = [0x1234_u64, 0x5678];
        assert_eq!(pow_rtl_2_34_2_39(a, 2), mul_2_33(&a, &a, &GCM_POLY, 128));
    }

    #[test]
    fn pow_fermat() {
        let a = [0xdeadbeefcafe1234_u64, 0xabcd1234];
        assert_eq!(pow_rtl_2_34_2_39(a, u128::MAX), [1, 0]);
    }

    #[test]
    fn pow_inverse_consistent() {
        let a = [0x1234_u64, 0x5678];
        assert_eq!(pow_rtl_2_34_2_39(a, u128::MAX - 1), invert_2_48(&a, &GCM_POLY, 128).unwrap());
    }
}
