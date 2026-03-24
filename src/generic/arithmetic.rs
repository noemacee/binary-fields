use super::reduce::reduce_2_40;

/// Algorithm 2.32 — Addition in GF(2^m): bitwise XOR, independent of the polynomial.
pub fn add_2_32<const N: usize>(a: [u64; N], b: [u64; N]) -> [u64; N] {
    let mut result = a;
    for i in 0..N {
        result[i] ^= b[i];
    }
    result
}

/// Algorithm 2.33 — Right-to-left shift-and-add multiplication in GF(2^degree).
///
/// Generic carry-less multiplication for any field degree and irreducible polynomial.
/// Uses bit-by-bit reduction (Algorithm 2.40) internally.
pub fn mul_2_33<const N: usize>(
    a: &[u64; N],
    b: &[u64; N],
    poly: &[u64; N],
    degree: usize,
) -> [u64; N] {
    let mut c = vec![0u64; 2 * N];

    for i in 0..degree {
        let word = i / 64;
        let bit_pos = i % 64;
        if (a[word] >> bit_pos) & 1 == 0 {
            continue;
        }
        let w_shift = i / 64;
        let b_shift = i % 64;
        for j in 0..N {
            let dst = j + w_shift;
            c[dst] ^= b[j] << b_shift;
            if b_shift > 0 && dst + 1 < 2 * N {
                c[dst + 1] ^= b[j] >> (64 - b_shift);
            }
        }
    }

    reduce_2_40(&mut c, poly, degree);

    let mut result = [0u64; N];
    result.copy_from_slice(&c[..N]);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    const GCM_POLY: [u64; 2] = [0x87, 0];
    const DEGREE_128: usize = 128;

    #[test]
    fn add_zero_is_identity() {
        let a = [0xdeadbeef_u64, 0xcafe];
        assert_eq!(add_2_32(a, [0, 0]), a);
    }

    #[test]
    fn add_self_is_zero() {
        let a = [0xdeadbeef_u64, 0xcafe];
        assert_eq!(add_2_32(a, a), [0, 0]);
    }

    #[test]
    fn mul_one_identity() {
        let a: [u64; 2] = [0xdeadbeef, 0xcafe];
        assert_eq!(mul_2_33(&a, &[1, 0], &GCM_POLY, DEGREE_128), a);
    }

    #[test]
    fn mul_zero_is_zero() {
        let a: [u64; 2] = [0xdeadbeef, 0xcafe];
        assert_eq!(mul_2_33(&a, &[0, 0], &GCM_POLY, DEGREE_128), [0, 0]);
    }

    #[test]
    fn mul_commutative() {
        let a: [u64; 2] = [0x1234, 0x5678];
        let b: [u64; 2] = [0xabcd, 0xef01];
        assert_eq!(mul_2_33(&a, &b, &GCM_POLY, DEGREE_128), mul_2_33(&b, &a, &GCM_POLY, DEGREE_128));
    }

    #[test]
    fn mul_reduction_z128() {
        // z^127 * z = z^128 ≡ z^7 + z^2 + z + 1 = 0x87
        let z127: [u64; 2] = [0, 1u64 << 63];
        let z1: [u64; 2] = [0b10, 0];
        assert_eq!(mul_2_33(&z127, &z1, &GCM_POLY, DEGREE_128), [0x87, 0]);
    }

    #[test]
    fn mul_known_vector() {
        let a: [u64; 2] = [0x1234, 0];
        let b: [u64; 2] = [0x5678, 0];
        assert_eq!(mul_2_33(&a, &b, &GCM_POLY, DEGREE_128), [0x5c58160, 0]);
    }
}
