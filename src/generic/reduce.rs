/// Algorithm 2.40 — Generic bit-by-bit modular reduction for GF(2^degree).
///
/// Works for any irreducible polynomial and any number of limbs N.
/// Parameterized by the field degree and the non-leading terms of f(z).
///
/// # Arguments
/// - `c`: mutable slice of exactly `2 * poly.len()` limbs (unreduced product).
/// - `poly`: non-leading terms of f(z) packed as `[u64; N]`.
///   E.g. for f = z^128 + z^7 + z^2 + z + 1: poly = [0x87, 0].
/// - `degree`: the extension degree m.
///
/// After the call `c[0..N]` holds the reduced result.
pub fn reduce_2_40(c: &mut [u64], poly: &[u64], degree: usize) {
    let n = poly.len();
    debug_assert_eq!(c.len(), 2 * n);

    for i in (degree..=(2 * degree - 2)).rev() {
        let word = i / 64;
        let bit = i % 64;
        if (c[word] >> bit) & 1 == 0 {
            continue;
        }
        let shift = i - degree;
        let w_shift = shift / 64;
        let b_shift = shift % 64;
        for j in 0..n {
            let dst = j + w_shift;
            if dst < c.len() {
                c[dst] ^= poly[j] << b_shift;
            }
            if b_shift > 0 {
                let dst1 = dst + 1;
                if dst1 < c.len() {
                    c[dst1] ^= poly[j] >> (64 - b_shift);
                }
            }
        }
    }

    // Mask dirty bits for non-word-boundary degrees.
    let top_bit = degree % 64;
    if top_bit != 0 {
        c[degree / 64] &= (1u64 << top_bit) - 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GCM_POLY: [u64; 2] = [0x87, 0];

    #[test]
    fn reduce_zero() {
        let mut c = [0u64; 4];
        reduce_2_40(&mut c, &GCM_POLY, 128);
        assert_eq!(c, [0, 0, 0, 0]);
    }

    #[test]
    fn reduce_one() {
        let mut c = [1u64, 0, 0, 0];
        reduce_2_40(&mut c, &GCM_POLY, 128);
        assert_eq!(&c[..2], &[1u64, 0]);
    }

    #[test]
    fn reduce_z128() {
        let mut c = [0u64, 0, 1, 0];
        reduce_2_40(&mut c, &GCM_POLY, 128);
        assert_eq!(&c[..2], &[0x87u64, 0]);
    }

    #[test]
    fn reduce_z192() {
        let mut c = [0u64, 0, 0, 1];
        reduce_2_40(&mut c, &GCM_POLY, 128);
        assert_eq!(&c[..2], &[0u64, 0x87]);
    }

    #[test]
    fn reduce_z254() {
        let mut c = [0u64, 0, 0, 1u64 << 62];
        reduce_2_40(&mut c, &GCM_POLY, 128);
        assert_eq!(&c[..2], &[0x1067u64, 0xc000000000000000u64]);
    }
}
