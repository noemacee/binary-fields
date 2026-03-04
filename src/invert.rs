use crate::GF2_128;

/// The irreducible polynomial f(z) = z^128 + z^7 + z^2 + z + 1
/// Stored as [u64; 3] for use in the EEA where we need degree 128.
/// Bit 128 is in word index 2, bit 0 of that word.
const F: [u64; 3] = [
    0x0000_0000_0000_0087, // z^7 + z^2 + z + 1  (low word)
    0x0000_0000_0000_0000, // (middle word)
    0x0000_0000_0000_0001, // z^128              (high word, bit 0)
];

/// Returns the degree of a polynomial stored as [u64; 3].
/// Degree = position of the highest set bit.
/// Returns -1 if the polynomial is zero.
fn deg(a: [u64; 3]) -> i32 {
    for i in (0..3).rev() {
        if a[i] != 0 {
            return (i as i32) * 64 + (63 - a[i].leading_zeros() as i32);
        }
    }
    -1
}

/// Left shift a [u64; 3] polynomial by j bits (multiply by z^j).
/// Assumes j < 192 (we never need more than that here).
fn poly_shl(a: [u64; 3], j: u32) -> [u64; 3] {
    if j == 0 {
        return a;
    }
    let word_shift = (j / 64) as usize;
    let bit_shift = j % 64;

    let mut result = [0u64; 3];
    for i in 0..3 {
        if i + word_shift < 3 {
            result[i + word_shift] |= a[i] << bit_shift;
        }
        if bit_shift > 0 && i + word_shift + 1 < 3 {
            result[i + word_shift + 1] |= a[i] >> (64 - bit_shift);
        }
    }
    result
}

/// XOR two [u64; 3] polynomials.
fn poly_xor(a: [u64; 3], b: [u64; 3]) -> [u64; 3] {
    [a[0] ^ b[0], a[1] ^ b[1], a[2] ^ b[2]]
}

/// Convert a GF2_128 element to a [u64; 3] polynomial.
fn to_poly(a: GF2_128) -> [u64; 3] {
    [a.0[0], a.0[1], 0]
}

/// Convert a [u64; 3] polynomial back to GF2_128.
/// The high word must be zero — i.e. degree < 128.
fn to_field(a: [u64; 3]) -> GF2_128 {
    GF2_128::new(a[0], a[1])
}

impl GF2_128 {
    /// Algorithm 2.48 — Inversion in GF(2^128) using the Extended Euclidean Algorithm.
    /// Returns None if self is zero (zero has no inverse).
    pub fn invert(self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let mut u = to_poly(self);
        let mut v = F;
        let mut g1 = to_poly(GF2_128::one());
        let mut g2 = to_poly(GF2_128::zero());

        while deg(u) != 0 {
            let j = deg(u) - deg(v);
            if j < 0 {
                std::mem::swap(&mut u, &mut v);
                std::mem::swap(&mut g1, &mut g2);
                let j = -j;
                u = poly_xor(u, poly_shl(v, j as u32));
                g1 = poly_xor(g1, poly_shl(g2, j as u32));
            } else {
                u = poly_xor(u, poly_shl(v, j as u32));
                g1 = poly_xor(g1, poly_shl(g2, j as u32));
            }
        }

        Some(to_field(g1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invert_zero_is_none() {
        assert_eq!(GF2_128::zero().invert(), None);
    }

    #[test]
    fn invert_one_is_one() {
        assert_eq!(GF2_128::one().invert(), Some(GF2_128::one()));
    }

    #[test]
    fn invert_roundtrip() {
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);
        let a_inv = a.invert().unwrap();
        assert_eq!(a.mul(a_inv), GF2_128::one());
    }

    #[test]
    fn invert_another() {
        let a = GF2_128::new(0x1234, 0x5678);
        let a_inv = a.invert().unwrap();
        assert_eq!(a.mul(a_inv), GF2_128::one());
    }
}
