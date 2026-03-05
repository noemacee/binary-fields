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

    /// Algorithm 2.49 — Binary algorithm for inversion in GF(2^128).
    ///
    /// # Comparison with Algorithm 2.48 (current invert)
    ///
    /// Algorithm 2.48 is the polynomial EEA: it works left-to-right, cancelling
    /// the highest-degree term of u at each step using a left-shift of v. This
    /// requires explicit degree computation and multi-word shifts (poly_shl).
    ///
    /// Algorithm 2.49 is the polynomial analogue of the binary GCD (Stein's
    /// algorithm). It works right-to-left, cancelling the lowest-degree term
    /// (the z^0 = constant term, i.e. bit 0) at each step. The operations are:
    ///
    ///   - "z divides u" ↔ bit 0 of u is 0  (u has no constant term)
    ///   - "u / z"       ↔ right-shift u by 1
    ///   - "deg(u) > deg(v)" can be replaced by u > v as integers (book §2.3.4)
    ///
    /// Because everything is right-shifts and bit-0 checks, no degree computation
    /// or variable-length shifts are needed — only fixed 1-bit shifts.
    ///
    /// # Representation
    ///
    /// u and v are polynomials that start at degree ≤ 128 (v starts as f of
    /// degree 128) and shrink. We store them as [u64; 3] to hold f initially.
    ///
    /// g1 and g2 are the Bezout coefficients, always kept reduced mod f, so
    /// they stay at degree < 128 and fit in [u64; 2].
    ///
    /// The invariants maintained throughout are:
    ///   g1 * a ≡ u  (mod f)
    ///   g2 * a ≡ v  (mod f)
    ///
    /// # The division step for g
    ///
    /// When we divide u by z, we must also divide g1 by z mod f.
    /// In GF(2)[z], dividing a polynomial g by z mod f:
    ///   - If bit 0 of g is 0: g/z is just g >> 1 (exact division, no reduction)
    ///   - If bit 0 of g is 1: (g + f) is divisible by z (since f has constant
    ///     term 1, g + f has bit 0 = 0), then (g + f)/z = (g + f) >> 1
    ///
    /// This keeps g1, g2 in the field (degree < 128, 2 words) throughout.
    pub fn invert_binary(self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // u and v: polynomials, stored as [u64; 3] since v starts as f (degree 128).
        // u starts as a (degree < 128), v starts as f (degree 128).
        let mut u = [self.0[0], self.0[1], 0u64];
        let mut v = F; // [lo, 0, 1] where lo = 0x87

        // g1, g2: Bezout coefficients, always degree < 128, stored as [u64; 2].
        // g1 starts as 1, g2 starts as 0.
        let mut g1 = [1u64, 0u64];
        let mut g2 = [0u64, 0u64];

        // Loop until u = 1 or v = 1.
        loop {
            // Check termination: u == 1 or v == 1
            if u[0] == 1 && u[1] == 0 && u[2] == 0 {
                return Some(GF2_128::new(g1[0], g1[1]));
            }
            if v[0] == 1 && v[1] == 0 && v[2] == 0 {
                return Some(GF2_128::new(g2[0], g2[1]));
            }

            // Step 3.1: while z divides u (bit 0 of u is 0), divide u and g1 by z.
            while u[0] & 1 == 0 {
                // u <- u / z: right shift u by 1 across 3 words
                u[0] = (u[0] >> 1) | (u[1] << 63);
                u[1] = (u[1] >> 1) | (u[2] << 63);
                u[2] >>= 1;

                // g1 <- g1 / z mod f
                if g1[0] & 1 == 0 {
                    // g1 is divisible by z: simple right shift
                    g1[0] = (g1[0] >> 1) | (g1[1] << 63);
                    g1[1] >>= 1;
                } else {
                    // g1 not divisible by z: add f first, then shift
                    // f mod z^128 = 0x87 (the low word of F; high word of g1+f
                    // may gain a bit from the carry, but degree stays < 128
                    // because after >> 1 the bit 128 of (g1+f) lands at bit 127)
                    let lo = g1[0] ^ F[0]; // XOR with low word of f (0x87)
                    let hi = g1[1] ^ F[1]; // XOR with middle word of f (0)
                    // Note: F[2]=1 represents z^128, but g1 has degree < 128,
                    // so g1 + f has degree exactly 128 (bit 128 = bit 0 of F[2]).
                    // After adding that z^128 term and shifting right by 1,
                    // it lands at bit 127 = the high bit of g1[1] after shift.
                    g1[0] = (lo >> 1) | (hi << 63);
                    g1[1] = (hi >> 1) | (1u64 << 63); // the z^128 bit lands here
                }
            }

            // Step 3.2: while z divides v (bit 0 of v is 0), divide v and g2 by z.
            while v[0] & 1 == 0 {
                // v <- v / z: right shift v by 1 across 3 words
                v[0] = (v[0] >> 1) | (v[1] << 63);
                v[1] = (v[1] >> 1) | (v[2] << 63);
                v[2] >>= 1;

                // g2 <- g2 / z mod f (same logic as g1 above)
                if g2[0] & 1 == 0 {
                    g2[0] = (g2[0] >> 1) | (g2[1] << 63);
                    g2[1] >>= 1;
                } else {
                    let lo = g2[0] ^ F[0];
                    let hi = g2[1] ^ F[1];
                    g2[0] = (lo >> 1) | (hi << 63);
                    g2[1] = (hi >> 1) | (1u64 << 63);
                }
            }

            // Step 3.3: compare u and v as integers (valid proxy for degree comparison).
            // If u > v: u <- u + v, g1 <- g1 + g2
            // Else:     v <- v + u, g2 <- g2 + g1
            //
            // Compare from most significant word down.
            let u_gt_v = if u[2] != v[2] {
                u[2] > v[2]
            } else if u[1] != v[1] {
                u[1] > v[1]
            } else {
                u[0] > v[0]
            };

            if u_gt_v {
                u[0] ^= v[0];
                u[1] ^= v[1];
                u[2] ^= v[2];
                g1[0] ^= g2[0];
                g1[1] ^= g2[1];
            } else {
                v[0] ^= u[0];
                v[1] ^= u[1];
                v[2] ^= u[2];
                g2[0] ^= g1[0];
                g2[1] ^= g1[1];
            }
        }
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

    #[test]
    fn invert_binary_zero_is_none() {
        assert_eq!(GF2_128::zero().invert_binary(), None);
    }

    #[test]
    fn invert_binary_one_is_one() {
        assert_eq!(GF2_128::one().invert_binary(), Some(GF2_128::one()));
    }

    #[test]
    fn invert_binary_roundtrip() {
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);
        assert_eq!(a.mul(a.invert_binary().unwrap()), GF2_128::one());
    }

    #[test]
    fn invert_binary_agrees_with_2_48() {
        let cases = [
            GF2_128::new(0x1234, 0x5678),
            GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234),
            GF2_128::new(0xffffffffffffffff, 0xffffffffffffffff),
            GF2_128::new(0x1, 0x0),
        ];
        for a in cases {
            assert_eq!(a.invert_binary(), a.invert(), "mismatch for {:?}", a);
        }
    }
}
