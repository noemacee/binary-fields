use crate::GF2_128;

/// The irreducible polynomial f(z) = z^128 + z^7 + z^2 + z + 1 : The one used in Galois Counter Mode (GCM)
/// Since f(z) = 0 in the field, z^128 = z^7 + z^2 + z + 1.
/// r(z) = z^7 + z^2 + z + 1 = 0b10000111 = 0x87
/// is what we XOR with when the high bit (z^128) spills over.
const R: u64 = 0x87;

/// Algorithm 2.39 — Polynomial squaring via byte expansion table.
/// Precomputed table T: for each byte d, T[d] is the 16-bit expansion
/// obtained by inserting a 0 bit between each bit of d.
/// T(d) = (0,d7,0,d6,0,d5,0,d4,0,d3,0,d2,0,d1,0,d0)
const SQUARE_TABLE: [u16; 256] = {
    let mut table = [0u16; 256];
    let mut d = 0usize;
    while d < 256 {
        let mut result = 0u16;
        let mut i = 0;
        while i < 8 {
            if (d >> i) & 1 == 1 {
                result |= 1 << (2 * i);
            }
            i += 1;
        }
        table[d] = result;
        d += 1;
    }
    table
};

/// Expand a single 64-bit word into a 128-bit result by inserting
/// a 0 bit between each pair of consecutive bits.
/// Each byte is expanded via SQUARE_TABLE to a 16-bit value.
fn expand_word(w: u64) -> u128 {
    let mut result = 0u128;
    for i in 0..8 {
        let byte = (w >> (i * 8)) as u8;
        let expanded = SQUARE_TABLE[byte as usize] as u128;
        result |= expanded << (i * 16);
    }
    result
}

impl GF2_128 {
    /// Algorithm 2.32 — Addition in GF(2^128)
    /// Addition of field elements is performed bitwise (XOR), word by word.
    pub fn add_2_32(self, rhs: Self) -> Self {
        GF2_128::new(self.0[0] ^ rhs.0[0], self.0[1] ^ rhs.0[1])
    }

    /// Subtraction is identical to addition in characteristic 2.
    pub fn sub(self, rhs: Self) -> Self {
        self.add_2_32(rhs)
    }

    /// Negation is the identity in characteristic 2.
    pub fn neg(self) -> Self {
        self
    }

    /// Multiply self by z mod f(z).
    /// Algorithm: left-shift by 1. If the high bit (a127) was 1,
    /// XOR the low word with R to reduce mod f(z).
    fn mul_by_z(self) -> Self {
        // Extract the highest bit of the high word before shifting
        let high_bit = self.0[1] >> 63;

        // Carry bit 63 of A[0] into bit 0 of A[1]
        // A[0] holds LOW degree terms, A[1] holds HIGH degree terms
        // shifting left (toward higher degree) means:
        //   bit 63 of A[0] carries into bit 0 of A[1]
        let carry = self.0[0] >> 63;
        let lo = self.0[0] << 1;
        let hi = (self.0[1] << 1) | carry;

        // If a127 (the old high bit of A[1]) was 1, reduce mod f(z)
        // z^128 = z^7 + z^2 + z + 1 = 0x87
        let lo = if high_bit == 1 { lo ^ R } else { lo };

        GF2_128::new(lo, hi)
    }

    /// Algorithm 2.33 — Right-to-left shift-and-add field multiplication.
    /// Computes a(z) * b(z) mod f(z).
    pub fn mul_2_33(self, rhs: Self) -> Self {
        let a = self;
        let mut b = rhs;
        let mut c = if a.0[0] & 1 == 1 {
            b // a0 = 1: initialize accumulator to b
        } else {
            GF2_128::zero()
        };

        for i in 1..128 {
            b = b.mul_by_z();
            // Extract bit i of a: which word and which bit within that word
            let word = i / 64;
            let bit = i % 64;
            if (a.0[word] >> bit) & 1 == 1 {
                c = c.add_2_32(b);
            }
        }

        c
    }

    /// Squaring via mul_2_33: computes self * self using Algorithm 2.33.
    /// Serves as a baseline to compare against the dedicated square_2_39.
    pub fn square_via_mul_2_33(self) -> Self {
        self.mul_2_33(self)
    }

    /// Algorithm 2.39 — Squaring in GF(2^128).
    /// Expands bits then reduces mod f(z) via Algorithm 2.40.
    pub fn square_2_39(self) -> Self {
        // Step 1: expand each word into 128 bits
        // A[0] (low word) expands into the low 128 bits of c
        // A[1] (high word) expands into the high 128 bits of c
        let lo = expand_word(self.0[0]);
        let hi = expand_word(self.0[1]);

        // Assemble into [u64; 4] intermediate
        let c: [u64; 4] = [lo as u64, (lo >> 64) as u64, hi as u64, (hi >> 64) as u64];

        // Step 2: reduce mod f(z)
        let reduced = crate::reduce::reduce_2_41(c);
        GF2_128::new(reduced[0], reduced[1])
    }

    /// Algorithm 2.34 — Right-to-left comb method for polynomial multiplication.
    ///
    /// # Difference from Algorithm 2.33 (current mul)
    ///
    /// Algorithm 2.33 processes A one BIT at a time from right to left:
    ///   for each bit i of A (0..127):
    ///     if bit i is set: C ^= B
    ///     B <- B * z
    /// That's 127 single-bit shifts of B, each touching 2 words.
    ///
    /// Algorithm 2.34 instead groups bits by their position WITHIN a word (column k),
    /// then processes all words of A for that column in the inner loop.
    /// B is still shifted one bit per outer iteration, but there are only W=64
    /// outer iterations instead of m=128, because the inner loop handles the
    /// word-offset positioning via C{j} (adding B into C starting at word j).
    ///
    /// # The C{j} notation
    ///
    /// C{j} means "the array C viewed starting at index j". So adding B into C{j}
    /// means: C[j] ^= B[0], C[j+1] ^= B[1], C[j+2] ^= B[2].
    /// This is how the word offset z^(W*j) is handled without an actual shift.
    ///
    /// # Cost comparison
    ///
    /// Algorithm 2.33: 127 B-shifts (2 words each) + up to 128 XOR-adds (2 words each)
    /// Algorithm 2.34:  63 B-shifts (3 words each) + up to 128 XOR-adds (3 words each)
    ///
    /// The shift count drops from 127 to 63 because W=64 outer iterations
    /// replace m=128 bit iterations. The XOR-add count stays the same (one per
    /// set bit of A) but each touches 3 words instead of 2 due to B's growth.
    /// Net effect: fewer iterations, better pipeline behavior, explicit reduction
    /// at the end via reduce_128_gf2p128.
    ///
    /// # Parameters
    ///
    /// m = 128, W = 64, t = 2.
    /// B grows to at most t+1 = 3 words during the shifts.
    /// C accumulates as [u64; 4] (degree ≤ 2m-2 = 254), then reduced.
    pub fn mul_2_34(self, rhs: Self) -> Self {
        // B stored as 3 words: the product B*z^k grows beyond 2 words as k increases.
        // Starts as [b.lo, b.hi, 0].
        let mut b = [rhs.0[0], rhs.0[1], 0u64];

        // C accumulates the unreduced product, needs 4 words (degree up to 254).
        let mut c = [0u64; 4];

        // Outer loop: k from 0 to W-1 = 0 to 63.
        // On iteration k, B holds rhs * z^k.
        for k in 0..64usize {
            // Inner loop: j = 0, 1 (t = 2 words in A).
            // If bit k of A[j] is set, add B into C starting at word j.
            // This accounts for the z^(64j) positional weight of A[j].
            for j in 0..2usize {
                if (self.0[j] >> k) & 1 == 1 {
                    c[j] ^= b[0];
                    c[j + 1] ^= b[1];
                    c[j + 2] ^= b[2];
                }
            }

            // B <- B * z (shift B left by 1 bit across 3 words), except after last iteration.
            // Carry bit 63 of word i into bit 0 of word i+1.
            if k != 63 {
                let carry_01 = b[0] >> 63;
                let carry_12 = b[1] >> 63;
                b[0] <<= 1;
                b[1] = (b[1] << 1) | carry_01;
                b[2] = (b[2] << 1) | carry_12;
            }
        }

        // Reduce the 4-word intermediate mod f(z) using the fast word-level reduction.
        let r = crate::reduce::reduce_2_41(c);
        GF2_128::new(r[0], r[1])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_zero_is_identity() {
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd);
        assert_eq!(a.add_2_32(GF2_128::zero()), a);
    }

    #[test]
    fn add_self_is_zero() {
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd);
        assert_eq!(a.add_2_32(a), GF2_128::zero());
    }

    #[test]
    fn add_is_commutative() {
        let a = GF2_128::new(0x1234, 0x5678);
        let b = GF2_128::new(0xabcd, 0xef01);
        assert_eq!(a.add_2_32(b), b.add_2_32(a));
    }

    #[test]
    fn sub_equals_add() {
        let a = GF2_128::new(0x1234, 0x5678);
        let b = GF2_128::new(0xabcd, 0xef01);
        assert_eq!(a.sub(b), a.add_2_32(b));
    }

    #[test]
    fn neg_is_identity() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe);
        assert_eq!(a.neg(), a);
    }

    #[test]
    fn mul_by_one_is_identity() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe);
        assert_eq!(a.mul_2_33(GF2_128::one()), a);
    }

    #[test]
    fn mul_by_zero_is_zero() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe);
        assert_eq!(a.mul_2_33(GF2_128::zero()), GF2_128::zero());
    }

    #[test]
    fn mul_is_commutative() {
        let a = GF2_128::new(0x1234, 0x5678);
        let b = GF2_128::new(0xabcd, 0xef01);
        assert_eq!(a.mul_2_33(b), b.mul_2_33(a));
    }

    #[test]
    fn mul_distributive_over_add() {
        let a = GF2_128::new(0x1234, 0x5678);
        let b = GF2_128::new(0xabcd, 0xef01);
        let c = GF2_128::new(0xdead, 0xbeef);
        assert_eq!(a.mul_2_33(b.add_2_32(c)), a.mul_2_33(b).add_2_32(a.mul_2_33(c)));
    }

    #[test]
    fn mul_known_value() {
        // In GF(2^128), (x) * (x) = x^2
        // i.e. element 0b10 * element 0b10 = element 0b100
        let x = GF2_128::new(0b10, 0);
        let x_squared = GF2_128::new(0b100, 0);
        assert_eq!(x.mul_2_33(x), x_squared);
    }

    #[test]
    fn mul_reduction_happens() {
        // Multiply z^127 * z — this should trigger reduction
        // z^127 * z = z^128 = z^7 + z^2 + z + 1 = 0x87
        let z127 = GF2_128::new(0, 1u64 << 63);
        let z1 = GF2_128::new(0b10, 0);
        assert_eq!(z127.mul_2_33(z1), GF2_128::new(0x87, 0));
    }

    #[test]
    fn square_zero_is_zero() {
        assert_eq!(GF2_128::zero().square_2_39(), GF2_128::zero());
    }

    #[test]
    fn square_one_is_one() {
        assert_eq!(GF2_128::one().square_2_39(), GF2_128::one());
    }

    #[test]
    fn square_consistent_with_mul() {
        let a = GF2_128::new(0x1234, 0x5678);
        assert_eq!(a.square_2_39(), a.mul_2_33(a));
    }

    #[test]
    fn square_z_is_z2() {
        // z^2 is just bit 2 set
        let z = GF2_128::new(0b10, 0);
        let z2 = GF2_128::new(0b100, 0);
        assert_eq!(z.square_2_39(), z2);
    }

    #[test]
    fn mul_comb_agrees_with_mul() {
        let pairs = [
            (GF2_128::new(0x1234, 0), GF2_128::new(0x5678, 0)),
            (
                GF2_128::new(0xdeadbeefcafe1234, 0),
                GF2_128::new(0xabcd1234, 0),
            ),
            (GF2_128::new(0, 1u64 << 63), GF2_128::new(0b10, 0)),
            (
                GF2_128::new(0xffffffffffffffff, 0xffffffffffffffff),
                GF2_128::one(),
            ),
        ];
        for (a, b) in pairs {
            assert_eq!(a.mul_2_34(b), a.mul_2_33(b), "mismatch for {:?} * {:?}", a, b);
        }
    }

    #[test]
    fn mul_comb_reduction_happens() {
        let z127 = GF2_128::new(0, 1u64 << 63);
        let z1 = GF2_128::new(0b10, 0);
        assert_eq!(z127.mul_2_34(z1), GF2_128::new(0x87, 0));
    }

    #[test]
    fn mul_associative() {
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let b = GF2_128::new(0x0102030405060708, 0xfedcba9876543210);
        let c = GF2_128::new(0x1111222233334444, 0x5555666677778888);
        // (a * b) * c == a * (b * c)  — fundamental field axiom
        assert_eq!(a.mul_2_33(b).mul_2_33(c), a.mul_2_33(b.mul_2_33(c)));
        assert_eq!(a.mul_2_34(b).mul_2_34(c), a.mul_2_34(b.mul_2_34(c)));
    }

    #[test]
    fn mul_2_33_both_limbs_nonzero() {
        // Both operands have both limbs non-zero; verified against mul_2_34.
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabad1dea12345678);
        let b = GF2_128::new(0x0102030405060708, 0xfedcba9876543210);
        assert_eq!(a.mul_2_33(b), a.mul_2_34(b));
    }
}
