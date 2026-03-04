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
    pub fn add(self, rhs: Self) -> Self {
        GF2_128::new(self.0[0] ^ rhs.0[0], self.0[1] ^ rhs.0[1])
    }

    /// Subtraction is identical to addition in characteristic 2.
    pub fn sub(self, rhs: Self) -> Self {
        self.add(rhs)
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
    pub fn mul(self, rhs: Self) -> Self {
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
                c = c.add(b);
            }
        }

        c
    }

    /// Algorithm 2.39 — Squaring in GF(2^128).
    /// Expands bits then reduces mod f(z) via Algorithm 2.40.
    pub fn square(self) -> Self {
        // Step 1: expand each word into 128 bits
        // A[0] (low word) expands into the low 128 bits of c
        // A[1] (high word) expands into the high 128 bits of c
        let lo = expand_word(self.0[0]);
        let hi = expand_word(self.0[1]);

        // Assemble into [u64; 4] intermediate
        let c: [u64; 4] = [lo as u64, (lo >> 64) as u64, hi as u64, (hi >> 64) as u64];

        // Step 2: reduce mod f(z)
        let reduced = crate::reduce::reduce_128_gf2p128(c);
        GF2_128::new(reduced[0], reduced[1])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_zero_is_identity() {
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd);
        assert_eq!(a.add(GF2_128::zero()), a);
    }

    #[test]
    fn add_self_is_zero() {
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd);
        assert_eq!(a.add(a), GF2_128::zero());
    }

    #[test]
    fn add_is_commutative() {
        let a = GF2_128::new(0x1234, 0x5678);
        let b = GF2_128::new(0xabcd, 0xef01);
        assert_eq!(a.add(b), b.add(a));
    }

    #[test]
    fn sub_equals_add() {
        let a = GF2_128::new(0x1234, 0x5678);
        let b = GF2_128::new(0xabcd, 0xef01);
        assert_eq!(a.sub(b), a.add(b));
    }

    #[test]
    fn neg_is_identity() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe);
        assert_eq!(a.neg(), a);
    }

    #[test]
    fn mul_by_one_is_identity() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe);
        assert_eq!(a.mul(GF2_128::one()), a);
    }

    #[test]
    fn mul_by_zero_is_zero() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe);
        assert_eq!(a.mul(GF2_128::zero()), GF2_128::zero());
    }

    #[test]
    fn mul_is_commutative() {
        let a = GF2_128::new(0x1234, 0x5678);
        let b = GF2_128::new(0xabcd, 0xef01);
        assert_eq!(a.mul(b), b.mul(a));
    }

    #[test]
    fn mul_distributive_over_add() {
        let a = GF2_128::new(0x1234, 0x5678);
        let b = GF2_128::new(0xabcd, 0xef01);
        let c = GF2_128::new(0xdead, 0xbeef);
        assert_eq!(a.mul(b.add(c)), a.mul(b).add(a.mul(c)));
    }

    #[test]
    fn mul_known_value() {
        // In GF(2^128), (x) * (x) = x^2
        // i.e. element 0b10 * element 0b10 = element 0b100
        let x = GF2_128::new(0b10, 0);
        let x_squared = GF2_128::new(0b100, 0);
        assert_eq!(x.mul(x), x_squared);
    }

    #[test]
    fn mul_reduction_happens() {
        // Multiply z^127 * z — this should trigger reduction
        // z^127 * z = z^128 = z^7 + z^2 + z + 1 = 0x87
        let z127 = GF2_128::new(0, 1u64 << 63);
        let z1 = GF2_128::new(0b10, 0);
        assert_eq!(z127.mul(z1), GF2_128::new(0x87, 0));
    }

    #[test]
    fn square_zero_is_zero() {
        assert_eq!(GF2_128::zero().square(), GF2_128::zero());
    }

    #[test]
    fn square_one_is_one() {
        assert_eq!(GF2_128::one().square(), GF2_128::one());
    }

    #[test]
    fn square_consistent_with_mul() {
        let a = GF2_128::new(0x1234, 0x5678);
        assert_eq!(a.square(), a.mul(a));
    }

    #[test]
    fn square_z_is_z2() {
        // z^2 is just bit 2 set
        let z = GF2_128::new(0b10, 0);
        let z2 = GF2_128::new(0b100, 0);
        assert_eq!(z.square(), z2);
    }
}
