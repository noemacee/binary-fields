use crate::GF2_128;

/// The irreducible polynomial f(z) = z^128 + z^7 + z^2 + z + 1 : The one used in Galois Counter Mode (GCM)
/// Since f(z) = 0 in the field, z^128 = z^7 + z^2 + z + 1.
/// r(z) = z^7 + z^2 + z + 1 = 0b10000111 = 0x87
/// is what we XOR with when the high bit (z^128) spills over.
const R: u64 = 0x87;

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
}
