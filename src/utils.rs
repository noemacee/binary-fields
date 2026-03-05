use crate::GF2_128;

impl GF2_128 {
    /// Construct a field element from a single u64.
    /// The value is placed in the low word, high word is zero.
    /// Useful for small constants like GF2_128::from_u64(0x1234).
    pub fn from_u64(val: u64) -> Self {
        GF2_128::new(val, 0)
    }

    /// Extract the raw [u64; 2] representation.
    /// Returns [lo, hi] where lo = a63..a0, hi = a127..a64.
    pub fn to_words(self) -> [u64; 2] {
        self.0
    }

    /// Exponentiation by square-and-multiply.
    /// Computes a^n in GF(2^128).
    /// a^0 = 1 for any a, including zero.
    pub fn pow(self, n: u128) -> Self {
        if n == 0 {
            return GF2_128::one();
        }

        let mut result = GF2_128::one();
        let mut base = self;
        let mut exp = n;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul_comb(base);
            }
            base = base.square();
            exp >>= 1;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_u64_roundtrip() {
        let a = GF2_128::from_u64(0xdeadbeef);
        assert_eq!(a.to_words()[0], 0xdeadbeef);
        assert_eq!(a.to_words()[1], 0);
    }

    #[test]
    fn to_words_roundtrip() {
        let a = GF2_128::new(0x1234, 0x5678);
        let words = a.to_words();
        assert_eq!(words, [0x1234, 0x5678]);
    }

    #[test]
    fn pow_zero_is_one() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe);
        assert_eq!(a.pow(0), GF2_128::one());
    }

    #[test]
    fn pow_one_is_self() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe);
        assert_eq!(a.pow(1), a);
    }

    #[test]
    fn pow_two_is_mul_self() {
        let a = GF2_128::new(0x1234, 0x5678);
        assert_eq!(a.pow(2), a.mul(a));
    }

    #[test]
    fn pow_fermat() {
        // Fermat's little theorem: a^(2^128 - 1) = 1 for any nonzero a
        // This is the strongest correctness check for the whole field
        let a = GF2_128::new(0xdeadbeefcafe1234, 0xabcd1234);
        let order = u128::MAX; // 2^128 - 1
        assert_eq!(a.pow(order), GF2_128::one());
    }

    #[test]
    fn pow_inverse_consistent() {
        // a^(2^128 - 2) should equal a^{-1}
        let a = GF2_128::new(0x1234, 0x5678);
        let exp = u128::MAX - 1; // 2^128 - 2
        assert_eq!(a.pow(exp), a.invert().unwrap());
    }
}
