/// An element of GF(2^128), the binary field with 2^128 elements.
///
/// Representation follows Figure 2.4 of Hankerson-Menezes-Vanstone:
///   - m = 128, W = 64, t = 2, s = 0 (no unused bits)
///   - A[0] holds coefficients a63..a0   (low word)
///   - A[1] holds coefficients a127..a64 (high word)
///   - The rightmost bit of A[0] is a0 (constant term)
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct GF2_128(pub [u64; 2]);

impl GF2_128 {
    /// Construct an element directly from its two 64-bit words.
    /// `lo` holds a63..a0, `hi` holds a127..a64.
    pub fn new(lo: u64, hi: u64) -> Self {
        GF2_128([lo, hi])
    }

    /// The additive identity: the zero polynomial.
    pub fn zero() -> Self {
        GF2_128([0, 0])
    }

    /// The multiplicative identity: the polynomial 1.
    pub fn one() -> Self {
        GF2_128([1, 0])
    }

    /// Returns true if this element is zero.
    pub fn is_zero(self) -> bool {
        self.0[0] == 0 && self.0[1] == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_is_zero() {
        assert!(GF2_128::zero().is_zero());
    }

    #[test]
    fn one_is_not_zero() {
        assert!(!GF2_128::one().is_zero());
    }

    #[test]
    fn new_roundtrip() {
        let a = GF2_128::new(0xdeadbeef, 0xcafe1234);
        assert_eq!(a.0[0], 0xdeadbeef);
        assert_eq!(a.0[1], 0xcafe1234);
    }
}
