/// Algorithm 2.40 — Modular reduction one bit at a time.
/// INPUT: c(z) of degree at most 2m-2 = 254, stored as [u64; 4].
/// OUTPUT: c(z) mod f(z), stored as [u64; 2].
///
/// f(z) = z^128 + z^7 + z^2 + z + 1
/// r(z) = z^7 + z^2 + z + 1 = 0x87
/// m = 128, W = 64, t = 2
///
/// Precomputation: u_k(z) = z^k * r(z) for k in 0..W-1
/// For each bit i from 2m-2=254 downto m=128:
///   if c_i = 1:
///     j = floor((i - m) / W)  — which word of C to add into
///     k = (i - m) - W*j       — which precomputed u_k to use
///     add u_k(z) to C{j}

/// r(z) = z^7 + z^2 + z + 1 = 0x87
const R: u64 = 0x87;

/// Precompute u_k(z) = z^k * r(z) for k in 0..64.
/// Each u_k fits in at most 2 words since deg(r) = 7 and k <= 63,
/// so deg(u_k) <= 70 which fits in 2 x 64-bit words.
fn precompute_uk() -> [(u64, u64); 64] {
    let mut uk = [(0u64, 0u64); 64];
    for k in 0..64 {
        // z^k * r(z): shift R left by k bits across two words
        let lo = R << k;
        let hi = if k == 0 { 0 } else { R >> (64 - k) };
        uk[k] = (lo, hi);
    }
    uk
}

/// Algorithm 2.40 — Modular reduction for GF(2^128).
/// Takes a [u64; 4] intermediate and reduces it mod f(z).
pub fn reduce_2_40(c: [u64; 4]) -> [u64; 2] {
    let mut c = c;
    let uk = precompute_uk();

    // Process bits from 254 downto 128
    // i is the bit position in c
    for i in (128..=254usize).rev() {
        // Extract bit i from c
        let word = i / 64;
        let bit = i % 64;
        if (c[word] >> bit) & 1 == 1 {
            // j = floor((i - 128) / 64) — word offset to add into
            // k = (i - 128) - 64*j     — which u_k to use
            let j = (i - 128) / 64;
            let k = (i - 128) % 64;

            // Add u_k(z) into C{j} — i.e. starting at word j
            let (lo, hi) = uk[k];
            c[j] ^= lo;
            if j + 1 < 4 {
                c[j + 1] ^= hi;
            }
        }
    }

    [c[0], c[1]]
}

/// Fast reduction modulo f(z) = z^128 + z^7 + z^2 + z + 1, with W = 64.
///
/// # Derivation (following HMV §2.3.5, word-at-a-time technique)
///
/// This is a specialization of Algorithm 2.40 to our concrete polynomial and
/// word size, analogous to Algorithms 2.41–2.45 in HMV but for W=64 instead
/// of W=32.
///
/// ## Setup
///
/// Parameters: m=128, W=64, t=2.
/// Input C = [C[0], C[1], C[2], C[3]], degree ≤ 254.
/// r(z) = z^7 + z^2 + z + 1  (the non-leading terms of f).
///
/// The key identity is:
///   z^128 ≡ z^7 + z^2 + z + 1  (mod f(z))
///
/// So for any k ≥ 0:
///   z^(128+k) ≡ z^k · r(z) = z^(k+7) + z^(k+2) + z^(k+1) + z^k
///
/// ## Reducing C[3] (bits 192..255)
///
/// C[3] represents bits z^192 .. z^255, i.e. z^(128 + 64+k) for k=0..63.
///
///   z^(192+k) ≡ z^(k+71) + z^(k+66) + z^(k+65) + z^(k+64)
///
/// These four families of bits land at:
///
///   z^(k+64) : bits 64..127  → entirely in C[1]  (shift C[3] left by 0)
///   z^(k+65) : bits 65..128  → bits 65..127 in C[1], bit 128 overflows to C[2]
///   z^(k+66) : bits 66..129  → bits 66..127 in C[1], bits 128..129 overflow to C[2]
///   z^(k+71) : bits 71..134  → bits 71..127 in C[1], bits 128..134 overflow to C[2]
///
/// Translating to word-level XOR/shift operations (relative to C[1] base = bit 64):
///
///   C[1] ^= C[3]        // z^(k+64): no shift needed
///   C[1] ^= C[3] << 1   // z^(k+65): shift left 1 within C[1]
///   C[1] ^= C[3] << 2   // z^(k+66): shift left 2 within C[1]
///   C[1] ^= C[3] << 7   // z^(k+71): shift left 7 within C[1]
///
/// Overflow into C[2] (bits that fall at position ≥ 128):
///
///   C[2] ^= C[3] >> 63  // overflow from << 1
///   C[2] ^= C[3] >> 62  // overflow from << 2
///   C[2] ^= C[3] >> 57  // overflow from << 7
///   (no overflow from << 0 since max bit is 64+63=127)
///
/// ## Reducing C[2] (bits 128..191), now updated
///
/// C[2] represents z^(128+k) for k=0..63.
///
///   z^(128+k) ≡ z^(k+7) + z^(k+2) + z^(k+1) + z^k
///
/// These four families land at:
///
///   z^(k+0) : bits 0..63    → entirely in C[0]
///   z^(k+1) : bits 1..64    → bits 1..63 in C[0], bit 64 overflows to C[1]
///   z^(k+2) : bits 2..65    → bits 2..63 in C[0], bits 64..65 overflow to C[1]
///   z^(k+7) : bits 7..70    → bits 7..63 in C[0], bits 64..70 overflow to C[1]
///
/// Word-level operations:
///
///   C[0] ^= C[2]        // z^(k+0): no shift
///   C[0] ^= C[2] << 1   // z^(k+1)
///   C[0] ^= C[2] << 2   // z^(k+2)
///   C[0] ^= C[2] << 7   // z^(k+7)
///
/// Overflow into C[1]:
///
///   C[1] ^= C[2] >> 63  // overflow from << 1
///   C[1] ^= C[2] >> 62  // overflow from << 2
///   C[1] ^= C[2] >> 57  // overflow from << 7
///
/// ## Result
///
/// After both steps, C[2] and C[3] are fully absorbed. The result is [C[1], C[0]].
/// No masking step is needed because m=128 falls exactly on a word boundary (128 = 2×64),
/// so there are no partial high-word bits to clear — unlike the NIST W=32 cases which
/// require an explicit mask in their final step.
///
/// Total cost: 14 XOR/shift operations, versus 127 loop iterations in Algorithm 2.40.
pub fn reduce_2_41(c: [u64; 4]) -> [u64; 2] {
    let mut c = c;

    // --- Step 1: Reduce C[3] into C[1] and C[2] ---
    let t = c[3];
    c[1] ^= t;
    c[1] ^= t << 1;
    c[1] ^= t << 2;
    c[1] ^= t << 7;
    c[2] ^= t >> 63;
    c[2] ^= t >> 62;
    c[2] ^= t >> 57;

    // --- Step 2: Reduce C[2] (now updated) into C[0] and C[1] ---
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduce_zero() {
        assert_eq!(reduce_2_40([0, 0, 0, 0]), [0, 0]);
    }

    #[test]
    fn reduce_one() {
        // 1 is already reduced
        assert_eq!(reduce_2_40([1, 0, 0, 0]), [1, 0]);
    }

    #[test]
    fn reduce_z128() {
        // z^128 = r(z) = 0x87
        // z^128 is bit 128 = bit 0 of word 2
        assert_eq!(reduce_2_40([0, 0, 1, 0]), [0x87, 0]);
    }

    #[test]
    fn reduce_z129() {
        // z^129 = z * r(z) = 0x87 << 1 = 0x10e
        assert_eq!(reduce_2_40([0, 0, 0b10, 0]), [0x10e, 0]);
    }
    #[test]
    fn fast_reduce_zero() {
        assert_eq!(reduce_2_41([0, 0, 0, 0]), [0, 0]);
    }

    #[test]
    fn fast_reduce_one() {
        assert_eq!(reduce_2_41([1, 0, 0, 0]), [1, 0]);
    }

    #[test]
    fn fast_reduce_z128() {
        // z^128 ≡ z^7 + z^2 + z + 1 = 0x87
        assert_eq!(reduce_2_41([0, 0, 1, 0]), [0x87, 0]);
    }

    #[test]
    fn fast_reduce_z129() {
        // z^129 ≡ z^8 + z^3 + z^2 + z = 0x10e
        assert_eq!(reduce_2_41([0, 0, 0b10, 0]), [0x10e, 0]);
    }

    #[test]
    fn fast_reduce_z192() {
        // z^192 = z^(128+64) ≡ z^71 + z^66 + z^65 + z^64
        // bits 64..71 relative to result start → C[1] = (1<<7)|(1<<2)|(1<<1)|(1<<0) = 0x87
        assert_eq!(reduce_2_41([0, 0, 0, 1]), [0, 0x87]);
    }

    #[test]
    fn fast_reduce_z254() {
        // z^254 = z^(128+126): highest valid degree
        // z^254 ≡ z^(126+7) + z^(126+2) + z^(126+1) + z^126
        //       = z^133 + z^128 + z^127 + z^126
        // z^133 = z^(128+5) ≡ z^12 + z^7 + z^6 + z^5
        // z^128 ≡ z^7 + z^2 + z + 1
        // Final: z^127 + z^126 + z^12 + (z^7+z^7) + z^6 + z^5 + z^2 + z + 1
        //      = z^127 + z^126 + z^12 + z^6 + z^5 + z^2 + z + 1
        let result = reduce_2_41([0, 0, 0, 1 << 62]); // C[3] bit 62 = z^254
        let expected_slow = reduce_2_40([0, 0, 0, 1 << 62]);
        assert_eq!(result, expected_slow);
    }

    /// Cross-check: fast and general reduction must agree on all valid inputs.
    /// Valid means degree ≤ 254, i.e. C[3] bit 63 = 0 (C[3] < 1<<63).
    #[test]
    fn fast_matches_general_z128() {
        let input = [0, 0, 1, 0];
        assert_eq!(reduce_2_41(input), reduce_2_40(input));
    }

    #[test]
    fn fast_matches_general_z192() {
        let input = [0, 0, 0, 1];
        assert_eq!(reduce_2_41(input), reduce_2_40(input));
    }

    #[test]
    fn fast_matches_general_z254() {
        // Highest valid input: C[3] bit 62
        let input = [0, 0, 0, 1 << 62];
        assert_eq!(reduce_2_41(input), reduce_2_40(input));
    }

    #[test]
    fn fast_matches_general_all_ones_valid() {
        // Valid all-ones: C[3] has bit 63 cleared (max valid C[3] = u64::MAX >> 1)
        let input = [u64::MAX, u64::MAX, u64::MAX, u64::MAX >> 1];
        assert_eq!(reduce_2_41(input), reduce_2_40(input));
    }

    #[test]
    fn fast_matches_general_c3_all_low_bits() {
        // C[3] with all bits set except bit 63
        let input = [0, 0, 0, u64::MAX >> 1];
        assert_eq!(reduce_2_41(input), reduce_2_40(input));
    }

    #[test]
    fn fast_matches_general_c2_all_bits() {
        // C[2] fully set is valid (z^128..z^191, all ≤ 254)
        let input = [0, 0, u64::MAX, 0];
        assert_eq!(reduce_2_41(input), reduce_2_40(input));
    }

    #[test]
    fn fast_matches_general_fixed_vectors() {
        // Fixed test vectors — all have C[3] bit 63 = 0 (valid inputs)
        let inputs: &[[u64; 4]] = &[
            [
                0xdeadbeefcafe1234,
                0x0102030405060708,
                0xfedcba9876543210,
                0x1111222233334444,
            ],
            [
                0x0000000000000001,
                0x8000000000000000,
                0x0000000000000001,
                0x0000000000000001,
            ],
            [
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0x7fffffffffffffff,
            ],
            [
                0x5a5a5a5a5a5a5a5a,
                0xa5a5a5a5a5a5a5a5,
                0x5a5a5a5a5a5a5a5a,
                0x255a5a5a5a5a5a5a,
            ],
            // Single bits at boundary positions
            [0, 0, 1 << 63, 0],       // z^191 (highest bit of C[2])
            [0, 0, 0, 1 << 62],       // z^254 (highest valid bit of C[3])
            [0, 0, 1 << 63, 1 << 62], // z^191 + z^254
        ];
        for &input in inputs {
            // Verify C[3] bit 63 is clear (valid input precondition)
            assert_eq!(input[3] >> 63, 0, "test vector has invalid degree-255 term");
            assert_eq!(
                reduce_2_41(input),
                reduce_2_40(input),
                "mismatch on input {:?}",
                input
            );
        }
    }

    // Ground-truth test vectors verified against the galois Python library.
    // These use unreduced products of field elements as inputs, so C[2]/C[3]
    // are genuinely non-zero and the reduction logic is fully exercised.
    //
    // Python:
    //   import galois
    //   GF = galois.GF(2**128, irreducible_poly="x^128 + x^7 + x^2 + x + 1")
    //   # see tests/generate_values.py for full derivation

    #[test]
    fn galois_reduce_z254() {
        // z^127 * z^127 = z^254, which requires reducing from C[3].
        // z^254 ≡ z^127 + z^126 + z^12 + z^6 + z^5 + z^2 + z + 1
        //       = [0x1067, 0xc000000000000000]
        let input = [0x0, 0x0, 0x0, 0x4000000000000000];
        let expected = [0x1067, 0xc000000000000000];
        assert_eq!(reduce_2_40(input), expected);
        assert_eq!(reduce_2_41(input), expected);
    }

    #[test]
    fn galois_reduce_dense() {
        // (2^128 - 1) * (2^128 - 1): both operands are all-ones (every bit set),
        // their unreduced product has C[2] and C[3] both fully populated.
        // Verified via galois: output = [0x555555555555402f, 0x5555555555555555]
        let input = [
            0x5555555555555555,
            0x5555555555555555,
            0x5555555555555555,
            0x5555555555555555,
        ];
        let expected = [0x555555555555402f, 0x5555555555555555];
        assert_eq!(reduce_2_40(input), expected);
        assert_eq!(reduce_2_41(input), expected);
    }
}
