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
pub fn reduce(c: [u64; 4]) -> [u64; 2] {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduce_zero() {
        assert_eq!(reduce([0, 0, 0, 0]), [0, 0]);
    }

    #[test]
    fn reduce_one() {
        // 1 is already reduced
        assert_eq!(reduce([1, 0, 0, 0]), [1, 0]);
    }

    #[test]
    fn reduce_z128() {
        // z^128 = r(z) = 0x87
        // z^128 is bit 128 = bit 0 of word 2
        assert_eq!(reduce([0, 0, 1, 0]), [0x87, 0]);
    }

    #[test]
    fn reduce_z129() {
        // z^129 = z * r(z) = 0x87 << 1 = 0x10e
        assert_eq!(reduce([0, 0, 0b10, 0]), [0x10e, 0]);
    }
}
