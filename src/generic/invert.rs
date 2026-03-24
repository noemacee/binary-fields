/// Algorithm 2.48 — Inversion in GF(2^degree) via the extended Euclidean algorithm.
///
/// Generic for any field degree and irreducible polynomial.
/// Returns `None` iff `a` is zero.
pub fn invert_2_48<const N: usize>(
    a: &[u64; N],
    poly: &[u64; N],
    degree: usize,
) -> Option<[u64; N]> {
    if a.iter().all(|&w| w == 0) {
        return None;
    }

    let work = 2 * N + 2;

    let mut u = vec![0u64; work];
    u[..N].copy_from_slice(a);

    let mut v = vec![0u64; work];
    v[..N].copy_from_slice(poly);
    v[degree / 64] |= 1u64 << (degree % 64);

    let mut g1 = vec![0u64; work];
    g1[0] = 1;
    let mut g2 = vec![0u64; work];

    let mut tmp = vec![0u64; work];

    while deg_slice(&u) > 0 {
        let du = deg_slice(&u);
        let dv = deg_slice(&v);
        let j = du - dv;

        if j < 0 {
            std::mem::swap(&mut u, &mut v);
            std::mem::swap(&mut g1, &mut g2);
            let j = (-j) as usize;
            poly_shl_into(&mut tmp, &v, j);
            xor_assign(&mut u, &tmp);
            poly_shl_into(&mut tmp, &g2, j);
            xor_assign(&mut g1, &tmp);
        } else {
            let j = j as usize;
            poly_shl_into(&mut tmp, &v, j);
            xor_assign(&mut u, &tmp);
            poly_shl_into(&mut tmp, &g2, j);
            xor_assign(&mut g1, &tmp);
        }
    }

    let mut result = [0u64; N];
    result.copy_from_slice(&g1[..N]);
    Some(result)
}

/// Algorithm 2.49 — Inversion in GF(2^degree) via the binary GCD algorithm.
///
/// Generic for any field degree and irreducible polynomial.
/// Returns `None` iff `a` is zero.
pub fn invert_2_49<const N: usize>(
    a: &[u64; N],
    poly: &[u64; N],
    degree: usize,
) -> Option<[u64; N]> {
    if a.iter().all(|&w| w == 0) {
        return None;
    }

    // u and v need N+1 limbs to hold the full polynomial (including the leading term z^degree).
    let work = N + 1;

    let mut u = vec![0u64; work];
    u[..N].copy_from_slice(a);

    let mut v = vec![0u64; work];
    v[..N].copy_from_slice(poly);
    v[degree / 64] |= 1u64 << (degree % 64);

    // g1, g2 stay at degree < m so N+1 limbs is sufficient (the +1 absorbs the
    // temporary degree-m intermediate when XORing with f before the right-shift).
    let mut g1 = vec![0u64; work];
    g1[0] = 1;
    let mut g2 = vec![0u64; work];

    loop {
        if u[0] == 1 && u[1..].iter().all(|&w| w == 0) {
            let mut result = [0u64; N];
            result.copy_from_slice(&g1[..N]);
            return Some(result);
        }
        if v[0] == 1 && v[1..].iter().all(|&w| w == 0) {
            let mut result = [0u64; N];
            result.copy_from_slice(&g2[..N]);
            return Some(result);
        }

        while u[0] & 1 == 0 {
            shr1(&mut u);
            adjust_g(&mut g1, poly, degree);
        }

        while v[0] & 1 == 0 {
            shr1(&mut v);
            adjust_g(&mut g2, poly, degree);
        }

        if cmp_gt(&u, &v) {
            xor_assign(&mut u, &v);
            xor_assign(&mut g1, &g2);
        } else {
            xor_assign(&mut v, &u);
            xor_assign(&mut g2, &g1);
        }
    }
}

/// Right-shift a limb slice by 1 bit.
fn shr1(a: &mut [u64]) {
    let n = a.len();
    for i in 0..n - 1 {
        a[i] = (a[i] >> 1) | (a[i + 1] << 63);
    }
    a[n - 1] >>= 1;
}

/// Divide g by z in the polynomial ring mod f: if g is even, right-shift;
/// if g is odd, XOR with f (adding 0 mod f makes it even) then right-shift.
fn adjust_g(g: &mut [u64], poly: &[u64], degree: usize) {
    if g[0] & 1 == 0 {
        shr1(g);
    } else {
        let n = poly.len();
        for i in 0..n {
            g[i] ^= poly[i];
        }
        g[degree / 64] ^= 1u64 << (degree % 64);
        shr1(g);
    }
}

/// Returns true if `a > b` (lexicographic on limbs, most-significant first).
fn cmp_gt(a: &[u64], b: &[u64]) -> bool {
    for i in (0..a.len()).rev() {
        if a[i] != b[i] {
            return a[i] > b[i];
        }
    }
    false
}

fn xor_assign(a: &mut [u64], b: &[u64]) {
    for (x, y) in a.iter_mut().zip(b.iter()) {
        *x ^= y;
    }
}

/// Returns the degree of a polynomial as a slice of u64 limbs.
/// Returns -1 if the polynomial is zero.
fn deg_slice(a: &[u64]) -> i64 {
    for i in (0..a.len()).rev() {
        if a[i] != 0 {
            return i as i64 * 64 + 63 - a[i].leading_zeros() as i64;
        }
    }
    -1
}

// Left shift a polynomial by a given number of bits.
/// Shifts the polynomial `src` by `shift` bits and stores the result in `dst`.
/// `dst` must be large enough to hold the result.
/// `shift` must be less than the degree of the polynomial.
fn poly_shl_into(dst: &mut [u64], src: &[u64], shift: usize) {
    for w in dst.iter_mut() {
        *w = 0;
    }
    if shift == 0 {
        let len = src.len().min(dst.len());
        dst[..len].copy_from_slice(&src[..len]);
        return;
    }
    let w_shift = shift / 64;
    let b_shift = shift % 64;
    for i in 0..src.len() {
        let d = i + w_shift;
        if d < dst.len() {
            dst[d] ^= src[i] << b_shift;
        }
        if b_shift > 0 {
            let d1 = d + 1;
            if d1 < dst.len() {
                dst[d1] ^= src[i] >> (64 - b_shift);
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::generic::arithmetic::mul_2_33;

    const GCM_POLY: [u64; 2] = [0x87, 0];
    const DEGREE_128: usize = 128;

    #[test]
    fn invert_2_48_zero_is_none() {
        assert!(invert_2_48::<2>(&[0, 0], &GCM_POLY, DEGREE_128).is_none());
    }

    #[test]
    fn invert_2_48_one_is_one() {
        assert_eq!(
            invert_2_48::<2>(&[1, 0], &GCM_POLY, DEGREE_128),
            Some([1u64, 0])
        );
    }

    #[test]
    fn invert_2_48_roundtrip() {
        let a: [u64; 2] = [0xdeadbeefcafe1234, 0xabcd1234];
        let inv = invert_2_48(&a, &GCM_POLY, DEGREE_128).unwrap();
        assert_eq!(mul_2_33(&a, &inv, &GCM_POLY, DEGREE_128), [1, 0]);
    }

    #[test]
    fn invert_2_48_known_vector() {
        let a: [u64; 2] = [0x1234, 0];
        assert_eq!(
            invert_2_48(&a, &GCM_POLY, DEGREE_128),
            Some([0xed291313806d0de5, 0x1a1bd097720be28e]),
        );
    }

    #[test]
    fn invert_2_49_zero_is_none() {
        assert!(invert_2_49::<2>(&[0, 0], &GCM_POLY, DEGREE_128).is_none());
    }

    #[test]
    fn invert_2_49_one_is_one() {
        assert_eq!(
            invert_2_49::<2>(&[1, 0], &GCM_POLY, DEGREE_128),
            Some([1u64, 0])
        );
    }

    #[test]
    fn invert_2_49_roundtrip() {
        let a: [u64; 2] = [0xdeadbeefcafe1234, 0xabcd1234];
        let inv = invert_2_49(&a, &GCM_POLY, DEGREE_128).unwrap();
        assert_eq!(mul_2_33(&a, &inv, &GCM_POLY, DEGREE_128), [1, 0]);
    }

    #[test]
    fn invert_2_49_agrees_with_2_48() {
        let cases: &[[u64; 2]] = &[
            [0x1234, 0x5678],
            [0xdeadbeefcafe1234, 0xabcd1234],
            [0xffffffffffffffff, 0xffffffffffffffff],
            [1, 0],
        ];
        for a in cases {
            assert_eq!(
                invert_2_49(a, &GCM_POLY, DEGREE_128),
                invert_2_48(a, &GCM_POLY, DEGREE_128)
            );
        }
    }
}
