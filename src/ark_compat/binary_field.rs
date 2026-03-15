//! Generic binary extension field `BinaryField<P, N>` and its config trait.
//!
//! `N` is the number of 64-bit limbs (e.g. N=2 for GF(2^128)).
//! `P` is a zero-size config struct that supplies the arithmetic.

use core::marker::PhantomData;

use ark_ff::{AdditiveGroup, Field, LegendreSymbol, SqrtPrecomputation, One, Zero};
use ark_serialize::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, Compress, EmptyFlags, Flags, SerializationError, Valid, Validate,
    buffer_byte_size,
};
use ark_std::{
    fmt::{Debug, Display, Formatter},
    hash::{Hash, Hasher},
    iter::{Product, Sum},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use zeroize::Zeroize;

use crate::ark_compat::Gf2;

// ─── BinaryFieldConfig ───────────────────────────────────────────────────────

/// Configuration for a binary extension field GF(2^DEGREE).
///
/// Implementors supply the DEGREE, identity constants, a primitive generator,
/// and the three performance-critical arithmetic operations.
pub trait BinaryFieldConfig<const N: usize>: Send + Sync + 'static + Sized {
    /// The extension degree over GF(2).
    ///
    /// Defaults to `N * 64` (all limbs fully used).  Override only for
    /// non-round degrees such as GF(2^127) where the top limb is partially
    /// used.  Must satisfy `DEGREE ≤ N * 64` and `DEGREE > (N-1) * 64`.
    const DEGREE: usize = N * 64;

    /// The additive identity [0, 0, ..., 0].
    const ZERO: BinaryField<Self, N>;

    /// The multiplicative identity, i.e. the polynomial 1.
    const ONE: BinaryField<Self, N>;

    /// Compute `a *= b` (carry-less multiplication mod the irreducible polynomial).
    fn mul_assign(a: &mut [u64; N], b: &[u64; N]);

    /// Compute `a *= a` (squaring, typically faster than general multiplication).
    fn square_in_place(a: &mut [u64; N]);

    /// Compute `a^{-1}`.  Returns `None` iff `a` is zero.
    fn inverse(a: &[u64; N]) -> Option<[u64; N]>;
}

// ─── BinaryField struct ───────────────────────────────────────────────────────

/// An element of the binary extension field GF(2^DEGREE) configured by `P`.
///
/// Limb layout: `limbs[0]` holds the low-degree coefficients (bits 0..63),
/// `limbs[N-1]` holds the high-degree coefficients.
pub struct BinaryField<P: BinaryFieldConfig<N>, const N: usize>(
    pub [u64; N],
    pub PhantomData<P>,
);

impl<P: BinaryFieldConfig<N>, const N: usize> BinaryField<P, N> {
    /// Construct directly from a limb array.
    #[inline]
    pub fn new(limbs: [u64; N]) -> Self {
        BinaryField(limbs, PhantomData)
    }
}

// ─── Copy / Clone / Default ──────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Copy for BinaryField<P, N> {}

impl<P: BinaryFieldConfig<N>, const N: usize> Clone for BinaryField<P, N> {
    fn clone(&self) -> Self { *self }
}

impl<P: BinaryFieldConfig<N>, const N: usize> Default for BinaryField<P, N> {
    fn default() -> Self { P::ZERO }
}

// ─── PartialEq / Eq ──────────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> PartialEq for BinaryField<P, N> {
    fn eq(&self, other: &Self) -> bool { self.0 == other.0 }
}

impl<P: BinaryFieldConfig<N>, const N: usize> Eq for BinaryField<P, N> {}

// ─── Hash ────────────────────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Hash for BinaryField<P, N> {
    fn hash<H: Hasher>(&self, state: &mut H) { self.0.hash(state); }
}

// ─── Ord ─────────────────────────────────────────────────────────────────────
// Lexicographic on limbs (most-significant limb first) — arbitrary but consistent.

impl<P: BinaryFieldConfig<N>, const N: usize> PartialOrd for BinaryField<P, N> {
    fn partial_cmp(&self, other: &Self) -> Option<ark_std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<P: BinaryFieldConfig<N>, const N: usize> Ord for BinaryField<P, N> {
    fn cmp(&self, other: &Self) -> ark_std::cmp::Ordering {
        use ark_std::cmp::Ordering;
        for i in (0..N).rev() {
            match self.0[i].cmp(&other.0[i]) {
                Ordering::Equal => continue,
                ord => return ord,
            }
        }
        Ordering::Equal
    }
}

// ─── Debug / Display ─────────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Debug for BinaryField<P, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        write!(f, "BinaryField(")?;
        for (i, w) in self.0.iter().enumerate().rev() {
            if i == N - 1 { write!(f, "{:#018x}", w)?; }
            else { write!(f, "_{:#018x}", w)?; }
        }
        write!(f, ")")
    }
}

impl<P: BinaryFieldConfig<N>, const N: usize> Display for BinaryField<P, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        Debug::fmt(self, f)
    }
}

// ─── Zeroize ─────────────────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Zeroize for BinaryField<P, N> {
    fn zeroize(&mut self) {
        for w in &mut self.0 { *w = 0; }
    }
}

// ─── Zero / One ──────────────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Zero for BinaryField<P, N> {
    fn zero() -> Self { P::ZERO }
    fn is_zero(&self) -> bool { self.0.iter().all(|&w| w == 0) }
}

impl<P: BinaryFieldConfig<N>, const N: usize> One for BinaryField<P, N> {
    fn one() -> Self { P::ONE }
    fn is_one(&self) -> bool { *self == P::ONE }
}

// ─── Neg ─────────────────────────────────────────────────────────────────────
// In characteristic 2, negation is the identity.

impl<P: BinaryFieldConfig<N>, const N: usize> Neg for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn neg(self) -> Self { self }
}

// ─── AddAssign / SubAssign (all variants) ────────────────────────────────────
// Addition = XOR in GF(2^m).  Subtraction is identical.

impl<P: BinaryFieldConfig<N>, const N: usize> AddAssign<&Self> for BinaryField<P, N> {
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        for i in 0..N { self.0[i] ^= other.0[i]; }
    }
}
impl<P: BinaryFieldConfig<N>, const N: usize> SubAssign<&Self> for BinaryField<P, N> {
    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        for i in 0..N { self.0[i] ^= other.0[i]; }
    }
}
impl<P: BinaryFieldConfig<N>, const N: usize> AddAssign<Self> for BinaryField<P, N> {
    #[inline] fn add_assign(&mut self, other: Self) { *self += &other; }
}
impl<P: BinaryFieldConfig<N>, const N: usize> SubAssign<Self> for BinaryField<P, N> {
    #[inline] fn sub_assign(&mut self, other: Self) { *self -= &other; }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> AddAssign<&'a mut Self> for BinaryField<P, N> {
    #[inline] fn add_assign(&mut self, other: &'a mut Self) { *self += &*other; }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> SubAssign<&'a mut Self> for BinaryField<P, N> {
    #[inline] fn sub_assign(&mut self, other: &'a mut Self) { *self -= &*other; }
}

// ─── MulAssign / DivAssign ───────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> MulAssign<&Self> for BinaryField<P, N> {
    #[inline] fn mul_assign(&mut self, other: &Self) { P::mul_assign(&mut self.0, &other.0); }
}
impl<P: BinaryFieldConfig<N>, const N: usize> DivAssign<&Self> for BinaryField<P, N> {
    #[inline]
    fn div_assign(&mut self, other: &Self) {
        let inv = other.inverse().expect("division by zero in BinaryField");
        P::mul_assign(&mut self.0, &inv.0);
    }
}
impl<P: BinaryFieldConfig<N>, const N: usize> MulAssign<Self> for BinaryField<P, N> {
    #[inline] fn mul_assign(&mut self, other: Self) { *self *= &other; }
}
impl<P: BinaryFieldConfig<N>, const N: usize> DivAssign<Self> for BinaryField<P, N> {
    #[inline] fn div_assign(&mut self, other: Self) { *self /= &other; }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> MulAssign<&'a mut Self> for BinaryField<P, N> {
    #[inline] fn mul_assign(&mut self, other: &'a mut Self) { *self *= &*other; }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> DivAssign<&'a mut Self> for BinaryField<P, N> {
    #[inline] fn div_assign(&mut self, other: &'a mut Self) { *self /= &*other; }
}

// ─── Binary ops (reference forms) ────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Add<&Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn add(mut self, other: &Self) -> Self { self += other; self }
}
impl<P: BinaryFieldConfig<N>, const N: usize> Sub<&Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn sub(mut self, other: &Self) -> Self { self -= other; self }
}
impl<P: BinaryFieldConfig<N>, const N: usize> Mul<&Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn mul(mut self, other: &Self) -> Self { self *= other; self }
}
impl<P: BinaryFieldConfig<N>, const N: usize> Div<&Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn div(mut self, other: &Self) -> Self { self /= other; self }
}

// ─── Binary ops (&mut Self) ──────────────────────────────────────────────────

impl<'a, P: BinaryFieldConfig<N>, const N: usize> Add<&'a mut Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn add(self, other: &'a mut Self) -> Self { self + &*other }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> Sub<&'a mut Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn sub(self, other: &'a mut Self) -> Self { self - &*other }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> Mul<&'a mut Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn mul(self, other: &'a mut Self) -> Self { self * &*other }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> Div<&'a mut Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn div(self, other: &'a mut Self) -> Self { self / &*other }
}

// ─── Binary ops (owned) ──────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Add<Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn add(self, other: Self) -> Self { self + &other }
}
impl<P: BinaryFieldConfig<N>, const N: usize> Sub<Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn sub(self, other: Self) -> Self { self - &other }
}
impl<P: BinaryFieldConfig<N>, const N: usize> Mul<Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn mul(self, other: Self) -> Self { self * &other }
}
impl<P: BinaryFieldConfig<N>, const N: usize> Div<Self> for BinaryField<P, N> {
    type Output = Self;
    #[inline] fn div(self, other: Self) -> Self { self / &other }
}

// ─── &Self op &Self ──────────────────────────────────────────────────────────

impl<'a, 'b, P: BinaryFieldConfig<N>, const N: usize> Add<&'b BinaryField<P, N>>
    for &'a BinaryField<P, N>
{
    type Output = BinaryField<P, N>;
    fn add(self, other: &'b BinaryField<P, N>) -> BinaryField<P, N> { *self + other }
}
impl<'a, 'b, P: BinaryFieldConfig<N>, const N: usize> Sub<&'b BinaryField<P, N>>
    for &'a BinaryField<P, N>
{
    type Output = BinaryField<P, N>;
    fn sub(self, other: &'b BinaryField<P, N>) -> BinaryField<P, N> { *self - other }
}
impl<'a, 'b, P: BinaryFieldConfig<N>, const N: usize> Mul<&'b BinaryField<P, N>>
    for &'a BinaryField<P, N>
{
    type Output = BinaryField<P, N>;
    fn mul(self, other: &'b BinaryField<P, N>) -> BinaryField<P, N> { *self * other }
}
impl<'a, 'b, P: BinaryFieldConfig<N>, const N: usize> Div<&'b BinaryField<P, N>>
    for &'a BinaryField<P, N>
{
    type Output = BinaryField<P, N>;
    fn div(self, other: &'b BinaryField<P, N>) -> BinaryField<P, N> { *self / other }
}

// ─── Sum / Product ───────────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Sum<Self> for BinaryField<P, N> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |a, b| a + b)
    }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> Sum<&'a Self> for BinaryField<P, N> {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |a, b| a + b)
    }
}
impl<P: BinaryFieldConfig<N>, const N: usize> Product<Self> for BinaryField<P, N> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |a, b| a * b)
    }
}
impl<'a, P: BinaryFieldConfig<N>, const N: usize> Product<&'a Self> for BinaryField<P, N> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |a, b| a * b)
    }
}

// ─── From<integer> ───────────────────────────────────────────────────────────
// `From<n>` is the ring homomorphism Z → GF(2^m): every integer maps to n mod 2,
// i.e. the constant polynomial 0 or 1.  Use `BinaryField::new(limbs)` directly
// to construct non-trivial polynomial elements from their limb representation.

macro_rules! impl_from_uint {
    ($uint:ty) => {
        impl<P: BinaryFieldConfig<N>, const N: usize> From<$uint> for BinaryField<P, N> {
            fn from(n: $uint) -> Self {
                // Characteristic-2: n mod 2 = low bit.
                let mut limbs = [0u64; N];
                limbs[0] = (n & 1) as u64;
                BinaryField::new(limbs)
            }
        }
    };
}

impl_from_uint!(u8);
impl_from_uint!(u16);
impl_from_uint!(u32);
impl_from_uint!(u64);
impl_from_uint!(u128);

impl<P: BinaryFieldConfig<N>, const N: usize> From<bool> for BinaryField<P, N> {
    fn from(b: bool) -> Self {
        let mut limbs = [0u64; N];
        limbs[0] = b as u64;
        BinaryField::new(limbs)
    }
}

// Signed: in characteristic 2, -n ≡ n mod 2 (same parity as n).
macro_rules! impl_from_int {
    ($int:ty) => {
        impl<P: BinaryFieldConfig<N>, const N: usize> From<$int> for BinaryField<P, N> {
            fn from(n: $int) -> Self {
                let mut limbs = [0u64; N];
                limbs[0] = (n & 1) as u64;
                BinaryField::new(limbs)
            }
        }
    };
}
impl_from_int!(i8);
impl_from_int!(i16);
impl_from_int!(i32);
impl_from_int!(i64);
impl_from_int!(i128);

// ─── UniformRand ─────────────────────────────────────────────────────────────
// For GF(2^DEGREE), every N-limb value with the top (N*64 - DEGREE) bits zeroed
// is a valid field element.

impl<P: BinaryFieldConfig<N>, const N: usize>
    ark_std::rand::distributions::Distribution<BinaryField<P, N>>
    for ark_std::rand::distributions::Standard
{
    fn sample<R: ark_std::rand::Rng + ?Sized>(&self, rng: &mut R) -> BinaryField<P, N> {
        let mut limbs = [0u64; N];
        for w in &mut limbs { *w = rng.next_u64(); }
        // Mask the top limb if DEGREE is not a multiple of 64.
        // For Gf128Config, DEGREE = 128 = 2*64, so no masking needed.
        let top_bits = P::DEGREE % 64;
        if top_bits != 0 {
            let top_idx = P::DEGREE / 64;
            if top_idx < N {
                limbs[top_idx] &= (1u64 << top_bits) - 1;
                for i in (top_idx + 1)..N { limbs[i] = 0; }
            }
        }
        BinaryField::new(limbs)
    }
}

// ─── Serialization ───────────────────────────────────────────────────────────
// Format: ceil(DEGREE/8) bytes, limbs in little-endian order.

impl<P: BinaryFieldConfig<N>, const N: usize> Valid for BinaryField<P, N> {
    fn check(&self) -> Result<(), SerializationError> { Ok(()) }
}

impl<P: BinaryFieldConfig<N>, const N: usize> CanonicalSerializeWithFlags for BinaryField<P, N> {
    fn serialize_with_flags<W: ark_std::io::Write, F: Flags>(
        &self,
        mut writer: W,
        flags: F,
    ) -> Result<(), SerializationError> {
        if F::BIT_SIZE > 8 { return Err(SerializationError::NotEnoughSpace); }
        let output_byte_size = buffer_byte_size(P::DEGREE + F::BIT_SIZE);
        let mut bytes = ark_std::vec::Vec::with_capacity(output_byte_size);
        bytes.resize(output_byte_size, 0u8);
        for i in 0..N {
            let start = i * 8;
            let end = (start + 8).min(output_byte_size);
            bytes[start..end].copy_from_slice(&self.0[i].to_le_bytes()[..(end - start)]);
        }
        bytes[output_byte_size - 1] |= flags.u8_bitmask();
        writer.write_all(&bytes)
            .map_err(SerializationError::IoError)
    }

    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        buffer_byte_size(P::DEGREE + F::BIT_SIZE)
    }
}

impl<P: BinaryFieldConfig<N>, const N: usize> CanonicalSerialize for BinaryField<P, N> {
    fn serialize_with_mode<W: ark_std::io::Write>(
        &self, writer: W, _: Compress,
    ) -> Result<(), SerializationError> {
        self.serialize_with_flags(writer, EmptyFlags)
    }

    fn serialized_size(&self, _: Compress) -> usize {
        self.serialized_size_with_flags::<EmptyFlags>()
    }
}

impl<P: BinaryFieldConfig<N>, const N: usize> CanonicalDeserializeWithFlags for BinaryField<P, N> {
    fn deserialize_with_flags<R: ark_std::io::Read, F: Flags>(
        mut reader: R,
    ) -> Result<(Self, F), SerializationError> {
        if F::BIT_SIZE > 8 { return Err(SerializationError::NotEnoughSpace); }
        let output_byte_size = Self::zero().serialized_size_with_flags::<F>();
        let mut bytes = ark_std::vec::Vec::with_capacity(output_byte_size);
        bytes.resize(output_byte_size, 0u8);
        reader.read_exact(&mut bytes)
            .map_err(SerializationError::IoError)?;
        let flag = F::from_u8_remove_flags(&mut bytes[output_byte_size - 1])
            .ok_or(SerializationError::UnexpectedFlags)?;
        let mut limbs = [0u64; N];
        for i in 0..N {
            let start = i * 8;
            let end = start + 8;
            if end <= output_byte_size {
                let mut buf = [0u8; 8];
                buf.copy_from_slice(&bytes[start..end]);
                limbs[i] = u64::from_le_bytes(buf);
            } else if start < output_byte_size {
                let mut buf = [0u8; 8];
                buf[..output_byte_size - start].copy_from_slice(&bytes[start..output_byte_size]);
                limbs[i] = u64::from_le_bytes(buf);
            }
        }
        Ok((BinaryField::new(limbs), flag))
    }
}

impl<P: BinaryFieldConfig<N>, const N: usize> CanonicalDeserialize for BinaryField<P, N> {
    fn deserialize_with_mode<R: ark_std::io::Read>(
        reader: R, _: Compress, _: Validate,
    ) -> Result<Self, SerializationError> {
        Self::deserialize_with_flags::<R, EmptyFlags>(reader).map(|(r, _)| r)
    }
}

// ─── AdditiveGroup ───────────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> AdditiveGroup for BinaryField<P, N> {
    type Scalar = Self;
    const ZERO: Self = P::ZERO;

    #[inline]
    fn double_in_place(&mut self) -> &mut Self {
        // x + x = 0 in characteristic 2
        for w in &mut self.0 { *w = 0; }
        self
    }

    #[inline]
    fn neg_in_place(&mut self) -> &mut Self {
        // -x = x in characteristic 2
        self
    }
}

// ─── Field ───────────────────────────────────────────────────────────────────

impl<P: BinaryFieldConfig<N>, const N: usize> Field for BinaryField<P, N> {
    type BasePrimeField = Gf2;

    // sqrt in GF(2^m) always exists; we override the method instead of using SQRT_PRECOMP.
    const SQRT_PRECOMP: Option<SqrtPrecomputation<Self>> = None;
    const ONE: Self = P::ONE;

    fn extension_degree() -> u64 { P::DEGREE as u64 }

    fn from_base_prime_field(elem: Gf2) -> Self {
        // Embed a GF(2) element as the constant term (bit 0) of the polynomial.
        let mut limbs = [0u64; N];
        limbs[0] = elem.0.0[0];
        BinaryField::new(limbs)
    }

    fn to_base_prime_field_elements(&self) -> impl Iterator<Item = Gf2> {
        // Yield one Gf2 element per bit of the field element (LSB first).
        let mut bits = ark_std::vec::Vec::with_capacity(P::DEGREE);
        for i in 0..N {
            for bit in 0..64usize {
                if i * 64 + bit < P::DEGREE {
                    bits.push(Gf2::from(((self.0[i] >> bit) & 1) as u64));
                }
            }
        }
        bits.into_iter()
    }

    fn from_base_prime_field_elems(
        elems: impl IntoIterator<Item = Gf2>,
    ) -> Option<Self> {
        let mut limbs = [0u64; N];
        let mut count = 0usize;
        for (idx, elem) in elems.into_iter().enumerate() {
            if idx >= P::DEGREE { return None; } // too many elements
            limbs[idx / 64] |= elem.0.0[0] << (idx % 64);
            count += 1;
        }
        if count != P::DEGREE { return None; }
        Some(BinaryField::new(limbs))
    }

    fn characteristic() -> &'static [u64] {
        Gf2::characteristic()
    }

    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        if F::BIT_SIZE > 8 { return None; }
        let output_byte_size = buffer_byte_size(P::DEGREE + F::BIT_SIZE);
        if bytes.len() < output_byte_size { return None; }

        let mut buf = ark_std::vec::Vec::with_capacity(output_byte_size);
        buf.extend_from_slice(&bytes[..output_byte_size]);

        // Extract flags from the last byte.
        let flags_mask = u8::MAX
            .checked_shl(8u32.saturating_sub(F::BIT_SIZE as u32))
            .unwrap_or(0);
        let flags_byte = buf[output_byte_size - 1] & flags_mask;
        buf[output_byte_size - 1] &= !flags_mask;
        let flag = F::from_u8(flags_byte)?;

        // All bit patterns are valid field elements, so no rejection needed.
        let mut limbs = [0u64; N];
        for i in 0..N {
            let start = i * 8;
            let end = start + 8;
            if end <= output_byte_size {
                let mut b = [0u8; 8];
                b.copy_from_slice(&buf[start..end]);
                limbs[i] = u64::from_le_bytes(b);
            } else if start < output_byte_size {
                let mut b = [0u8; 8];
                b[..output_byte_size - start].copy_from_slice(&buf[start..output_byte_size]);
                limbs[i] = u64::from_le_bytes(b);
            }
        }
        Some((BinaryField::new(limbs), flag))
    }

    fn legendre(&self) -> LegendreSymbol {
        // In GF(2^m), squaring is a bijection: every nonzero element is a square.
        if self.is_zero() { LegendreSymbol::Zero } else { LegendreSymbol::QuadraticResidue }
    }

    fn sqrt(&self) -> Option<Self> {
        // sqrt(a) = a^(2^(m-1))  [since (a^(2^(m-1)))^2 = a^(2^m) = a in GF(2^m)]
        // Computed by repeated squaring P::DEGREE - 1 times.
        let mut result = *self;
        for _ in 0..(P::DEGREE - 1) {
            result.square_in_place();
        }
        Some(result)
    }

    #[inline]
    fn square(&self) -> Self {
        let mut tmp = *self;
        tmp.square_in_place();
        tmp
    }

    #[inline]
    fn square_in_place(&mut self) -> &mut Self {
        P::square_in_place(&mut self.0);
        self
    }

    #[inline]
    fn inverse(&self) -> Option<Self> {
        P::inverse(&self.0).map(BinaryField::new)
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        match P::inverse(&self.0) {
            Some(inv) => { self.0 = inv; Some(self) }
            None => None,
        }
    }

    fn frobenius_map_in_place(&mut self, power: usize) {
        // In GF(2^m), the Frobenius automorphism is x → x^2.
        // Applying it `power` times gives x → x^(2^power).
        for _ in 0..power {
            P::square_in_place(&mut self.0);
        }
    }

    fn mul_by_base_prime_field(&self, elem: &Gf2) -> Self {
        // Multiply by 0 → zero; multiply by 1 → self.
        if elem.is_zero() { Self::zero() } else { *self }
    }
}

#[cfg(test)]
mod tests {
    // Tests live in gf128_config.rs where we have a concrete type to work with.
}
