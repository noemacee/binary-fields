use hekate_math::{Block8, Block16, Block32, Block64, Block128, Flat, HardwareField, TowerField};

// ── Tower squaring — 2 sub-squarings + 1 mul-by-τ, no cross term ─────────────
//
// In char 2, for a tower field element (a + b·x) where x² = x + τ:
//   (a + b·x)² = a² + b²·x² = a² + b²·(x + τ) = (a² + b²·τ) + b²·x
//
// Cost: 2 recursive squarings + 1 mul-by-τ (constant mul at next level down),
// vs Karatsuba general mul which costs 3 recursive muls. Savings compound across
// the 4 tower levels: 2^4=16 base muls for squaring vs 3^4=81 for general mul.
//
// Block8 is the base: no public split(), so we fall back to hekate's mul table.
// At Block16 and above, we exploit the tower structure.

// Precomputed square table for GF(2^8) under hekate's tower basis.
// Computed by squaring all 256 elements via Block8 * Block8 at compile time
// — but Block8::mul is not const, so we use a runtime-initialized approach
// and just call a * a inline. The tower savings start at Block16.

#[inline(always)]
fn square_8(a: Block8) -> Block8 { a * a }

#[inline(always)]
fn square_16(a: Block16) -> Block16 {
    let (lo, hi) = a.split();
    let lo2 = square_8(lo);
    let hi2 = square_8(hi);
    // lo2 + hi2 * τ,  hi2
    Block16::new(lo2 + hi2 * Block8::EXTENSION_TAU, hi2)
}

#[inline(always)]
fn square_32(a: Block32) -> Block32 {
    let (lo, hi) = a.split();
    let lo2 = square_16(lo);
    let hi2 = square_16(hi);
    Block32::new(lo2 + hi2 * Block16::TAU, hi2)
}

#[inline(always)]
fn square_64(a: Block64) -> Block64 {
    let (lo, hi) = a.split();
    let lo2 = square_32(lo);
    let hi2 = square_32(hi);
    Block64::new(lo2 + hi2 * Block32::TAU, hi2)
}

#[inline(always)]
pub(crate) fn square_128(a: Block128) -> Block128 {
    let (lo, hi) = a.split();
    let lo2 = square_64(lo);
    let hi2 = square_64(hi);
    Block128::new(lo2 + hi2 * Block64::TAU, hi2)
}
use crate::ark::Gf2;
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

/// Wrapper around hekate's `Block128` that implements `ark_ff::Field`.
///
/// All arithmetic delegates to hekate's optimised implementation.
/// `From<integer>` uses the ring homomorphism (n mod 2), not the
/// bit-pattern embedding — use `Block128Ark(Block128(i as u128))` directly
/// when distinct evaluation points are needed.
#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Block128Ark(pub Block128);

impl Debug for Block128Ark {
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        write!(f, "Block128Ark({:#034x})", self.0.0)
    }
}
impl Display for Block128Ark {
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        Debug::fmt(self, f)
    }
}

impl Hash for Block128Ark {
    fn hash<H: Hasher>(&self, state: &mut H) { self.0.0.hash(state); }
}

impl PartialOrd for Block128Ark {
    fn partial_cmp(&self, other: &Self) -> Option<ark_std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Block128Ark {
    fn cmp(&self, other: &Self) -> ark_std::cmp::Ordering { self.0.0.cmp(&other.0.0) }
}

impl Zeroize for Block128Ark {
    fn zeroize(&mut self) { self.0.0 = 0; }
}

impl Zero for Block128Ark {
    fn zero() -> Self { Self(Block128::ZERO) }
    fn is_zero(&self) -> bool { self.0.0 == 0 }
}
impl One for Block128Ark {
    fn one() -> Self { Self(Block128::ONE) }
    fn is_one(&self) -> bool { self.0.0 == 1 }
}

impl Neg for Block128Ark {
    type Output = Self;
    #[inline] fn neg(self) -> Self { self }
}

// ── AddAssign / SubAssign ────────────────────────────────────────────────────

impl AddAssign<&Self> for Block128Ark {
    #[inline] fn add_assign(&mut self, rhs: &Self) { self.0 += rhs.0; }
}
impl SubAssign<&Self> for Block128Ark {
    #[inline] fn sub_assign(&mut self, rhs: &Self) { self.0 -= rhs.0; }
}
impl AddAssign<Self> for Block128Ark {
    #[inline] fn add_assign(&mut self, rhs: Self) { *self += &rhs; }
}
impl SubAssign<Self> for Block128Ark {
    #[inline] fn sub_assign(&mut self, rhs: Self) { *self -= &rhs; }
}
impl<'a> AddAssign<&'a mut Self> for Block128Ark {
    #[inline] fn add_assign(&mut self, rhs: &'a mut Self) { *self += &*rhs; }
}
impl<'a> SubAssign<&'a mut Self> for Block128Ark {
    #[inline] fn sub_assign(&mut self, rhs: &'a mut Self) { *self -= &*rhs; }
}

// ── MulAssign / DivAssign ────────────────────────────────────────────────────

impl MulAssign<&Self> for Block128Ark {
    #[inline] fn mul_assign(&mut self, rhs: &Self) { self.0 *= rhs.0; }
}
impl DivAssign<&Self> for Block128Ark {
    #[inline]
    fn div_assign(&mut self, rhs: &Self) {
        let inv = rhs.inverse().expect("division by zero in Block128Ark");
        self.0 *= inv.0;
    }
}
impl MulAssign<Self> for Block128Ark {
    #[inline] fn mul_assign(&mut self, rhs: Self) { *self *= &rhs; }
}
impl DivAssign<Self> for Block128Ark {
    #[inline] fn div_assign(&mut self, rhs: Self) { *self /= &rhs; }
}
impl<'a> MulAssign<&'a mut Self> for Block128Ark {
    #[inline] fn mul_assign(&mut self, rhs: &'a mut Self) { *self *= &*rhs; }
}
impl<'a> DivAssign<&'a mut Self> for Block128Ark {
    #[inline] fn div_assign(&mut self, rhs: &'a mut Self) { *self /= &*rhs; }
}

// ── Binary ops (&Self) ───────────────────────────────────────────────────────

impl Add<&Self> for Block128Ark { type Output = Self; #[inline] fn add(mut self, rhs: &Self) -> Self { self += rhs; self } }
impl Sub<&Self> for Block128Ark { type Output = Self; #[inline] fn sub(mut self, rhs: &Self) -> Self { self -= rhs; self } }
impl Mul<&Self> for Block128Ark { type Output = Self; #[inline] fn mul(mut self, rhs: &Self) -> Self { self *= rhs; self } }
impl Div<&Self> for Block128Ark { type Output = Self; #[inline] fn div(mut self, rhs: &Self) -> Self { self /= rhs; self } }

impl<'a> Add<&'a mut Self> for Block128Ark { type Output = Self; #[inline] fn add(self, rhs: &'a mut Self) -> Self { self + &*rhs } }
impl<'a> Sub<&'a mut Self> for Block128Ark { type Output = Self; #[inline] fn sub(self, rhs: &'a mut Self) -> Self { self - &*rhs } }
impl<'a> Mul<&'a mut Self> for Block128Ark { type Output = Self; #[inline] fn mul(self, rhs: &'a mut Self) -> Self { self * &*rhs } }
impl<'a> Div<&'a mut Self> for Block128Ark { type Output = Self; #[inline] fn div(self, rhs: &'a mut Self) -> Self { self / &*rhs } }

impl Add<Self> for Block128Ark { type Output = Self; #[inline] fn add(self, rhs: Self) -> Self { self + &rhs } }
impl Sub<Self> for Block128Ark { type Output = Self; #[inline] fn sub(self, rhs: Self) -> Self { self - &rhs } }
impl Mul<Self> for Block128Ark { type Output = Self; #[inline] fn mul(self, rhs: Self) -> Self { self * &rhs } }
impl Div<Self> for Block128Ark { type Output = Self; #[inline] fn div(self, rhs: Self) -> Self { self / &rhs } }

impl<'a, 'b> Add<&'b Block128Ark> for &'a Block128Ark { type Output = Block128Ark; fn add(self, rhs: &'b Block128Ark) -> Block128Ark { *self + rhs } }
impl<'a, 'b> Sub<&'b Block128Ark> for &'a Block128Ark { type Output = Block128Ark; fn sub(self, rhs: &'b Block128Ark) -> Block128Ark { *self - rhs } }
impl<'a, 'b> Mul<&'b Block128Ark> for &'a Block128Ark { type Output = Block128Ark; fn mul(self, rhs: &'b Block128Ark) -> Block128Ark { *self * rhs } }
impl<'a, 'b> Div<&'b Block128Ark> for &'a Block128Ark { type Output = Block128Ark; fn div(self, rhs: &'b Block128Ark) -> Block128Ark { *self / rhs } }

// ── Sum / Product ────────────────────────────────────────────────────────────

impl Sum<Self> for Block128Ark {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self { iter.fold(Self::zero(), |a, b| a + b) }
}
impl<'a> Sum<&'a Self> for Block128Ark {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self { iter.fold(Self::zero(), |a, b| a + b) }
}
impl Product<Self> for Block128Ark {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self { iter.fold(Self::one(), |a, b| a * b) }
}
impl<'a> Product<&'a Self> for Block128Ark {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self { iter.fold(Self::one(), |a, b| a * b) }
}

// ── From<integer> — ring homomorphism (n mod 2) ──────────────────────────────

macro_rules! impl_from_uint {
    ($t:ty) => {
        impl From<$t> for Block128Ark {
            #[inline] fn from(n: $t) -> Self { Self(Block128((n & 1) as u128)) }
        }
    };
}
macro_rules! impl_from_int {
    ($t:ty) => {
        impl From<$t> for Block128Ark {
            #[inline] fn from(n: $t) -> Self { Self(Block128((n & 1) as u128)) }
        }
    };
}
impl_from_uint!(u8);
impl_from_uint!(u16);
impl_from_uint!(u32);
impl_from_uint!(u64);
impl_from_uint!(u128);
impl_from_int!(i8);
impl_from_int!(i16);
impl_from_int!(i32);
impl_from_int!(i64);
impl_from_int!(i128);
impl From<bool> for Block128Ark {
    #[inline] fn from(b: bool) -> Self { Self(Block128(b as u128)) }
}

// ── UniformRand ──────────────────────────────────────────────────────────────

impl ark_std::rand::distributions::Distribution<Block128Ark>
    for ark_std::rand::distributions::Standard
{
    fn sample<R: ark_std::rand::Rng + ?Sized>(&self, rng: &mut R) -> Block128Ark {
        let lo = rng.next_u64() as u128;
        let hi = rng.next_u64() as u128;
        Block128Ark(Block128(lo | (hi << 64)))
    }
}

// ── Serialization ────────────────────────────────────────────────────────────

impl Valid for Block128Ark {
    fn check(&self) -> Result<(), SerializationError> { Ok(()) }
}

impl CanonicalSerializeWithFlags for Block128Ark {
    fn serialize_with_flags<W: ark_std::io::Write, F: Flags>(
        &self, mut writer: W, flags: F,
    ) -> Result<(), SerializationError> {
        if F::BIT_SIZE > 8 { return Err(SerializationError::NotEnoughSpace); }
        let byte_size = buffer_byte_size(128 + F::BIT_SIZE);
        let mut bytes = ark_std::vec![0u8; byte_size];
        bytes[..16].copy_from_slice(&self.0.0.to_le_bytes());
        bytes[byte_size - 1] |= flags.u8_bitmask();
        writer.write_all(&bytes).map_err(SerializationError::IoError)
    }
    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        buffer_byte_size(128 + F::BIT_SIZE)
    }
}

impl CanonicalSerialize for Block128Ark {
    fn serialize_with_mode<W: ark_std::io::Write>(
        &self, writer: W, _: Compress,
    ) -> Result<(), SerializationError> {
        self.serialize_with_flags(writer, EmptyFlags)
    }
    fn serialized_size(&self, _: Compress) -> usize {
        self.serialized_size_with_flags::<EmptyFlags>()
    }
}

impl CanonicalDeserializeWithFlags for Block128Ark {
    fn deserialize_with_flags<R: ark_std::io::Read, F: Flags>(
        mut reader: R,
    ) -> Result<(Self, F), SerializationError> {
        if F::BIT_SIZE > 8 { return Err(SerializationError::NotEnoughSpace); }
        let byte_size = Self::zero().serialized_size_with_flags::<F>();
        let mut bytes = ark_std::vec![0u8; byte_size];
        reader.read_exact(&mut bytes).map_err(SerializationError::IoError)?;
        let flag = F::from_u8_remove_flags(&mut bytes[byte_size - 1])
            .ok_or(SerializationError::UnexpectedFlags)?;
        let mut buf = [0u8; 16];
        buf.copy_from_slice(&bytes[..16]);
        Ok((Self(Block128(u128::from_le_bytes(buf))), flag))
    }
}

impl CanonicalDeserialize for Block128Ark {
    fn deserialize_with_mode<R: ark_std::io::Read>(
        reader: R, _: Compress, _: Validate,
    ) -> Result<Self, SerializationError> {
        Self::deserialize_with_flags::<R, EmptyFlags>(reader).map(|(r, _)| r)
    }
}

// ── AdditiveGroup ────────────────────────────────────────────────────────────

impl AdditiveGroup for Block128Ark {
    type Scalar = Self;
    const ZERO: Self = Self(Block128(0));
    #[inline]
    fn double_in_place(&mut self) -> &mut Self {
        self.0.0 = 0;
        self
    }
    #[inline]
    fn neg_in_place(&mut self) -> &mut Self { self }
}

// ── Field ────────────────────────────────────────────────────────────────────

impl Field for Block128Ark {
    type BasePrimeField = Gf2;
    const SQRT_PRECOMP: Option<SqrtPrecomputation<Self>> = None;
    const ONE: Self = Self(Block128(1));
    const NEG_ONE: Self = Self(Block128(1)); // -1 = 1 in characteristic 2

    fn extension_degree() -> u64 { 128 }

    fn characteristic() -> &'static [u64] { Gf2::characteristic() }

    fn from_base_prime_field(elem: Gf2) -> Self {
        Self(Block128(elem.0.0[0] as u128))
    }

    fn to_base_prime_field_elements(&self) -> impl Iterator<Item = Gf2> {
        let val = self.0.0;
        (0..128usize).map(move |i| Gf2::from(((val >> i) & 1) as u64))
    }

    fn from_base_prime_field_elems(elems: impl IntoIterator<Item = Gf2>) -> Option<Self> {
        let mut val = 0u128;
        let mut count = 0usize;
        for (i, e) in elems.into_iter().enumerate() {
            if i >= 128 { return None; }
            val |= (e.0.0[0] as u128) << i;
            count += 1;
        }
        if count != 128 { return None; }
        Some(Self(Block128(val)))
    }

    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        if F::BIT_SIZE > 8 { return None; }
        let byte_size = buffer_byte_size(128 + F::BIT_SIZE);
        if bytes.len() < byte_size { return None; }
        let flag = F::from_u8(0)?;
        let mut buf = [0u8; 16];
        buf.copy_from_slice(&bytes[..16]);
        Some((Self(Block128(u128::from_le_bytes(buf))), flag))
    }

    fn legendre(&self) -> LegendreSymbol {
        if self.is_zero() { LegendreSymbol::Zero } else { LegendreSymbol::QuadraticResidue }
    }

    fn sqrt(&self) -> Option<Self> {
        // sqrt(a) = a^(2^127) in GF(2^128)
        let mut r = *self;
        for _ in 0..127 { r.square_in_place(); }
        Some(r)
    }

    #[inline]
    fn square(&self) -> Self { Self(square_128(self.0)) }

    #[inline]
    fn square_in_place(&mut self) -> &mut Self {
        self.0 = square_128(self.0);
        self
    }

    #[inline]
    fn inverse(&self) -> Option<Self> {
        if self.is_zero() { return None; }
        Some(Self(self.0.invert()))
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        if self.is_zero() { return None; }
        self.0 = self.0.invert();
        Some(self)
    }

    fn frobenius_map_in_place(&mut self, power: usize) {
        for _ in 0..power { self.square_in_place(); }
    }

    fn mul_by_base_prime_field(&self, elem: &Gf2) -> Self {
        if elem.is_zero() { Self::zero() } else { *self }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Block128FlatArk — hardware/flat-basis wrapper
// ─────────────────────────────────────────────────────────────────────────────

/// Wrapper around hekate's `Flat<Block128>` that implements `ark_ff::Field`.
///
/// Elements are stored in the polynomial (flat/hardware) basis.
/// Multiplication dispatches to `HardwareField::mul_hardware`, which on
/// aarch64 uses a single `vmull_p64` (PMULL) Karatsuba round instead of
/// the recursive tower Karatsuba used by `Block128Ark`.
///
/// `inverse` and `sqrt` round-trip through the tower basis since no direct
/// flat-basis inversion is implemented yet.
#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Block128FlatArk(pub Flat<Block128>);

impl Block128FlatArk {
    #[inline]
    fn from_flat(f: Flat<Block128>) -> Self { Self(f) }
}

impl Debug for Block128FlatArk {
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        write!(f, "Block128FlatArk({:#034x})", self.0.into_raw().0)
    }
}
impl Display for Block128FlatArk {
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result { Debug::fmt(self, f) }
}

impl Hash for Block128FlatArk {
    fn hash<H: Hasher>(&self, state: &mut H) { self.0.into_raw().0.hash(state); }
}

impl PartialOrd for Block128FlatArk {
    fn partial_cmp(&self, other: &Self) -> Option<ark_std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Block128FlatArk {
    fn cmp(&self, other: &Self) -> ark_std::cmp::Ordering {
        self.0.into_raw().0.cmp(&other.0.into_raw().0)
    }
}

impl Zeroize for Block128FlatArk {
    fn zeroize(&mut self) {
        let mut raw = self.0.into_raw();
        raw.0 = 0;
        self.0 = Flat::from_raw(raw);
    }
}

impl Zero for Block128FlatArk {
    fn zero() -> Self { Self::from_flat(Block128::ZERO.to_hardware()) }
    fn is_zero(&self) -> bool { self.0.into_raw().0 == 0 }
}
impl One for Block128FlatArk {
    fn one() -> Self { Self::from_flat(Block128::ONE.to_hardware()) }
    fn is_one(&self) -> bool { *self == Self::one() }
}

impl Neg for Block128FlatArk {
    type Output = Self;
    #[inline] fn neg(self) -> Self { self }
}

impl AddAssign<&Self> for Block128FlatArk {
    #[inline] fn add_assign(&mut self, rhs: &Self) { self.0 = self.0 + rhs.0; }
}
impl SubAssign<&Self> for Block128FlatArk {
    #[inline] fn sub_assign(&mut self, rhs: &Self) { self.0 = self.0 + rhs.0; }
}
impl AddAssign<Self> for Block128FlatArk {
    #[inline] fn add_assign(&mut self, rhs: Self) { *self += &rhs; }
}
impl SubAssign<Self> for Block128FlatArk {
    #[inline] fn sub_assign(&mut self, rhs: Self) { *self -= &rhs; }
}
impl<'a> AddAssign<&'a mut Self> for Block128FlatArk {
    #[inline] fn add_assign(&mut self, rhs: &'a mut Self) { *self += &*rhs; }
}
impl<'a> SubAssign<&'a mut Self> for Block128FlatArk {
    #[inline] fn sub_assign(&mut self, rhs: &'a mut Self) { *self -= &*rhs; }
}

impl MulAssign<&Self> for Block128FlatArk {
    #[inline] fn mul_assign(&mut self, rhs: &Self) { self.0 = self.0 * rhs.0; }
}
impl DivAssign<&Self> for Block128FlatArk {
    #[inline]
    fn div_assign(&mut self, rhs: &Self) {
        let inv = rhs.inverse().expect("division by zero in Block128FlatArk");
        self.0 = self.0 * inv.0;
    }
}
impl MulAssign<Self> for Block128FlatArk {
    #[inline] fn mul_assign(&mut self, rhs: Self) { *self *= &rhs; }
}
impl DivAssign<Self> for Block128FlatArk {
    #[inline] fn div_assign(&mut self, rhs: Self) { *self /= &rhs; }
}
impl<'a> MulAssign<&'a mut Self> for Block128FlatArk {
    #[inline] fn mul_assign(&mut self, rhs: &'a mut Self) { *self *= &*rhs; }
}
impl<'a> DivAssign<&'a mut Self> for Block128FlatArk {
    #[inline] fn div_assign(&mut self, rhs: &'a mut Self) { *self /= &*rhs; }
}

impl Add<&Self> for Block128FlatArk { type Output = Self; #[inline] fn add(mut self, rhs: &Self) -> Self { self += rhs; self } }
impl Sub<&Self> for Block128FlatArk { type Output = Self; #[inline] fn sub(mut self, rhs: &Self) -> Self { self -= rhs; self } }
impl Mul<&Self> for Block128FlatArk { type Output = Self; #[inline] fn mul(mut self, rhs: &Self) -> Self { self *= rhs; self } }
impl Div<&Self> for Block128FlatArk { type Output = Self; #[inline] fn div(mut self, rhs: &Self) -> Self { self /= rhs; self } }

impl<'a> Add<&'a mut Self> for Block128FlatArk { type Output = Self; #[inline] fn add(self, rhs: &'a mut Self) -> Self { self + &*rhs } }
impl<'a> Sub<&'a mut Self> for Block128FlatArk { type Output = Self; #[inline] fn sub(self, rhs: &'a mut Self) -> Self { self - &*rhs } }
impl<'a> Mul<&'a mut Self> for Block128FlatArk { type Output = Self; #[inline] fn mul(self, rhs: &'a mut Self) -> Self { self * &*rhs } }
impl<'a> Div<&'a mut Self> for Block128FlatArk { type Output = Self; #[inline] fn div(self, rhs: &'a mut Self) -> Self { self / &*rhs } }

impl Add<Self> for Block128FlatArk { type Output = Self; #[inline] fn add(self, rhs: Self) -> Self { self + &rhs } }
impl Sub<Self> for Block128FlatArk { type Output = Self; #[inline] fn sub(self, rhs: Self) -> Self { self - &rhs } }
impl Mul<Self> for Block128FlatArk { type Output = Self; #[inline] fn mul(self, rhs: Self) -> Self { self * &rhs } }
impl Div<Self> for Block128FlatArk { type Output = Self; #[inline] fn div(self, rhs: Self) -> Self { self / &rhs } }

impl<'a, 'b> Add<&'b Block128FlatArk> for &'a Block128FlatArk { type Output = Block128FlatArk; fn add(self, rhs: &'b Block128FlatArk) -> Block128FlatArk { *self + rhs } }
impl<'a, 'b> Sub<&'b Block128FlatArk> for &'a Block128FlatArk { type Output = Block128FlatArk; fn sub(self, rhs: &'b Block128FlatArk) -> Block128FlatArk { *self - rhs } }
impl<'a, 'b> Mul<&'b Block128FlatArk> for &'a Block128FlatArk { type Output = Block128FlatArk; fn mul(self, rhs: &'b Block128FlatArk) -> Block128FlatArk { *self * rhs } }
impl<'a, 'b> Div<&'b Block128FlatArk> for &'a Block128FlatArk { type Output = Block128FlatArk; fn div(self, rhs: &'b Block128FlatArk) -> Block128FlatArk { *self / rhs } }

impl Sum<Self> for Block128FlatArk {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self { iter.fold(Self::zero(), |a, b| a + b) }
}
impl<'a> Sum<&'a Self> for Block128FlatArk {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self { iter.fold(Self::zero(), |a, b| a + b) }
}
impl Product<Self> for Block128FlatArk {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self { iter.fold(Self::one(), |a, b| a * b) }
}
impl<'a> Product<&'a Self> for Block128FlatArk {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self { iter.fold(Self::one(), |a, b| a * b) }
}

macro_rules! flat_from_uint {
    ($t:ty) => {
        impl From<$t> for Block128FlatArk {
            #[inline] fn from(n: $t) -> Self {
                Self::from_flat(Block128((n & 1) as u128).to_hardware())
            }
        }
    };
}
macro_rules! flat_from_int {
    ($t:ty) => {
        impl From<$t> for Block128FlatArk {
            #[inline] fn from(n: $t) -> Self {
                Self::from_flat(Block128((n & 1) as u128).to_hardware())
            }
        }
    };
}
flat_from_uint!(u8);
flat_from_uint!(u16);
flat_from_uint!(u32);
flat_from_uint!(u64);
flat_from_uint!(u128);
flat_from_int!(i8);
flat_from_int!(i16);
flat_from_int!(i32);
flat_from_int!(i64);
flat_from_int!(i128);
impl From<bool> for Block128FlatArk {
    #[inline] fn from(b: bool) -> Self {
        Self::from_flat(Block128(b as u128).to_hardware())
    }
}

impl ark_std::rand::distributions::Distribution<Block128FlatArk>
    for ark_std::rand::distributions::Standard
{
    fn sample<R: ark_std::rand::Rng + ?Sized>(&self, rng: &mut R) -> Block128FlatArk {
        let lo = rng.next_u64() as u128;
        let hi = rng.next_u64() as u128;
        Block128FlatArk::from_flat(Flat::from_raw(Block128(lo | (hi << 64))))
    }
}

impl Valid for Block128FlatArk {
    fn check(&self) -> Result<(), SerializationError> { Ok(()) }
}

impl CanonicalSerializeWithFlags for Block128FlatArk {
    fn serialize_with_flags<W: ark_std::io::Write, F: Flags>(
        &self, mut writer: W, flags: F,
    ) -> Result<(), SerializationError> {
        if F::BIT_SIZE > 8 { return Err(SerializationError::NotEnoughSpace); }
        let byte_size = buffer_byte_size(128 + F::BIT_SIZE);
        let mut bytes = ark_std::vec![0u8; byte_size];
        bytes[..16].copy_from_slice(&self.0.into_raw().0.to_le_bytes());
        bytes[byte_size - 1] |= flags.u8_bitmask();
        writer.write_all(&bytes).map_err(SerializationError::IoError)
    }
    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        buffer_byte_size(128 + F::BIT_SIZE)
    }
}

impl CanonicalSerialize for Block128FlatArk {
    fn serialize_with_mode<W: ark_std::io::Write>(
        &self, writer: W, _: Compress,
    ) -> Result<(), SerializationError> {
        self.serialize_with_flags(writer, EmptyFlags)
    }
    fn serialized_size(&self, _: Compress) -> usize {
        self.serialized_size_with_flags::<EmptyFlags>()
    }
}

impl CanonicalDeserializeWithFlags for Block128FlatArk {
    fn deserialize_with_flags<R: ark_std::io::Read, F: Flags>(
        mut reader: R,
    ) -> Result<(Self, F), SerializationError> {
        if F::BIT_SIZE > 8 { return Err(SerializationError::NotEnoughSpace); }
        let byte_size = Self::zero().serialized_size_with_flags::<F>();
        let mut bytes = ark_std::vec![0u8; byte_size];
        reader.read_exact(&mut bytes).map_err(SerializationError::IoError)?;
        let flag = F::from_u8_remove_flags(&mut bytes[byte_size - 1])
            .ok_or(SerializationError::UnexpectedFlags)?;
        let mut buf = [0u8; 16];
        buf.copy_from_slice(&bytes[..16]);
        Ok((Self::from_flat(Flat::from_raw(Block128(u128::from_le_bytes(buf)))), flag))
    }
}

impl CanonicalDeserialize for Block128FlatArk {
    fn deserialize_with_mode<R: ark_std::io::Read>(
        reader: R, _: Compress, _: Validate,
    ) -> Result<Self, SerializationError> {
        Self::deserialize_with_flags::<R, EmptyFlags>(reader).map(|(r, _)| r)
    }
}

impl AdditiveGroup for Block128FlatArk {
    type Scalar = Self;
    const ZERO: Self = Self(Flat::from_raw(Block128(0)));
    #[inline]
    fn double_in_place(&mut self) -> &mut Self {
        self.0 = Flat::from_raw(Block128(0));
        self
    }
    #[inline]
    fn neg_in_place(&mut self) -> &mut Self { self }
}

impl Field for Block128FlatArk {
    type BasePrimeField = Gf2;
    const SQRT_PRECOMP: Option<SqrtPrecomputation<Self>> = None;
    const ONE: Self = Self(Flat::from_raw(Block128(1)));
    const NEG_ONE: Self = Self(Flat::from_raw(Block128(1))); // -1 = 1 in characteristic 2

    fn extension_degree() -> u64 { 128 }

    fn characteristic() -> &'static [u64] { Gf2::characteristic() }

    fn from_base_prime_field(elem: Gf2) -> Self {
        Self::from_flat(Block128(elem.0.0[0] as u128).to_hardware())
    }

    fn to_base_prime_field_elements(&self) -> impl Iterator<Item = Gf2> {
        let tower_val = self.0.to_tower().0;
        (0..128usize).map(move |i| Gf2::from(((tower_val >> i) & 1) as u64))
    }

    fn from_base_prime_field_elems(elems: impl IntoIterator<Item = Gf2>) -> Option<Self> {
        let mut val = 0u128;
        let mut count = 0usize;
        for (i, e) in elems.into_iter().enumerate() {
            if i >= 128 { return None; }
            val |= (e.0.0[0] as u128) << i;
            count += 1;
        }
        if count != 128 { return None; }
        Some(Self::from_flat(Block128(val).to_hardware()))
    }

    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        if F::BIT_SIZE > 8 { return None; }
        let byte_size = buffer_byte_size(128 + F::BIT_SIZE);
        if bytes.len() < byte_size { return None; }
        let flag = F::from_u8(0)?;
        let mut buf = [0u8; 16];
        buf.copy_from_slice(&bytes[..16]);
        Some((Self::from_flat(Flat::from_raw(Block128(u128::from_le_bytes(buf)))), flag))
    }

    fn legendre(&self) -> LegendreSymbol {
        if self.is_zero() { LegendreSymbol::Zero } else { LegendreSymbol::QuadraticResidue }
    }

    fn sqrt(&self) -> Option<Self> {
        let mut r = *self;
        for _ in 0..127 { r.square_in_place(); }
        Some(r)
    }

    #[inline]
    fn square(&self) -> Self { Self::from_flat(self.0 * self.0) }

    #[inline]
    fn square_in_place(&mut self) -> &mut Self {
        self.0 = self.0 * self.0;
        self
    }

    #[inline]
    fn inverse(&self) -> Option<Self> {
        if self.is_zero() { return None; }
        Some(Self::from_flat(self.0.to_tower().invert().to_hardware()))
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        if self.is_zero() { return None; }
        self.0 = self.0.to_tower().invert().to_hardware();
        Some(self)
    }

    fn frobenius_map_in_place(&mut self, power: usize) {
        for _ in 0..power { self.square_in_place(); }
    }

    fn mul_by_base_prime_field(&self, elem: &Gf2) -> Self {
        if elem.is_zero() { Self::zero() } else { *self }
    }
}

#[cfg(test)]
mod flat_ark_tests {
    use super::*;
    use ark_ff::{Field, One, Zero};
    use ark_std::UniformRand;

    fn tower_to_flat(a: Block128Ark) -> Block128FlatArk {
        Block128FlatArk::from_flat(a.0.to_hardware())
    }

    #[test]
    fn zero_add_identity() {
        let mut rng = ark_std::test_rng();
        let a = Block128FlatArk::rand(&mut rng);
        assert_eq!(a + Block128FlatArk::zero(), a);
    }

    #[test]
    fn one_mul_identity() {
        let mut rng = ark_std::test_rng();
        let a = Block128FlatArk::rand(&mut rng);
        assert_eq!(a * Block128FlatArk::one(), a);
    }

    #[test]
    fn add_self_is_zero() {
        let mut rng = ark_std::test_rng();
        let a = Block128FlatArk::rand(&mut rng);
        assert!((a + a).is_zero());
    }

    #[test]
    fn mul_zero_is_zero() {
        let mut rng = ark_std::test_rng();
        let a = Block128FlatArk::rand(&mut rng);
        assert!((a * Block128FlatArk::zero()).is_zero());
    }

    #[test]
    fn inverse_roundtrip() {
        let mut rng = ark_std::test_rng();
        for _ in 0..100 {
            let a = Block128FlatArk::rand(&mut rng);
            if !a.is_zero() {
                assert_eq!(a * a.inverse().unwrap(), Block128FlatArk::one());
            }
        }
    }

    #[test]
    fn square_consistent_with_mul() {
        let mut rng = ark_std::test_rng();
        for _ in 0..100 {
            let a = Block128FlatArk::rand(&mut rng);
            assert_eq!(a.square(), a * a);
        }
    }

    #[test]
    fn consistent_with_tower() {
        let mut rng = ark_std::test_rng();
        for _ in 0..100 {
            let a_tower = Block128Ark::rand(&mut rng);
            let b_tower = Block128Ark::rand(&mut rng);

            let expected = tower_to_flat(a_tower * b_tower);
            let actual   = tower_to_flat(a_tower) * tower_to_flat(b_tower);

            assert_eq!(actual, expected, "flat mul inconsistent with tower mul");
        }
    }
}
