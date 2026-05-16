# GF(2) as SmallFp — Implementation Record

This document covers every change made to integrate GF(2) as a first-class
`SmallFp` field. All paths are relative to `arkworks/algebra/` for fork changes
and `binary-fields/` for the new module.

---

## Fork changes

### 1. `MODULUS_MINUS_ONE_DIV_TWO` — p=2 special case in `PrimeField for SmallFp<P>`

**File:** `ff/src/fields/models/small_fp/small_fp_backend.rs`

**Background:** `PrimeField` requires `MODULUS_MINUS_ONE_DIV_TWO` = `(p-1)/2`.
Used in `legendre()` as the Euler criterion exponent.

**Problem:** The existing computation:

```rust
const MODULUS_MINUS_ONE_DIV_TWO = Self::MODULUS.divide_by_2_round_down();
```

`divide_by_2_round_down`: if input is odd subtract 1 then shift, else just shift.

- p=251 (odd): `(251-1) >> 1 = 125` ✓
- p=2 (even): skip subtract, just shift → `2 >> 1 = 1` ✗ (should be 0)

**Fix:** Keep `divide_by_2_round_down` for all odd primes, add a special case
for the one even prime p=2:

```rust
const MODULUS_MINUS_ONE_DIV_TWO: Self::BigInt = if P::MODULUS_U128 == 2 {
    BigInt([0u64])  // (2-1)/2 = 0
} else {
    Self::MODULUS.divide_by_2_round_down()  // unchanged for all odd primes
};
```

No changes to `divide_by_2_round_down` itself. No new trait constant. No config
needs to do anything.

-> Keep in mind it will also never be used since the Legendre symbol is undefined for p=2, but it's cleaner to have it be correct anyway.

---

### 2. `TRACE` / `TRACE_MINUS_ONE_DIV_TWO` — p=2 special case in `PrimeField for SmallFp<P>`

**File:** `ff/src/fields/models/small_fp/small_fp_backend.rs`

**Problem:** `PrimeField for SmallFp<P>` computed `TRACE` by calling
`two_adic_coefficient` on `BigInt([MODULUS])`. This function asserts
`self.const_is_odd()` before doing anything. For `p=2`, `MODULUS = 2` is even
— the assert panics. It's lazy so it only fires if `Gf2::TRACE` is accessed,
but it's a landmine.

**Fix:** Same pattern as `MODULUS_MINUS_ONE_DIV_TWO` — keep `two_adic_coefficient`
for all odd primes, special-case p=2:

```rust
const TRACE: Self::BigInt = if P::MODULUS_U128 == 2 {
    BigInt([1u64])  // odd part of (p-1) = odd part of 1 = 1
} else {
    Self::MODULUS.two_adic_coefficient()  // unchanged for all odd primes
};

const TRACE_MINUS_ONE_DIV_TWO: Self::BigInt = if P::MODULUS_U128 == 2 {
    BigInt([0u64])  // (1-1)/2 = 0
} else {
    Self::TRACE.divide_by_2_round_down()  // unchanged for all odd primes
};
```

No new trait constant. No config needs to do anything.

Note: neither constant will ever be read for GF(2) in practice — `TRACE` feeds
Tonelli-Shanks sqrt (we use `CharTwo` precomp instead) and NTT roots of unity
(`TWO_ADICITY = 0` means no NTT). They are set correctly for trait compliance
only.

---

### 3. `fn legendre` — overridable method in `SmallFpConfig`

**File:** `ff/src/fields/models/small_fp/small_fp_backend.rs` (trait) and
`ff/src/fields/models/small_fp/field.rs` (delegation)

**Problem:** `Field::legendre` in `field.rs` inlined the Euler criterion
directly. There was no way to override it per config. For GF(2), the Legendre
symbol is mathematically undefined (the Euler criterion only applies to odd
primes), so GF(2) needs to `unimplemented!()` it.

**Fix:** Move the implementation into `SmallFpConfig` as a default method,
delegate from `Field::legendre`:

```rust
// Added to SmallFpConfig trait — default Euler criterion:
fn legendre(a: &SmallFp<Self>) -> LegendreSymbol {
    let s = a.pow([((Self::MODULUS_U128 - 1) / 2) as u64]);
    if s.is_zero() { LegendreSymbol::Zero }
    else if s.is_one() { LegendreSymbol::QuadraticResidue }
    else { LegendreSymbol::QuadraticNonResidue }
}

// field.rs — delegate instead of inline:
fn legendre(&self) -> LegendreSymbol {
    P::legendre(self)
}
```

**Result for each field:**

| Config                    | Explicit override? | Gets                              |
| ------------------------- | ------------------ | --------------------------------- |
| `SmallFp8Config` (p=251)  | No                 | trait default → Euler criterion ✓ |
| `SmallFp32BabybearConfig` | No                 | trait default → Euler criterion ✓ |
| `Gf2SmallFpConfig`        | **Yes**            | `unimplemented!()`                |

The macro-generated configs don't emit a `legendre` method — they inherit the
default. Only GF(2) overrides it. Zero changes to existing generated code.

---

### 4. `UniformRand` — removed `val > 0` guard

**File:** `ff/src/fields/models/small_fp/small_fp_backend.rs`

**Problem:** The rejection sampler had an extra guard:

```rust
if val > 0 && u128::from(val) < P::MODULUS_U128 { ... }
```

The `val > 0` condition means `0` is never returned. For GF(2) this is fatal —
half the field (element `0`) is unreachable. For other SmallFp fields it's a
latent bug (1/p probability of never sampling 0).

`Fp<P, N>` has no such guard — its sampler uses `!is_geq_modulus()` which
correctly allows 0. The guard is a SmallFp-only oversight.

**Fix:**

```rust
// Before:
if val > 0 && u128::from(val) < P::MODULUS_U128 {
// After:
if u128::from(val) < P::MODULUS_U128 {
```

**Result:** 0 is now reachable for all SmallFp fields. Consistent with `Fp`.

**Proof:** `test_fp_samples_zero` in `test-curves/src/smallfp.rs` verifies that
`SmallFp8` (p=251) samples 0 in 10k draws (~1/251 probability per draw).
For large primes like Goldilocks (p≈2^64) the probability per draw is ~1/2^64
— not practical to assert in a test, but the code path is identical and there
is no guard blocking it.

**Why not use `bool` as `T` for GF(2)?**

`bool` (`false = 0`, `true = 1`) seems like a natural fit, but it does not
satisfy the `SmallFpConfig::type T` bounds:

```rust
type T: Unsigned              // ✗ bool is not a numeric type
      + Add<Output = Self::T> // ✗ bool has no Add
      + Sub<Output = Self::T> // ✗ bool has no Sub
      + Mul<Output = Self::T> // ✗ bool has no Mul
      + Div<Output = Self::T> // ✗ bool has no Div
      + Rem<Output = Self::T> // ✗ bool has no Rem
      + Into<u128>            // ✗ bool does not impl Into<u128>
      + TryFrom<u128>         // ✗ bool does not impl TryFrom<u128>
```

The bounds are designed for integer primitives only. `u8` is the minimum.

**Can we encode in 1 bit instead of 8?**

No. Rust has no 1-bit integer primitive — `u8` is the smallest addressable
type. Every GF(2) element stored as `u8` wastes 7 bits. Packing 8 elements
into one `u8` (bitpacking) is possible but requires a completely different
architecture: `SmallFp` stores exactly one element per `T`, so packing would
break the model entirely. If dense bitpacking matters, it needs its own type
outside of `SmallFp`.

---

### 5. Test template fixes

**File:** `test-templates/src/fields.rs`

The `val > 0` fix exposed two latent bugs in the test templates — they assumed
`UniformRand` never returns 0, which was only true due to the now-fixed bug.

**`test_mul_properties` (line 248):**

```rust
// Before — panics when a == 0:
assert_eq!(a * a.inverse().unwrap(), one);

// After — inverse only defined for nonzero:
if !a.is_zero() { assert_eq!(a * a.inverse().unwrap(), one); }
```

**`test_sqrt` (line 355):**

```rust
// Before — legendre(0) = Zero, not QuadraticResidue:
let b = a.square();
assert_eq!(b.legendre(), LegendreSymbol::QuadraticResidue);

// After — 0² = 0, skip:
let b = a.square();
if !b.is_zero() { assert_eq!(b.legendre(), LegendreSymbol::QuadraticResidue); }
```

Both fixes are mathematically correct: inverse and Legendre symbol are only
defined for nonzero elements.

---

## New module: `src/ark/smallfp/`

### `src/ark/smallfp/gf2.rs`

Manual `SmallFpConfig` implementation for GF(2). No macro involved — all
arithmetic is direct XOR/AND with no Montgomery layer.

| Constant / method          | Value                          | Why                             |
| -------------------------- | ------------------------------ | ------------------------------- |
| `type T`                   | `u8`                           | smallest type fitting modulus 2 |
| `MODULUS`                  | `2u8`                          |                                 |
| `GENERATOR`                | `from_raw(1)`                  | GF(2)\* = {1}                   |
| `ZERO` / `ONE` / `NEG_ONE` | `0`, `1`, `1`                  | -1 = 1 in char 2                |
| `TWO_ADICITY`              | `0`                            | \|GF(2)\*\| = 1 = 2⁰ × 1        |
| `TWO_ADIC_ROOT_OF_UNITY`   | `from_raw(1)`                  | unique 1st root of unity        |
| `SQRT_PRECOMP`             | `Some(CharTwo)`                | sqrt(x) = x in GF(2)            |
| `add_assign`               | `a.value ^= b.value`           | XOR                             |
| `sub_assign`               | `a.value ^= b.value`           | same as add in char 2           |
| `double_in_place`          | `a.value = 0`                  | x+x=0 in char 2                 |
| `neg_in_place`             | no-op                          | -x=x in char 2                  |
| `mul_assign`               | `a.value &= b.value`           | AND                             |
| `square_in_place`          | no-op                          | x²=x in GF(2)                   |
| `inverse`                  | `None` if 0, `Some(*a)` if 1   |                                 |
| `new(val)`                 | `from_raw(val % 2)`            | no Montgomery encoding          |
| `from_bigint`              | check `< 2`, return `from_raw` |                                 |
| `into_bigint`              | `BigInt([a.value as u64])`     |                                 |
| `legendre`                 | `unimplemented!()`             | undefined in char 2             |

**`pub type Gf2SmallFp = SmallFp<Gf2SmallFpConfig>`**

---

## Test coverage

11 tests in `src/ark/smallfp/gf2.rs`:

| Test                        | Verifies                                                                   |
| --------------------------- | -------------------------------------------------------------------------- |
| `add_is_xor`                | full XOR truth table                                                       |
| `mul_is_and`                | full AND truth table                                                       |
| `neg_is_identity`           | -x = x in char 2                                                           |
| `square_is_identity`        | x² = x in GF(2)                                                            |
| `inverse`                   | 0→None, 1→Some(1)                                                          |
| `double_is_zero`            | x+x=0                                                                      |
| `from_int_mod2`             | 0,1,2,3 → 0,1,0,1                                                          |
| `characteristic_is_2`       | `[2u64]`                                                                   |
| `prime_field_constants`     | MODULUS=2, MODULUS_MINUS_ONE_DIV_TWO=0, TRACE=1, TRACE_MINUS_ONE_DIV_TWO=0 |
| `uniform_rand_samples_zero` | 0 is reachable in 1000 samples                                             |
| `serialization_roundtrip`   | 1 byte, 0↔0 and 1↔1                                                        |

---

## Summary of files changed

| File                                                | Type     | What changed                                                                                                                              |
| --------------------------------------------------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `ff/src/fields/models/small_fp/small_fp_backend.rs` | Fork     | Added `MODULUS_MINUS_ONE_DIV_TWO_CONST`, `TRACE_CONST`, `fn legendre` to trait; routed `PrimeField` consts through `P::`; fixed `val > 0` |
| `ff/src/fields/models/small_fp/field.rs`            | Fork     | Delegated `legendre()` to `P::legendre(self)`                                                                                             |
| `test-templates/src/fields.rs`                      | Fork     | Fixed two latent zero-element bugs in test templates                                                                                      |
| `src/ark/smallfp/gf2.rs`                            | New      | `Gf2SmallFpConfig` + `Gf2SmallFp` type + 11 tests                                                                                         |
| `src/ark/smallfp/mod.rs`                            | New      | Module exports                                                                                                                            |
| `src/ark/mod.rs`                                    | Modified | Added `pub mod smallfp` and re-exports                                                                                                    |
