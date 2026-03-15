# `ark_compat` — arkworks Field Integration for Binary Fields

This module integrates the binary field arithmetic from the rest of this crate
into the [arkworks](https://github.com/arkworks-rs) `Field` trait ecosystem.
It exposes two concrete types (`Gf2` and `Gf128`) and two generic building
blocks (`BinaryFieldConfig` / `BinaryField`) that can instantiate any binary
extension field GF(2^m).

---

## Files

| File | Contents |
|------|----------|
| `gf2.rs` | GF(2) base field — `Gf2Config` + `type Gf2 = Fp<Gf2Config, 1>` |
| `binary_field.rs` | Generic GF(2^m) — `BinaryFieldConfig<N>` trait + `BinaryField<P, N>` struct |
| `gf128_config.rs` | GF(2^128) concrete config — `Gf128Config` + `type Gf128 = BinaryField<Gf128Config, 2>` |
| `mod.rs` | Re-exports: `Gf2`, `Gf128`, `BinaryField`, `BinaryFieldConfig` |

---

## Architecture

The design mirrors the pattern arkworks uses for prime fields:

```
Prime fields               Binary extension fields
───────────────────────    ──────────────────────────────────
MontConfig                 BinaryFieldConfig<N>   (your trait)
     │                           │
MontgomeryBackend          BinaryField<P, N>      (generic shell)
     │                           │
Fp<P, N>   ──── implements ──►  Field
     │
 PrimeField
```

For GF(2) specifically, since it is itself a prime field:

```
Gf2Config  ──── implements ──►  FpConfig<1>
                                      │
                           Fp<Gf2Config, 1>  ──── implements ──►  PrimeField
                                      │
                           type Gf2 = Fp<Gf2Config, 1>
```

---

## `Gf2` — the base field GF(2)

**File:** `gf2.rs`

```rust
pub type Gf2 = Fp<Gf2Config, 1>;
```

`Gf2Config` implements arkworks' `FpConfig<1>` with:

| Operation | Implementation |
|-----------|---------------|
| Addition  | XOR (`^`) |
| Subtraction | XOR (same as addition in char 2) |
| Multiplication | AND (`&`) |
| Negation | Identity (`-x = x` in char 2) |
| Doubling | Always 0 (`x + x = 0` in char 2) |
| Squaring | Identity (`x² = x` in char 2, Frobenius) |
| Inverse | `0 → None`, `1 → Some(1)` |

Because `Gf2` is a proper `Fp<P, 1>`, all of the following come for free
from the `Fp<P, N>` machinery without any additional code:

- All arithmetic operator overloads (`+`, `-`, `*`, `/`, `+=`, …)
- `From<u8/u16/u32/u64/u128/i8/…/bool>` — ring homomorphism, maps `n` to `n mod 2`
- `CanonicalSerialize` / `CanonicalDeserialize` — 1 byte on the wire
- `UniformRand`
- `PrimeField`, `FftField`, `Field`, `AdditiveGroup`, `Zero`, `One`
- `Ord`, `Hash`, `Display`, `Debug`, `Zeroize`, `FromStr`

**Key constants for `FpConfig<1>`:**

```
MODULUS                      = 2
ZERO                         = Fp([0])
ONE                          = Fp([1])
GENERATOR                    = 1   (trivial multiplicative group {1})
TWO_ADICITY                  = 0   (|GF(2)*| = 1 = 2⁰ × 1)
TWO_ADIC_ROOT_OF_UNITY       = 1
SQRT_PRECOMP                 = TonelliShanks { two_adicity: 0, … }
```

**Why `FpConfig<1>` instead of a standalone struct?**

A standalone `struct Gf2(u8)` would require manually implementing every
trait in the hierarchy (`AdditiveGroup`, `Field`, `FftField`, `PrimeField`,
`CanonicalSerialize`, all operator overloads, …). Using `FpConfig<1>` reduces
the implementation to just the 9 arithmetic hooks, with everything else
inherited from `Fp<P, N>`.

---

## `BinaryFieldConfig<N>` and `BinaryField<P, N>` — generic GF(2^m)

**File:** `binary_field.rs`

### `BinaryFieldConfig<N>` trait

The config trait a concrete binary extension field must implement:

```rust
pub trait BinaryFieldConfig<const N: usize>: Send + Sync + 'static + Sized {
    /// Extension degree over GF(2). Defaults to N * 64.
    /// Override only for non-round degrees (e.g. GF(2^127)).
    const DEGREE: usize = N * 64;

    const ZERO:    BinaryField<Self, N>;
    const ONE:     BinaryField<Self, N>;
    const NEG_ONE: BinaryField<Self, N>;   // = ONE in char 2
    const GENERATOR: BinaryField<Self, N>;

    fn mul_assign(a: &mut [u64; N], b: &[u64; N]);
    fn square_in_place(a: &mut [u64; N]);
    fn inverse(a: &[u64; N]) -> Option<[u64; N]>;
}
```

`N` is the number of 64-bit limbs: `N = 2` for GF(2^128), `N = 4` for
GF(2^256), etc. The config only needs to supply three arithmetic methods and
four constants; everything else is derived generically.

### `BinaryField<P, N>` struct

The generic shell that turns a `BinaryFieldConfig` into a full arkworks `Field`:

```rust
pub struct BinaryField<P: BinaryFieldConfig<N>, const N: usize>(
    pub [u64; N],
    pub PhantomData<P>,
);
```

Limb layout: `limbs[0]` holds coefficients for powers z^0…z^63 (LSB = z^0),
`limbs[N-1]` holds the highest-degree coefficients.

**Traits implemented generically for all `BinaryField<P, N>`:**

| Trait | Notes |
|-------|-------|
| `AdditiveGroup` | `double_in_place → 0`, `neg_in_place → identity` |
| `Field` | `BasePrimeField = Gf2`, extension degree = `P::DEGREE` |
| `CanonicalSerialize` / `Deserialize` | `ceil(DEGREE/8)` bytes, little-endian limbs |
| `UniformRand` | Uniform over all valid field elements |
| All arithmetic operators | `+`, `-`, `*`, `/`, `+=`, … (all reference variants) |
| `Sum`, `Product` | Iterator folding |
| `From<integer>` | Ring homomorphism `n → n mod 2` (0 or 1) |
| `Ord`, `Hash`, `Debug`, `Display`, `Zeroize` | — |

**Field methods and their meaning for binary fields:**

| Method | Behavior |
|--------|----------|
| `add` / `sub` | XOR of limbs |
| `neg` / `neg_in_place` | Identity (−x = x in char 2) |
| `double` / `double_in_place` | Always zero |
| `mul` | Delegates to `P::mul_assign` (carry-less multiplication mod irreducible poly) |
| `square` | Delegates to `P::square_in_place` |
| `inverse` | Delegates to `P::inverse` |
| `sqrt` | Computes a^(2^(m-1)) via m−1 squarings (always exists in GF(2^m)) |
| `frobenius_map_in_place(k)` | x → x^(2^k) via k squarings |
| `legendre` | Always `QuadraticResidue` for nonzero elements (Frobenius is bijective) |
| `characteristic` | Returns `&[2]` |
| `from_base_prime_field(g: Gf2)` | Embeds 0 or 1 as the constant polynomial term |
| `to_base_prime_field_elements` | Yields 128 `Gf2` bits (LSB first) |
| `from_base_prime_field_elems` | Reconstructs element from exactly `DEGREE` bits |

**`DEGREE` constant:**

`DEGREE` defaults to `N * 64`, which is correct for GF(2^128) with N=2,
GF(2^256) with N=4, etc. It can be overridden for non-round degrees:

```rust
// GF(2^127) would set:
const DEGREE: usize = 127;  // override; top bit of limbs[1] is always 0
```

---

## `Gf128Config` and `Gf128` — GF(2^128)

**File:** `gf128_config.rs`

```rust
pub type Gf128 = BinaryField<Gf128Config, 2>;
```

Implements `BinaryFieldConfig<2>` for the NIST GCM polynomial:

```
f(z) = z^128 + z^7 + z^2 + z + 1  (reduction constant 0x87)
```

| Method | Algorithm | Source |
|--------|-----------|--------|
| `mul_assign` | Comb method (Alg. 2.34) + fast reduction (Alg. 2.41) | `crate::arithmetic::mul_2_34`, `crate::reduce` |
| `square_in_place` | Bit-expansion table (Alg. 2.39) + fast reduction (Alg. 2.41) | `crate::arithmetic::square_2_39` |
| `inverse` | Binary GCD (Alg. 2.49) | `crate::invert::invert_2_49` |

All three delegate to the low-level algorithms in the parent crate, which are
independently tested against NIST and cross-validated reference vectors.

**Wire format:** 16 bytes, little-endian (limbs[0] first).

---

## Adding a new binary extension field

To instantiate GF(2^m) for a new m:

1. Choose N = ⌈m/64⌉ limbs.
2. Implement `BinaryFieldConfig<N>` for a new config struct, supplying:
   - `mul_assign` — carry-less multiplication mod your irreducible polynomial
   - `square_in_place` — squaring (can use generic squaring if no fast path)
   - `inverse` — extended Euclidean / binary GCD
   - `ZERO`, `ONE`, `NEG_ONE`, `GENERATOR`
   - Optionally override `DEGREE` if m is not a multiple of 64
3. Define the type alias:
   ```rust
   pub type MyField = BinaryField<MyFieldConfig, N>;
   ```

All trait impls (`Field`, serialization, operators, iterators, …) are inherited
automatically.

---

## Test coverage

Tests are located in each file's `#[cfg(test)]` module:

- **`gf2.rs`** — `add_is_xor`, `mul_is_and`, `neg_is_identity`,
  `square_is_identity`, `inverse`, `from_int_mod2`, `double_is_zero`,
  `serialization_roundtrip`, `prime_field_roundtrip`

- **`gf128_config.rs`** — field axioms (`zero`, `one`, `add_self_is_zero`,
  `mul_by_one`, `mul_by_zero`, `mul_associative`), `inverse_roundtrip`,
  `square_consistent_with_mul` (both single and two-limb variants),
  `frobenius_is_squaring`, `frobenius_power_3`, `sqrt_roundtrip`,
  `from_integer_is_mod2`, `characteristic_is_2`, `extension_degree_is_128`,
  `serialization_roundtrip` (16 bytes), `from_random_bytes_roundtrip`,
  `to_base_prime_field_roundtrip`, `from_base_prime_field_elems_wrong_count_is_none`,
  `mul_by_base_prime_field_{zero,one}`, `legendre_{zero,nonzero}`,
  `double_in_place_is_zero`, `neg_in_place_is_identity`,
  `inverse_in_place_{roundtrip,zero}`, `pow_is_consistent`,
  `sum_of_elements`, `product_of_elements`, `uniform_rand_not_always_zero`
