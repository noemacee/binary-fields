# binary-fields

Arithmetic in binary extension fields GF(2^m), implemented in Rust with an
[arkworks](https://arkworks.rs) `Field` trait integration.

Algorithms follow *Guide to Elliptic Curve Cryptography* (Hankerson, Menezes,
Vanstone), Chapter 2.3. All implementations are cross-validated against the
[`galois`](https://github.com/mhostetter/galois) Python library.

---

## Architecture

The crate is organized in three layers:

```
src/
├── polynomials.rs          — catalogue of irreducible polynomials
├── generic/                — polynomial-agnostic algorithms
│   ├── arithmetic.rs       — addition (2.32), multiplication (2.33)
│   ├── reduce.rs           — bit-at-a-time reduction (2.40)
│   └── invert.rs           — EEA (2.48) and binary GCD (2.49)
├── fields/                 — polynomial-specific optimized implementations
│   ├── z128_z7_z2_z1/      — GF(2^128) GCM polynomial, fully optimized
│   └── z233_z74_z1/        — GF(2^233) NIST polynomial (no overrides yet)
└── ark/                    — arkworks Field trait integration
    ├── gf2.rs              — GF(2) base field
    ├── binary_field.rs     — generic BinaryField<Config, N> and BinaryFieldConfig trait
    └── configs/            — one config per implemented field
        ├── gf128.rs        — Gf128Config / Gf128 (optimized)
        ├── gf128_generic.rs — Gf128GenericConfig / Gf128Generic (generic baseline)
        └── gf233.rs        — Gf233Config / Gf233
```

The **generic layer** works for any field degree and irreducible polynomial.
The **field-specific layer** provides optimized fast paths for concrete
polynomials. The **arkworks layer** exposes fields as `Field` trait
implementations via a config trait that selects which algorithms to use.

---

## Polynomial Catalogue

`src/polynomials.rs` is the single source of truth for all irreducible
polynomials used in the crate. Each `IrreduciblePoly` constant carries the
degree, non-leading exponents, the `ALPHA_POW_M` encoding (what α^m equals
after reduction), and a usage note.

| Constant       | Polynomial                          | Field size | Note                              |
|----------------|-------------------------------------|------------|-----------------------------------|
| `GF64`         | z^64 + z^4 + z^3 + z + 1           | GF(2^64)   | Lightweight crypto / LFSR         |
| `GF128_GCM`    | z^128 + z^7 + z^2 + z + 1          | GF(2^128)  | AES-GCM (NIST SP 800-38D)         |
| `GF163_NIST`   | z^163 + z^7 + z^6 + z^3 + 1        | GF(2^163)  | NIST B-163 / K-163                |
| `GF233_NIST`   | z^233 + z^74 + 1                    | GF(2^233)  | NIST B-233 / K-233                |
| `GF283_NIST`   | z^283 + z^12 + z^7 + z^5 + 1       | GF(2^283)  | NIST B-283 / K-283                |
| `GF409_NIST`   | z^409 + z^87 + 1                    | GF(2^409)  | NIST B-409 / K-409                |
| `GF571_NIST`   | z^571 + z^10 + z^5 + z^2 + 1       | GF(2^571)  | NIST B-571 / K-571                |

Field configs reference these constants directly — `ALPHA_POW_M` and `DEGREE`
are never hardcoded in the configs themselves.

---

## Usage

### Field arithmetic via arkworks

```rust
use binary_fields::ark::{Gf128, Gf233};
use ark_ff::{Field, One, Zero};

// Construct elements
let a = Gf128::new([0xdeadbeefcafe1234, 0xabad1dea12345678]);
let b = Gf128::new([0x1234567890abcdef, 0xfedcba9876543210]);

// Arithmetic
let sum     = a + b;          // XOR
let product = a * b;          // optimized comb multiplication
let square  = a.square();     // bit-expansion table squaring
let inv     = a.inverse();    // binary GCD inversion (None if a == 0)
let power   = a.pow([42u64]); // square-and-multiply
```

### Using the polynomial catalogue

```rust
use binary_fields::polynomials;

let poly = polynomials::GF128_GCM;
println!("{}", poly.name);   // "z^128 + z^7 + z^2 + z + 1"
println!("{}", poly.note);   // "GF(2^128) — GCM/AES-GCM ..."

// Iterate over all known polynomials
for p in polynomials::ALL {
    println!("degree {}: {}", p.degree, p.name);
}
```

### Using the generic algorithms directly

```rust
use binary_fields::generic::arithmetic::{add_2_32, mul_2_33};
use binary_fields::generic::invert::invert_2_49;

// Works for any field — supply the polynomial and degree
const POLY: [u64; 2] = [0x87, 0];  // GF(2^128) GCM
const DEG: usize = 128;

let a = [0xdeadbeefcafe1234u64, 0xabad1dea12345678];
let b = [0x1234567890abcdefu64, 0xfedcba9876543210];

let sum     = add_2_32(a, b);
let product = mul_2_33(&a, &b, &POLY, DEG);
let inv     = invert_2_49(&a, &POLY, DEG);
```

---

## Implemented Algorithms

### Generic (any field)

| Algorithm | Function           | Description                              |
|-----------|--------------------|------------------------------------------|
| 2.32      | `add_2_32`         | Addition — word-wise XOR                 |
| 2.33      | `mul_2_33`         | Right-to-left shift-and-add multiplication |
| 2.40      | `reduce_2_40`      | Bit-at-a-time modular reduction          |
| 2.48      | `invert_2_48`      | Inversion via extended Euclidean algorithm |
| 2.49      | `invert_2_49`      | Inversion via binary GCD                 |

### GF(2^128) optimized (`fields::z128_z7_z2_z1`)

| Algorithm | Function              | Description                                      |
|-----------|-----------------------|--------------------------------------------------|
| 2.34      | `mul_2_34`            | Right-to-left comb multiplication                |
| 2.39      | `square_2_39`         | Squaring via 512-byte bit-expansion table        |
| 2.41      | `reduce_2_41`         | Word-at-a-time reduction (14 XOR/shift ops)      |
| —         | `pow_rtl_2_34_2_39`   | Square-and-multiply exponentiation               |

The fast reduction (`reduce_2_41`) exploits two properties of the GCM
polynomial: m = 128 is word-aligned (W = 64), and deg(r) = 7 < W, so all
reduction terms stay within a single word boundary.

---

## Benchmarks

Run with `cargo bench`. Criterion generates an HTML report in
`target/criterion/`.

Latest results (Apple M-series, single core):

| Benchmark                  | Time     | Notes                              |
|----------------------------|----------|------------------------------------|
| `add_2_32`                 | ~0.9 ns  | Two XORs                           |
| `reduce_2_41`              | ~1.0 ns  | Word-at-a-time, 14 ops             |
| `reduce_2_40`              | ~175 ns  | Bit-at-a-time loop                 |
| `square_2_39`              | ~2.8 ns  | Lookup table + fast reduction      |
| `mul_2_34`                 | ~49 ns   | Comb + fast reduction              |
| `mul_2_33`                 | ~100 ns  | Baseline shift-and-add             |
| `invert_2_49`              | ~280 ns  | Binary GCD                         |
| `invert_2_48`              | ~553 ns  | Extended Euclidean                 |
| `pow_rtl_2_34_2_39`        | ~7 µs    | 128-bit exponent, square-and-mul   |
| `gf128/mul/optimized`      | ~49 ns   | Field-level, optimized config      |
| `gf128/mul/generic`        | ~100 ns  | Field-level, generic config        |
| `gf128/square/optimized`   | ~2.8 ns  |                                    |
| `gf128/square/generic`     | ~100 ns  |                                    |
| `gf128/inverse/optimized`  | ~280 ns  |                                    |
| `gf128/inverse/generic`    | ~553 ns  |                                    |

---

## Adding a New Field

1. Add an `IrreduciblePoly` constant to `src/polynomials.rs`.
2. Optionally add polynomial-specific optimizations under `src/fields/`.
3. Create a config in `src/ark/configs/`:

```rust
pub struct Gf64Config;

impl Gf64Config {
    pub const POLY: IrreduciblePoly = polynomials::GF64;
}

impl BinaryFieldConfig<1> for Gf64Config {
    const ALPHA_POW_M: [u64; 1] = [polynomials::GF64.alpha_pow_m[0]];
    const ZERO: BinaryField<Self, 1> = BinaryField([0u64; 1], core::marker::PhantomData);
    const ONE:  BinaryField<Self, 1> = BinaryField([1u64],    core::marker::PhantomData);
    // Override mul_assign / square_in_place / inverse for optimized variants.
}

pub type Gf64 = BinaryField<Gf64Config, 1>;
```

4. Re-export from `src/ark/configs/mod.rs` and `src/ark/mod.rs`.

---

## Reference

Hankerson, D., Menezes, A., Vanstone, S.
*Guide to Elliptic Curve Cryptography*, Chapter 2.3 — Binary Field Arithmetic.
Springer, 2004.
