# Binary Fields

Implementation of arithmetic in GF(2^128), following *Guide to Elliptic Curve Cryptography* (Hankerson, Menezes, Vanstone), Chapter 2.3.

Cross-validated against the [`galois`](https://github.com/mhostetter/galois) Python library.

---

## Target Parameters

- **Field:** GF(2^128)
- **Irreducible polynomial:** f(z) = z^128 + z^7 + z^2 + z + 1 (used in AES-GCM, NIST SP 800-38D)
- **Word size:** W = 64
- **Words per element:** t = 2, s = 0 (no unused bits)
- **Storage:** `[u64; 2]` — `A[0]` holds `a63..a0`, `A[1]` holds `a127..a64`

---

## Implementation

### ✅ Step 1 — Representation (Figure 2.4)
Define the `GF2_128` struct as `[u64; 2]`. Establish the bit layout convention and implement `zero()`, `one()`, `is_zero()`, and `new(lo, hi)`.

### ✅ Step 2 — Addition (Algorithm 2.32)
XOR each word independently. Also implement `sub()` (identical to addition in characteristic 2) and `neg()` (identity in characteristic 2).

### ✅ Step 3 — Multiplication, baseline (Algorithm 2.33)
Right-to-left shift-and-add field multiplication. Reduction mod f(z) is handled inline at each step via `mul_by_z()`: left-shift by 1, then XOR with `r(z) = 0x87` if the high bit spills over.

### ✅ Step 4 — Reduction, general (Algorithm 2.40)
General modular reduction one bit at a time. Takes a `[u64; 4]` intermediate polynomial of degree at most 254 and reduces it mod f(z) to `[u64; 2]`.

### ✅ Step 5 — Reduction, fast (`reduce_2_41`)
Word-at-a-time reduction specialized for f(z) = z^128 + z^7 + z^2 + z + 1 with W=64. Derived following the method of HMV §2.3.5 (analogous to Algorithms 2.41–2.45 but for W=64 instead of W=32).

The derivation exploits two structural properties of this polynomial:
- m = 128 is a multiple of W = 64: perfect word alignment, no masking step needed
- deg(r) = 7 < W = 64: all reduction terms land within one word of their source, keeping overflow shifts small (>>57, >>62, >>63)

Total cost: 14 XOR/shift operations versus 127 loop iterations in Algorithm 2.40.

### ✅ Step 6 — Squaring (Algorithm 2.39)
Exploit the fact that squaring is a linear map in characteristic 2: insert a zero bit between consecutive bits of `a`. Implemented via a precomputed 512-byte lookup table on bytes, followed by fast reduction.

### ✅ Step 7 — Multiplication, fast (Algorithm 2.34)
Right-to-left comb method. Processes all words of A column by column (one bit position k at a time across all words), shifting B by one bit per outer iteration instead of one bit per element bit. Followed by `reduce_128_gf2p128`.

Cost reduction vs Algorithm 2.33: outer loop shrinks from m=128 to W=64 iterations; B is shifted 63 times instead of 127 times. Followed by `reduce_2_41`.

### ✅ Step 8 — Inversion, EEA (Algorithm 2.48)
Extended Euclidean algorithm over GF(2)[x]. Polynomials are represented as `[u64; 3]` internally to accommodate `f` (degree 128) and intermediate shifts. Clears high-degree terms first (left-to-right). Returns `None` for zero.

### ✅ Step 9 — Inversion, binary (Algorithm 2.49)
Binary GCD algorithm, the polynomial analogue of Stein's algorithm. Clears low-degree terms first (right-to-left): checks divisibility by z via bit 0, divides by z via right-shift, and compares u vs v as integers rather than computing explicit degrees.

Division of a Bezout coefficient g by z mod f: if bit 0 of g is 0, right-shift directly; if bit 0 is 1, add f first (making it divisible by z since f has constant term 1), then right-shift. The z^128 term of f after shifting lands exactly at bit 127.

### ✅ Step 10 — Exponentiation
Square-and-multiply. Uses `square_2_39()` and `mul_2_34()`.

### ✅ Step 11 — Utilities
`from_u64()`, `to_words()` for constructing and inspecting field elements.

---

## Benchmarks

| Operation                    | Algorithm           | Time    | Notes                             |
|------------------------------|---------------------|---------|-----------------------------------|
| `add_2_32`                   | 2.32                | ~0.9 ns | Two XORs                          |
| `reduce_2_40`                | 2.40                | ~175 ns | Bit-at-a-time loop                |
| `reduce_2_41`                | 2.41                | ~1.0 ns | Word-at-a-time, 14 XOR/shifts     |
| `mul_2_33`                   | 2.33                | ~100 ns | Baseline, inline reduction        |
| `mul_2_34`                   | 2.34                | ~49 ns  | Comb method + fast reduction      |
| `square_2_39`                | 2.39                | ~2.8 ns | Lookup table + fast reduction     |
| `invert_2_48`                | 2.48 (EEA)          | ~553 ns | Left-to-right, variable shifts    |
| `invert_2_49`                | 2.49 (binary GCD)   | ~280 ns | Right-to-left, fixed 1-bit shifts |
| `pow_rtl_2_34_2_39`          | square-and-multiply | ~7 µs   | RTL binary, uses 2.34 + 2.39      |

Run benchmarks with `cargo bench`. Criterion generates an HTML report in `target/criterion/` including per-algorithm comparison plots.

---

## Functions

### Arithmetic Core
- ✅ `add_2_32(a, b)` — XOR word by word (Algorithm 2.32)
- ✅ `sub(a, b)` — identical to addition in characteristic 2
- ✅ `neg(a)` — identity in characteristic 2
- ✅ `mul_2_33(a, b)` — shift-and-add with inline reduction (Algorithm 2.33)
- ✅ `mul_2_34(a, b)` — right-to-left comb + fast reduction (Algorithm 2.34)
- ✅ `square_2_39(a)` — bit-expansion via lookup table + fast reduction (Algorithm 2.39)

### Reduction
- ✅ `reduce_2_40(c)` — general bit-at-a-time reduction (Algorithm 2.40)
- ✅ `reduce_2_41(c)` — word-at-a-time reduction specialized for W=64 (Algorithm 2.41)

### Inversion
- ✅ `invert_2_48(a)` — EEA over GF(2)[x], left-to-right (Algorithm 2.48)
- ✅ `invert_2_49(a)` — binary GCD, right-to-left (Algorithm 2.49)

### Field Structure
- ✅ `zero()` — additive identity
- ✅ `one()` — multiplicative identity
- ✅ `is_zero(a)` — true if `a` is the zero element

### Exponentiation
- ✅ `pow_rtl_2_34_2_39(a, n)` — right-to-left binary, uses `mul_2_34` + `square_2_39`

### Utilities
- ✅ `from_u64(val)` — construct a field element from a raw integer
- ✅ `to_words(a)` — extract the raw `[u64; 2]` representation

