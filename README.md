# Binary Fields

Implementation of arithmetic in GF(2^128), following *Guide to Elliptic Curve Cryptography* (Hankerson, Menezes, Vanstone), Chapter 2.3.

---

## Target Parameters

- **Field:** GF(2^128)
- **Word size:** W = 64
- **Words per element:** t = 2, s = 0 (no unused bits)
- **Storage:** `[u64; 2]` — `A[0]` holds `a63..a0`, `A[1]` holds `a127..a64`

---

## Implementation Roadmap

### Step 1 — Representation (Figure 2.4)
Define the `GF2_128` struct as `[u64; 2]`. Establish the bit layout convention and add `zero()`, `one()`, and `is_zero()`.

### Step 2 — Addition (Algorithm 2.32)
XOR each word independently. Also implement `sub()` (identical to addition in characteristic 2) and `neg()` (identity).

### Step 3 — Polynomial Multiplication (Algorithm 2.34)
Right-to-left comb method. Produces the unreduced product `a(z) * b(z)` as a polynomial of degree up to `2m - 2`, stored in a `[u64; 4]` intermediate.

### Step 4 — Reduction (Algorithm 2.40, then fast variant)
Reduce the intermediate polynomial modulo `f(z)`. Start with Algorithm 2.40 (one bit at a time, works for any `f`). Once correct, replace with the fast word-at-a-time reduction for the chosen irreducible polynomial.

### Step 5 — Squaring (Algorithm 2.39)
Exploit the fact that squaring is a linear map in characteristic 2: insert a zero bit between each pair of consecutive bits of `a`. Implemented via a precomputed 512-byte lookup table on bytes.

### Step 6 — Inversion (Algorithm 2.48)
Extended Euclidean algorithm over GF(2)[x]. Finds `a^{-1} mod f` directly. The binary algorithm (Algorithm 2.49) and almost inverse algorithm (Algorithm 2.50) are left as later optimizations.

---

## Functions to Implement

### Arithmetic Core
- `add(a, b)` — XOR word by word (Algorithm 2.32)
- `sub(a, b)` — identical to addition in characteristic 2
- `neg(a)` — identity in characteristic 2
- `mul(a, b)` — polynomial multiplication (Algorithm 2.34) followed by reduction (Algorithm 2.40)
- `square(a)` — bit-expansion via lookup table (Algorithm 2.39)

### Field Structure
- `zero()` — additive identity
- `one()` — multiplicative identity
- `is_zero(a)` — true if `a` is the zero element

### Exponentiation
- `pow(a, n)` — square-and-multiply
- `invert(a)` — extended Euclidean algorithm over GF(2)[x] (Algorithm 2.48)

### Utilities
- `from_u64(val)` — construct a field element from a raw integer
- `to_words(a)` — extract the raw `[u64; 2]` representation