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

## Implementation Roadmap

### ✅ Step 1 — Representation (Figure 2.4)
Define the `GF2_128` struct as `[u64; 2]`. Establish the bit layout convention and implement `zero()`, `one()`, `is_zero()`, and `new(lo, hi)`.

### ✅ Step 2 — Addition (Algorithm 2.32)
XOR each word independently. Also implement `sub()` (identical to addition in characteristic 2) and `neg()` (identity in characteristic 2).

### ✅ Step 3 — Multiplication (Algorithm 2.33)
Right-to-left shift-and-add field multiplication. Reduction mod f(z) is handled inline at each step via `mul_by_z()`: left-shift by 1, then XOR with `r(z) = 0x87` if the high bit spills over.

### ✅ Step 4 — Inversion (Algorithm 2.48)
Extended Euclidean algorithm over GF(2)[x]. Polynomials are represented as `[u64; 3]` internally to accommodate `f` (degree 128) and intermediate shifts. Returns `None` for zero. The binary algorithm (Algorithm 2.49) and almost inverse algorithm (Algorithm 2.50) are left as later optimizations.

### ⬜ Step 5 — Squaring (Algorithm 2.39)
Exploit the fact that squaring is a linear map in characteristic 2: insert a zero bit between consecutive bits of `a`. Implemented via a precomputed 512-byte lookup table on bytes.

### ⬜ Step 6 — Faster Multiplication (Algorithm 2.34 + Algorithm 2.40)
Replace Algorithm 2.33 with the right-to-left comb method (Algorithm 2.34) for polynomial multiplication, followed by separate reduction mod f(z) (Algorithm 2.40). Same result, significantly faster.

---

## Functions to Implement

### Arithmetic Core
- ✅ `add(a, b)` — XOR word by word (Algorithm 2.32)
- ✅ `sub(a, b)` — identical to addition in characteristic 2
- ✅ `neg(a)` — identity in characteristic 2
- ✅ `mul(a, b)` — shift-and-add with inline reduction (Algorithm 2.33)
- ⬜ `square(a)` — bit-expansion via lookup table (Algorithm 2.39)

### Field Structure
- ✅ `zero()` — additive identity
- ✅ `one()` — multiplicative identity
- ✅ `is_zero(a)` — true if `a` is the zero element

### Exponentiation
- ⬜ `pow(a, n)` — square-and-multiply
- ✅ `invert(a)` — extended Euclidean algorithm over GF(2)[x] (Algorithm 2.48)

### Utilities
- ⬜ `from_u64(val)` — construct a field element from a raw integer
- ⬜ `to_words(a)` — extract the raw `[u64; 2]` representation