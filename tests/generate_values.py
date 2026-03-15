import galois
import random

GF = galois.GF(2**128, irreducible_poly="x^128 + x^7 + x^2 + x + 1")
GF2 = galois.GF(2)

# Irreducible polynomial f(z) = z^128 + z^7 + z^2 + z + 1
F = galois.Poly([1] + [0]*120 + [1,0,0,0,0,0,1,1], field=GF2)


def int_to_poly(n, field=GF2):
    """Convert an integer to a GF(2) polynomial (LSB = lowest degree)."""
    if n == 0:
        return galois.Poly([0], field=field)
    bits = []
    while n:
        bits.append(n & 1)
        n >>= 1
    # galois.Poly expects coefficients from highest degree to lowest
    return galois.Poly(bits[::-1], field=field)


def poly_to_int(p):
    """Convert a GF(2) polynomial back to an integer (LSB = lowest degree)."""
    result = 0
    coeffs = p.coeffs.tolist()
    # coeffs[0] is the highest degree
    for i, c in enumerate(reversed(coeffs)):
        if c:
            result |= (1 << i)
    return result


def int_to_words2(n):
    """Convert a 128-bit integer to [lo, hi] u64 pair."""
    lo = n & 0xFFFFFFFFFFFFFFFF
    hi = (n >> 64) & 0xFFFFFFFFFFFFFFFF
    return lo, hi


def int_to_words4(n):
    """Convert a 256-bit integer to [c0, c1, c2, c3] u64 quadruple."""
    mask = 0xFFFFFFFFFFFFFFFF
    return (n & mask, (n >> 64) & mask, (n >> 128) & mask, (n >> 192) & mask)


def unreduced_mul(a_int, b_int):
    """Compute the unreduced product of two field elements as a 256-bit integer."""
    pa = int_to_poly(a_int)
    pb = int_to_poly(b_int)
    product = pa * pb  # polynomial multiplication over GF(2), no reduction
    return poly_to_int(product)


print("=" * 60)
print("MUL / INV / POW test vectors (galois ground truth)")
print("=" * 60)

pairs = [
    (0x1234, 0x5678),
    (0xDEADBEEFCAFE1234, 0xABCD1234),
    (0x1, 0x2),
    (0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF, 0x1),
]

for a_val, b_val in pairs:
    a = GF(a_val)
    b = GF(b_val)
    print(f"a={hex(a_val)}, b={hex(b_val)}")
    print(f"  mul = {hex(int(a * b))}")
    print(f"  inv = {hex(int(a ** -1))}")
    print()

a = GF(0xDEADBEEFCAFE1234)
print(f"pow(a, 0)  = {hex(int(a ** 0))}")
print(f"pow(a, 1)  = {hex(int(a ** 1))}")
print(f"pow(a, 3)  = {hex(int(a ** 3))}")
print(f"pow(a, 17) = {hex(int(a ** 17))}")

b = GF(0x1234567890ABCDEF1234567890ABCDEF)
print(f"pow(b, 5)  = {hex(int(b ** 5))}")

a = GF(0xDEADBEEFCAFE1234)
b = GF(0x1234567890ABCDEF1234567890ABCDEF)
print(f"square(a) = {hex(int(a ** 2))}")
print(f"square(b) = {hex(int(b ** 2))}")

print()
print("=" * 60)
print("REDUCTION test vectors: [u64; 4] input -> [u64; 2] output")
print("These are unreduced products reduced mod f(z).")
print("=" * 60)

reduce_cases = [
    (0x1234, 0x5678),
    (0xDEADBEEFCAFE1234, 0xABCD1234),
    (0x1234567890ABCDEF, 0xFEDCBA9876543210),
    (0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF),
    (0x1, 0x80000000000000000000000000000000),  # 1 * z^127
    (0x80000000000000000000000000000000,
     0x80000000000000000000000000000000),       # z^127 * z^127 = z^254
]

for a_int, b_int in reduce_cases:
    unreduced = unreduced_mul(a_int, b_int)
    c0, c1, c2, c3 = int_to_words4(unreduced)
    expected = int(GF(a_int) * GF(b_int))
    lo, hi = int_to_words2(expected)
    print(f"// a={hex(a_int)}, b={hex(b_int)}")
    print(f"// unreduced product = {hex(unreduced)}")
    print(f"input  = [{hex(c0)}, {hex(c1)}, {hex(c2)}, {hex(c3)}]")
    print(f"output = [{hex(lo)}, {hex(hi)}]")
    print()
