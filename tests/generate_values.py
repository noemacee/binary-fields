import galois
import random

GF = galois.GF(2**128, irreducible_poly="x^128 + x^7 + x^2 + x + 1")

# Generate some test vectors
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
