//! Catalogue of irreducible polynomials over GF(2) for binary extension fields.
//!
//! Each constant is an [`IrreduciblePoly`] carrying the polynomial degree, its
//! non-leading exponents, the `ALPHA_POW_M` encoding used by
//! [`BinaryFieldConfig`], a human-readable name, and a usage note.
//!
//! [`BinaryFieldConfig`]: crate::ark::BinaryFieldConfig

/// An irreducible polynomial over GF(2) that defines a binary extension field.
///
/// The polynomial is always monic (`z^degree + ...`). `terms` lists the
/// exponents of all remaining non-zero terms in descending order, including
/// the constant exponent 0 when the polynomial has a constant term.
///
/// `alpha_pow_m` encodes α^degree in the polynomial basis as little-endian
/// u64 limbs: `alpha_pow_m[i]` holds bits `[64i .. 64i+63]`.  This is exactly
/// the value required by [`BinaryFieldConfig::ALPHA_POW_M`].
///
/// [`BinaryFieldConfig::ALPHA_POW_M`]: crate::ark::BinaryFieldConfig::ALPHA_POW_M
#[derive(Debug, Clone, Copy)]
pub struct IrreduciblePoly {
    /// Extension degree: the polynomial is `z^degree + (lower terms)`.
    pub degree: usize,
    /// Exponents of all non-leading non-zero terms, descending.
    /// E.g. `z^128 + z^7 + z^2 + z + 1` → `&[7, 2, 1, 0]`.
    pub terms: &'static [usize],
    /// α^degree as little-endian u64 limbs (the `ALPHA_POW_M` encoding).
    pub alpha_pow_m: &'static [u64],
    /// Polynomial in mathematical notation.
    pub name: &'static str,
    /// Standard name and usage context.
    pub note: &'static str,
}

// ─── GF(2^64) ────────────────────────────────────────────────────────────────

/// z^64 + z^4 + z^3 + z + 1
///
/// Primitive pentanomial for GF(2^64), used in lightweight cryptographic
/// primitives and 64-bit LFSR constructions.
///
/// α^64 = α^4 + α^3 + α + 1 → `ALPHA_POW_M = [0x1B]`
pub const GF64: IrreduciblePoly = IrreduciblePoly {
    degree: 64,
    terms: &[4, 3, 1, 0],
    alpha_pow_m: &[0x1B],
    name: "z^64 + z^4 + z^3 + z + 1",
    note: "GF(2^64) — lightweight crypto, 64-bit LFSR / GHASH variants",
};

// ─── GF(2^128) ───────────────────────────────────────────────────────────────

/// z^128 + z^7 + z^2 + z + 1
///
/// The GCM polynomial (NIST SP 800-38D).  Used in AES-GCM, GHASH
/// authentication, and binary-field ZK/MPC protocols over 128-bit fields.
/// This is the polynomial implemented with optimized algorithms in
/// [`crate::fields::z128_z7_z2_z1`] and exposed as [`crate::ark::Gf128`].
///
/// α^128 = α^7 + α^2 + α + 1 → `ALPHA_POW_M = [0x87, 0]`
pub const GF128_GCM: IrreduciblePoly = IrreduciblePoly {
    degree: 128,
    terms: &[7, 2, 1, 0],
    alpha_pow_m: &[0x87, 0],
    name: "z^128 + z^7 + z^2 + z + 1",
    note: "GF(2^128) — GCM/AES-GCM (NIST SP 800-38D), binary-field MPC/ZK",
};

// ─── NIST binary curve fields (FIPS 186-3, Appendix D) ───────────────────────
// https://csrc.nist.gov/files/pubs/fips/186-3/final/docs/fips_186-3.pdf

// The use of binary fields in curves have been deprecated
// Proof of deprecation : https://csrc.nist.gov/news/2023/nist-releases-fips-186-5-and-sp-800-186
/// z^163 + z^7 + z^6 + z^3 + 1
///
/// Field polynomial for NIST curves B-163 and K-163.
///
/// α^163 = α^7 + α^6 + α^3 + 1 → `ALPHA_POW_M = [0xC9, 0, 0]`
pub const GF163_NIST: IrreduciblePoly = IrreduciblePoly {
    degree: 163,
    terms: &[7, 6, 3, 0],
    alpha_pow_m: &[0xC9, 0, 0],
    name: "z^163 + z^7 + z^6 + z^3 + 1",
    note: "GF(2^163) — NIST B-163 / K-163 (FIPS 186-4 §D.1.2)",
};

/// z^233 + z^74 + 1
///
/// Field polynomial for NIST curves B-233 and K-233.
/// Implemented as [`crate::ark::Gf233`].
///
/// α^233 = α^74 + 1 → `ALPHA_POW_M = [1, 1<<10, 0, 0]`
pub const GF233_NIST: IrreduciblePoly = IrreduciblePoly {
    degree: 233,
    terms: &[74, 0],
    alpha_pow_m: &[1, 1 << 10, 0, 0],
    name: "z^233 + z^74 + 1",
    note: "GF(2^233) — NIST B-233 / K-233 (FIPS 186-4 §D.1.3)",
};

/// z^283 + z^12 + z^7 + z^5 + 1
///
/// Field polynomial for NIST curves B-283 and K-283.
///
/// α^283 = α^12 + α^7 + α^5 + 1 → `ALPHA_POW_M = [0x10A1, 0, 0, 0, 0]`
pub const GF283_NIST: IrreduciblePoly = IrreduciblePoly {
    degree: 283,
    terms: &[12, 7, 5, 0],
    alpha_pow_m: &[0x10A1, 0, 0, 0, 0],
    name: "z^283 + z^12 + z^7 + z^5 + 1",
    note: "GF(2^283) — NIST B-283 / K-283 (FIPS 186-4 §D.1.4)",
};

/// z^409 + z^87 + 1
///
/// Field polynomial for NIST curves B-409 and K-409.
///
/// α^409 = α^87 + 1 → `ALPHA_POW_M = [1, 0x80_0000, 0, 0, 0, 0, 0]`
pub const GF409_NIST: IrreduciblePoly = IrreduciblePoly {
    degree: 409,
    terms: &[87, 0],
    alpha_pow_m: &[1, 0x80_0000, 0, 0, 0, 0, 0],
    name: "z^409 + z^87 + 1",
    note: "GF(2^409) — NIST B-409 / K-409 (FIPS 186-4 §D.1.5)",
};

/// z^571 + z^10 + z^5 + z^2 + 1
///
/// Field polynomial for NIST curves B-571 and K-571.
///
/// α^571 = α^10 + α^5 + α^2 + 1 → `ALPHA_POW_M = [0x425, 0, 0, 0, 0, 0, 0, 0, 0]`
pub const GF571_NIST: IrreduciblePoly = IrreduciblePoly {
    degree: 571,
    terms: &[10, 5, 2, 0],
    alpha_pow_m: &[0x425, 0, 0, 0, 0, 0, 0, 0, 0],
    name: "z^571 + z^10 + z^5 + z^2 + 1",
    note: "GF(2^571) — NIST B-571 / K-571 (FIPS 186-4 §D.1.6)",
};

// ─── Full catalogue ───────────────────────────────────────────────────────────

/// Every polynomial in this catalogue, ordered by degree.
pub const ALL: &[IrreduciblePoly] = &[
    GF64, GF128_GCM, GF163_NIST, GF233_NIST, GF283_NIST, GF409_NIST, GF571_NIST,
];
