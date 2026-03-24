//! Config for GF(2^128) using generic algorithms — for benchmarking.

use crate::ark::{BinaryField, BinaryFieldConfig};
use crate::polynomials::{self, IrreduciblePoly};

/// Config for GF(2^128) using only generic algorithms.
///
/// Same polynomial as [`Gf128Config`] but does not override `mul_assign`,
/// `square_in_place`, or `inverse`.  Useful for benchmarking the generic
/// baseline against the optimized implementation.
///
/// [`Gf128Config`]: crate::ark::Gf128Config
pub struct Gf128GenericConfig;

impl Gf128GenericConfig {
    /// The irreducible polynomial that defines this field.
    pub const POLY: IrreduciblePoly = polynomials::GF128_GCM;
}

impl BinaryFieldConfig<2> for Gf128GenericConfig {
    // f(z) = z^128 + z^7 + z^2 + z + 1  (same polynomial as Gf128Config)
    const ALPHA_POW_M: [u64; 2] = [
        polynomials::GF128_GCM.alpha_pow_m[0],
        polynomials::GF128_GCM.alpha_pow_m[1],
    ];

    const ZERO: BinaryField<Self, 2> = BinaryField([0u64; 2], core::marker::PhantomData);
    const ONE: BinaryField<Self, 2> = BinaryField([1u64, 0u64], core::marker::PhantomData);

    // No overrides: all three operations use the generic defaults
    // (Algorithms 2.33 / 2.33 / 2.48).
}

/// GF(2^128) with generic algorithms, for benchmarking against [`Gf128`].
///
/// [`Gf128`]: crate::ark::Gf128
pub type Gf128Generic = BinaryField<Gf128GenericConfig, 2>;
