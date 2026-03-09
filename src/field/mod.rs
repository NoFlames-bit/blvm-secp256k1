//! Field element arithmetic for secp256k1.
//!
//! Field modulus: p = 2^256 - 2^32 - 977 (SEC2 secp256k1)
//!
//! Platform layout:
//! - x86_64, aarch64: 5x52 (u64 limbs), pure Rust with u128 wide mul
//! - arm (32-bit): 10x26 (u32 limbs), ASM for mul/sqr

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
mod layout_5x52;

#[cfg(target_arch = "arm")]
mod layout_10x26;

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
pub use layout_5x52::{FeStorage, FieldElement};

#[cfg(target_arch = "arm")]
pub use layout_10x26::FieldElement;
