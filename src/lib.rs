//! # GBWT: Graph BWT
//!
//! This is a Rust reimplementation of parts of the [GBWT](https://github.com/jltsiren/gbwt) and the [GBWTGraph](https://github.com/jltsiren/gbwtgraph).
//! It is based on the [Simple-SDS](https://github.com/jltsiren/simple-sds) library.
//!
//! # Notes
//!
//! * See [Simple-SDS](https://github.com/jltsiren/simple-sds) for assumptions on the environment.
//! * This implementation supports the Simple-SDS file formats for [GBWT](https://github.com/jltsiren/gbwt/SERIALIZATION.md) and [GBZ](https://github.com/jltsiren/gbwtgraph/SERIALIZATION.md).

pub mod bwt;
pub mod support;

// mod gbwt: GBWT
// mod gbz: GBZ
