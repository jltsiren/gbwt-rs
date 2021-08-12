//! # GBWT: Graph BWT
//!
//! This is a Rust reimplementation of parts of the [GBWT](https://github.com/jltsiren/gbwt) and the [GBWTGraph](https://github.com/jltsiren/gbwtgraph).
//! It is based on the [Simple-SDS](https://github.com/jltsiren/simple-sds) library.
//!
//! # References
//!
//! ### GBWT
//!
//! Jouni Sirén, Erik Garrison, Adam M. Novak, Benedict Paten, and Richard Durbin: **Haplotype-aware graph indexes**.\
//! Bioinformatics 36(2):400-407, 2020.
//! DOI: [10.1093/bioinformatics/btz575](https://doi.org/10.1093/bioinformatics/btz575)
//!
//! # Notes
//!
//! * See [Simple-SDS](https://github.com/jltsiren/simple-sds) for assumptions on the environment.
//! * This implementation supports the Simple-SDS file formats for [GBWT](https://github.com/jltsiren/gbwt/SERIALIZATION.md) and [GBZ](https://github.com/jltsiren/gbwtgraph/SERIALIZATION.md).

pub mod bwt;
pub mod gbwt;
pub mod headers;
pub mod support;

pub use crate::gbwt::{GBWT, SearchState, BidirectionalState};

// mod gbz: GBZ

//-----------------------------------------------------------------------------

/// Node identifier `0` is used for technical purposes and does not exist in the graph.
pub const ENDMARKER: usize = 0;

/// Key of the source tag.
pub const SOURCE_KEY: &str = "source";

/// Value of the source tag.
pub const SOURCE_VALUE: &str = "jltsiren/gbwt-rs";

//-----------------------------------------------------------------------------
