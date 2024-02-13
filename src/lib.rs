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
//! ### GBWTGraph
//!
//! Jouni Sirén, Jean Monlong, Xian Chang, Adam M. Novak, Jordan M. Eizenga, Charles Markello, Jonas A. Sibbesen, Glenn Hickey, Pi-Chuan Chang, Andrew Carroll, Namrata Gupta, Stacey Gabriel, Thomas W. Blackwell, Aakrosh Ratan, Kent D. Taylor, Stephen S. Rich, Jerome I. Rotter, David Haussler, Erik Garrison, and Benedict Paten:\
//! **Pangenomics enables genotyping of known structural variants in 5202 diverse genomes**.\
//! Science 374(6574):abg8871, 2021.
//! DOI: [10.1126/science.abg8871](https://doi.org/10.1126/science.abg8871)
//!
//! ### GBZ
//!
//! Jouni Sirén and Benedict Paten: **GBZ file format for pangenome graphs**.\
//! Bioinformatics 38(22):5012-5018, 2022.
//! DOI: [10.1093/bioinformatics/btac656](https://doi.org/10.1093/bioinformatics/btac656)
//!
//! # Notes
//!
//! * See [Simple-SDS](https://github.com/jltsiren/simple-sds) for assumptions on the environment.
//! * This implementation supports the Simple-SDS file formats for [GBWT](https://github.com/jltsiren/gbwt/blob/master/SERIALIZATION.md) and [GBZ](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md).
//! * GBWT / GBZ files written by this library can be identified by `source` tag value `jltsiren/gbwt-rs`.

pub mod bwt;
pub mod gbwt;
pub mod gbz;
pub mod graph;
pub mod headers;
pub mod support;

// Shared internal code for the binaries.
#[cfg(feature = "binaries")]
#[doc(hidden)]
pub mod internal;

//-----------------------------------------------------------------------------

pub use crate::bwt::Pos;
pub use crate::gbwt::{GBWT, SearchState, BidirectionalState, Metadata, PathName, FullPathName};
pub use crate::gbz::GBZ;
pub use crate::graph::Segment;
pub use crate::support::Orientation;

//-----------------------------------------------------------------------------

/// Node identifier `0` is used for technical purposes and does not exist in the graph.
pub const ENDMARKER: usize = 0;

/// Key of the source tag.
pub const SOURCE_KEY: &str = "source";

/// Value of the source tag.
pub const SOURCE_VALUE: &str = "jltsiren/gbwt-rs";

/// Sample name for generic named paths.
///
/// This sample name is used for GFA P-lines with opaque string names. The actual name
/// is stored as a contig name.
pub const REF_SAMPLE: &str = "_gbwt_ref";

/// Key for the tag listing the names of reference samples.
pub const REFERENCE_SAMPLES_KEY: &str = "reference_samples";

//-----------------------------------------------------------------------------
