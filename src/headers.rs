//! File format headers.

use std::path::Path;

use simple_sds::serialize::{Serialize, Serializable};

//-----------------------------------------------------------------------------

/// Common functionality for file format headers.
///
/// This struct contains the following fields: `tag`, `version`, and `flags`.
/// The payload type contains the remaining fields.
///
/// # Examples
///
/// ```
/// use gbz::headers::{Header, Payload};
/// use simple_sds::serialize::Serialize;
///
/// #[derive(Copy, Clone, Default, PartialEq, Eq)]
/// struct Example {
///     data: u64,
/// }
///
/// impl Example {
///     const FLAG: u64 = 0x1;
/// }
/// 
/// impl Payload for Example {
///     const NAME: &'static str = "Example";
///     const TAG: u32 = 1234567890;
///     const VERSION: u32 = 1;
///     const MIN_VERSION: u32 = 1;
///     const DEFAULT_FLAGS: u64 = 0;
///
///     fn update(&mut self) {}
///
///     fn mask(_: u32) -> u64 {
///         0x1
///     }
///
///     fn validate(_: &Header<Self>) -> Result<(), String> {
///         Ok(())
///     }
/// }
///
/// let mut header = Header::<Example>::default();
/// assert_eq!(header.size_in_elements(), 3);
/// header.set(Example::FLAG);
/// assert!(header.validate().is_ok());
/// assert!(header.is_set(Example::FLAG));
/// header.unset(Example::FLAG);
/// assert!(!header.is_set(Example::FLAG));
/// ```
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Header<T: Payload> {
    tag: u32,
    version: u32,
    payload: T,
    flags: u64,
}

impl<T: Payload> Header<T> {
    /// Creates a default header.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the file format version in the header.
    #[inline]
    pub fn version(&self) -> u32 {
        self.version
    }

    /// Updates the header to the latest version.
    pub fn update(&mut self) {
        self.version = T::VERSION;
        self.payload.update()
    }

    /// Returns `true` if the specified binary flag is set.
    #[inline]
    pub fn is_set(&self, flag: u64) -> bool {
        (self.flags & flag) != 0
    }

    /// Sets the specified binary flag.
    #[inline]
    pub fn set(&mut self, flag: u64) {
        self.flags |= flag;
    }

    /// Unsets the specified binary flag.
    #[inline]
    pub fn unset(&mut self, flag: u64) {
        self.flags &= !flag;
    }

    /// Validates the header and returns an error message if the header is invalid.
    pub fn validate(&self) -> Result<(), String> {
        if self.tag != T::TAG {
            return Err(format!("{}: Invalid tag {:X}", T::NAME, self.tag));
        }
        for v in T::MIN_VERSION..T::VERSION + 1 {
            if self.version == v {
                if (self.flags & T::mask(v)) == self.flags {
                    return T::validate(self);
                } else {
                    return Err(format!("{}: Invalid flags {:X} for version {}", T::NAME, self.flags, self.version));
                }
            }
        }
        Err(format!("{}: Invalid version {} (expected {} to {})", T::NAME, self.version, T::MIN_VERSION, T::VERSION))
    }

    /// Returns `true` if the given file starts with a header of this type.
    pub fn found_in<P: AsRef<Path>>(filename: P) -> bool {
        if let Ok(mut file) = std::fs::File::open(filename)
            && let Ok(header) = Self::load(&mut file) {
            return header.tag == T::TAG;
        }
        false
    }

    /// Returns a reference to the payload.
    #[inline]
    pub fn payload(&self) -> &T {
        &self.payload
    }

    /// Returns a mutable reference to the payload.
    #[inline]
    pub fn payload_mut(&mut self) -> &mut T {
        &mut self.payload
    }
}

impl<T: Payload> Default for Header<T> {
    fn default() -> Self {
        Header {
            tag: T::TAG,
            version: T::VERSION,
            payload: T::default(),
            flags: T::DEFAULT_FLAGS,
        }
    }
}

impl<T: Payload> Serializable for Header<T> {}

//-----------------------------------------------------------------------------

/// Format-specific payload stored in a file format header.
///
/// The implementing type must be either empty or `#[repr(C)]`.
/// If not empty, the size must be a multiple of 8 bytes.
/// See [`Header`] for an example.
pub trait Payload: Copy + Eq + Default {
    /// User-friendly type name for the header.
    const NAME: &'static str;

    /// The first four bytes of the header as an unsigned little-endian integer.
    const TAG: u32;

    /// The latest supported version.
    const VERSION: u32;

    /// The earliest supported version.
    const MIN_VERSION: u32;

    /// Binary flags that should be set by default.
    const DEFAULT_FLAGS: u64;

    /// Updates the header to the latest version.
    fn update(&mut self);

    /// Returns the binary mask corresponding to valid flags in the specified version.
    fn mask(version: u32) -> u64;

    /// Performs type-specific validation and returns an error message if the header is invalid.
    fn validate(header: &Header<Self>) -> Result<(), String>;
}

//-----------------------------------------------------------------------------

/// Payload for the GBWT header.
#[repr(C)]
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct GBWTPayload {
    /// Number of sequences in the GBWT.
    pub sequences: usize,

    /// Total length of the sequences, including the endmarkers.
    pub size: usize,

    /// Alphabet offset: node identifiers in `1..offset + 1` are not used.
    pub offset: usize,

    /// Alphabet size: all node identifiers are in `1..alphabet_size`.
    pub alphabet_size: usize,
}

impl GBWTPayload {
    /// The GBWT index is bidirectional.
    pub const FLAG_BIDIRECTIONAL: u64 = 0x0001;

    /// The GBWT index contains metadata.
    pub const FLAG_METADATA: u64      = 0x0002;

    /// The serialized data is in the simple-sds format.
    pub const FLAG_SIMPLE_SDS: u64    = 0x0004;
}

impl Payload for GBWTPayload {
    const NAME: &'static str = "GBWTHeader";
    const TAG: u32 = 0x6B376B37;
    const VERSION: u32 = 5;
    const MIN_VERSION: u32 = 5;
    const DEFAULT_FLAGS: u64 = Self::FLAG_SIMPLE_SDS;

    fn update(&mut self) {}

    fn mask(_: u32) -> u64 {
        Self::FLAG_BIDIRECTIONAL | Self::FLAG_METADATA | Self::FLAG_SIMPLE_SDS
    }

    fn validate(header: &Header<Self>) -> Result<(), String> {
        if !header.is_set(Self::FLAG_SIMPLE_SDS) {
            return Err(format!("{}: SDSL format is not supported", Self::NAME));
        }
        Ok(())
    }
}

//-----------------------------------------------------------------------------

/// Payload for the GBWT metadata header.
#[repr(C)]
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct MetadataPayload {
    /// Number of samples in the GBWT.
    pub sample_count: usize,

    /// Number of haplotypes in the GBWT.
    pub haplotype_count: usize,

    /// Number of contigs in the graph.
    pub contig_count: usize,
}

impl MetadataPayload {
    /// The metadata contains path names.
    pub const FLAG_PATH_NAMES: u64   = 0x0001;

    /// The metadata contains sample names.
    pub const FLAG_SAMPLE_NAMES: u64 = 0x0002;

    /// The metadata contains contig names.
    pub const FLAG_CONTIG_NAMES: u64 = 0x0004;
}

impl Payload for MetadataPayload {
    const NAME: &'static str = "MetadataHeader";
    const TAG: u32 = 0x6B375E7A;
    const VERSION: u32 = 2;
    const MIN_VERSION: u32 = 2;
    const DEFAULT_FLAGS: u64 = 0;

    fn update(&mut self) {}

    fn mask(_: u32) -> u64 {
        Self::FLAG_PATH_NAMES | Self::FLAG_SAMPLE_NAMES | Self::FLAG_CONTIG_NAMES
    }

    fn validate(_: &Header<Self>) -> Result<(), String> {
        Ok(())
    }
}

//-----------------------------------------------------------------------------

/// Payload for the GBWTGraph header.
#[repr(C)]
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct GraphPayload {
    /// Number of nodes in the original graph.
    pub nodes: usize,
}

impl GraphPayload {
    /// The graph contains a node-to-segment translation.
    pub const FLAG_TRANSLATION: u64 = 0x0001;

    /// The serialized data is in the simple-sds format.
    pub const FLAG_SIMPLE_SDS: u64  = 0x0002;

    /// First version with zstd-compressed sequences.
    pub const ZSTD_VERSION: u32 = 4;
}

impl Payload for GraphPayload {
    const NAME: &'static str = "GraphHeader";
    const TAG: u32 = 0x6B3764AF;
    const VERSION: u32 = 4;
    const MIN_VERSION: u32 = 3;
    const DEFAULT_FLAGS: u64 = Self::FLAG_SIMPLE_SDS;

    fn update(&mut self) {}

    fn mask(_: u32) -> u64 {
        Self::FLAG_TRANSLATION | Self::FLAG_SIMPLE_SDS
    }

    fn validate(header: &Header<Self>) -> Result<(), String> {
        if !header.is_set(Self::FLAG_SIMPLE_SDS) {
            return Err(format!("{}: SDSL format is not supported", Self::NAME));
        }
        Ok(())
    }
}

//-----------------------------------------------------------------------------

/// Payload for the GBZ header.
#[repr(C)]
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct GBZPayload {
}

impl Payload for GBZPayload {
    const NAME: &'static str = "GBZHeader";
    const TAG: u32 = 0x205A4247;
    const VERSION: u32 = 2;
    const MIN_VERSION: u32 = 1;
    const DEFAULT_FLAGS: u64 = 0;

    fn update(&mut self) {}

    fn mask(_: u32) -> u64 {
        0
    }

    fn validate(_: &Header<Self>) -> Result<(), String> {
        Ok(())
    }
}

//-----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    use std::fs;

    use simple_sds::serialize;

    #[test]
    fn gbwt_header() {
        let header = Header::<GBWTPayload>::new();
        if let Err(msg) = header.validate() {
            panic!("{}", msg);
        }
        assert!(!header.is_set(GBWTPayload::FLAG_BIDIRECTIONAL), "Default: Bidirectional flag is set");
        assert!(!header.is_set(GBWTPayload::FLAG_METADATA), "Default: Metadata flag is set");
        assert!(header.is_set(GBWTPayload::FLAG_SIMPLE_SDS), "Default: Simple-SDS flag is not set");
        serialize::test(&header, "gbwt-header", Some(6), true);

        // We only have to test setting / unsetting flags for one payload type.
        let mut header = header;
        header.set(GBWTPayload::FLAG_BIDIRECTIONAL);
        header.set(GBWTPayload::FLAG_METADATA);
        if let Err(msg) = header.validate() {
            panic!("{}", msg);
        }
        assert!(header.is_set(GBWTPayload::FLAG_BIDIRECTIONAL), "Modified: Bidirectional flag could not be set");
        assert!(header.is_set(GBWTPayload::FLAG_METADATA), "Modified: Metadata flag could not be set");
        assert!(header.is_set(GBWTPayload::FLAG_SIMPLE_SDS), "Modified: Simple-SDS flag is not set");
        serialize::test(&header, "modified-gbwt-header", Some(6), true);

        header.unset(GBWTPayload::FLAG_METADATA);
        assert!(!header.is_set(GBWTPayload::FLAG_METADATA), "Modified: Metadata flag could not be unset");
    }

    #[test]
    fn metadata_header() {
        let header = Header::<MetadataPayload>::new();
        if let Err(msg) = header.validate() {
            panic!("{}", msg);
        }
        assert!(!header.is_set(MetadataPayload::FLAG_PATH_NAMES), "Default: Path name flag is set");
        assert!(!header.is_set(MetadataPayload::FLAG_SAMPLE_NAMES), "Default: Sample name flag is set");
        assert!(!header.is_set(MetadataPayload::FLAG_CONTIG_NAMES), "Default: Contig name flag is set");
        serialize::test(&header, "metadata-header", Some(5), true);
    }

    #[test]
    fn graph_header() {
        let header = Header::<GraphPayload>::new();
        if let Err(msg) = header.validate() {
            panic!("{}", msg);
        }
        assert!(!header.is_set(GraphPayload::FLAG_TRANSLATION), "Default: Translation flag is set");
        assert!(header.is_set(GraphPayload::FLAG_SIMPLE_SDS), "Default: Simple-SDS flag is not set");
        serialize::test(&header, "graph-header", Some(3), true);
    }

    #[test]
    fn gbz_header() {
        let header = Header::<GBZPayload>::new();
        if let Err(msg) = header.validate() {
            panic!("{}", msg);
        }
        serialize::test(&header, "gbz-header", Some(2), true);
    }

    #[test]
    fn found_in() {
        let header = Header::<GBZPayload>::new();
        let name = "found-in";
        let filename = serialize::temp_file_name(name);
        serialize::serialize_to(&header, &filename).unwrap();
        assert!(Header::<GBZPayload>::found_in(&filename), "The file does not start with a GBZ header");
        assert!(!Header::<GBWTPayload>::found_in(&filename), "The file starts with a GBWT header");
        fs::remove_file(&filename).unwrap();
        assert!(!Header::<GBZPayload>::found_in(&filename), "Deleted file starts with a GBZ header");
    }
}

//-----------------------------------------------------------------------------
