//! Support structures for GBWT and GBZ.

use crate::Pos;
use crate::bwt::SmallPos;

use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Vector, Access, Push, BitVec, Select};
use simple_sds::serialize::Serialize;
use simple_sds::sparse_vector::SparseVector;
use simple_sds::bits;

use zstd::stream::Encoder as ZstdEncoder;
use zstd::stream::Decoder as ZstdDecoder;

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, BTreeSet};
use std::collections::btree_map::Iter as TagIter;
use std::collections::hash_map::Entry;
use std::convert::TryFrom;
use std::io::{Error, ErrorKind, Write, Read};
use std::iter::FusedIterator;
use std::ops::Range;
use std::path::PathBuf;
use std::str::Utf8Error;
use std::{cmp, fmt, io, mem};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Orientation of a node or a path in a bidirected sequence graph.
#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Orientation {
    /// Forward orientation.
    Forward = 0,
    /// Reverse or reverse complement orientation.
    Reverse = 1,
}

impl Orientation {
    /// Returns the other orientation.
    #[inline]
    pub fn flip(&self) -> Orientation {
        match *self {
            Self::Forward => Self::Reverse,
            Self::Reverse => Self::Forward,
        }
    }
}

impl fmt::Display for Orientation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Orientation::Forward => write!(f, "forward"),
            Orientation::Reverse => write!(f, "reverse"),
        }
    }
}

/// Side of a node.
///
/// Some graph algorithms use combinations (node id, node side) as nodes.
#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum NodeSide {
    /// Left side.
    ///
    /// This is the entry side for [`Orientation::Forward`] and the exit side for [`Orientation::Reverse`].
    Left = 0,
    /// Right side.
    ///
    /// This is the exit side for [`Orientation::Forward`] and the entry side for [`Orientation::Reverse`].
    Right = 1,
}

impl NodeSide {
    /// Returns the other side.
    #[inline]
    pub fn flip(&self) -> NodeSide {
        match *self {
            Self::Left => Self::Right,
            Self::Right => Self::Left,
        }
    }
}

impl fmt::Display for NodeSide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NodeSide::Left => write!(f, "left"),
            NodeSide::Right => write!(f, "right"),
        }
    }
}

/// Position in a bidirected sequence graph.
#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct GraphPosition {
    /// Identifier of the node.
    pub node: usize,
    /// Orientation of the node.
    pub orientation: Orientation,
    /// Offset in the node.
    pub offset: usize,
}

impl GraphPosition {
    /// Creates a new graph position.
    #[inline]
    pub fn new(node: usize, orientation: Orientation, offset: usize) -> Self {
        GraphPosition {
            node, orientation, offset,
        }
    }

    /// Returns the GBWT node identifier corresponding to the position.
    #[inline]
    pub fn to_gbwt(&self) -> usize {
        encode_node(self.node, self.orientation)
    }
}

//-----------------------------------------------------------------------------

const fn generate_complement_table() -> [u8; 256] {
    let mut result: [u8; 256] = [b'N'; 256];
    result[b'A' as usize] = b'T'; result[b'a' as usize] = b'T';
    result[b'C' as usize] = b'G'; result[b'c' as usize] = b'G';
    result[b'G' as usize] = b'C'; result[b'g' as usize] = b'C';
    result[b'T' as usize] = b'A'; result[b't' as usize] = b'A';
    result
}

/// Complement table for DNA bases, normalized to upper case.
///
/// Invalid characters are mapped to `N`.
pub const COMPLEMENT: [u8; 256] = generate_complement_table();

/// Returns the reverse complement of the sequence, normalized to upper case.
///
/// Invalid characters are mapped to `N`.
pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::with_capacity(sequence.len());
    for &c in sequence.iter().rev() {
        result.push(COMPLEMENT[c as usize]);
    }
    result
}

//-----------------------------------------------------------------------------

/// A run as a (value, length) pair.
///
/// See [`SmallRun`] for a more compact version.
#[derive(Copy, Clone, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Run {
    /// Value in the run.
    pub value: usize,
    /// Length of the run.
    pub len: usize,
}

impl Run {
    /// Creates a new run.
    #[inline]
    pub fn new(value: usize, len: usize) -> Self {
        Run {
            value, len,
        }
    }
}

impl From<(usize, usize)> for Run {
    #[inline]
    fn from(run: (usize, usize)) -> Self {
        Self::new(run.0, run.1)
    }
}

/// A [`Run`] encoded using 32-bit integers to save space.
///
/// This is intended for situations, where we store a large number of runs explicitly.
/// The C++ implementation also uses 32-bit integers in similar situations.
#[derive(Copy, Clone, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct SmallRun {
    // Value in the run.
    pub value: u32,
    // Length of the run.
    pub len: u32,
}

impl SmallRun {
    /// Creates a new run.
    pub fn new(value: usize, len: usize) -> Self {
        SmallRun {
            value: value as u32,
            len: len as u32,

        }
    }
}

impl From<(usize, usize)> for SmallRun {
    fn from(run: (usize, usize)) -> Self {
        Self::new(run.0, run.1)
    }
}

impl From<SmallRun> for Run {
    fn from(run: SmallRun) -> Self {
        Run::new(run.value as usize, run.len as usize)
    }
}

impl From<Run> for SmallRun {
    fn from(run: Run) -> Self {
        Self::new(run.value, run.len)
    }
}

//-----------------------------------------------------------------------------

/// Returns the GBWT node identifier (handle) corresponding to the given original node and orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
///
/// # Arguments
///
/// * `id`: Identifier of the original node.
/// * `orientation`: Orientation of the node.
///
/// # Panics
///
/// May panic if `id > usize::MAX / 2`.
#[inline]
pub fn encode_node(id: usize, orientation: Orientation) -> usize {
    2 * id + (orientation as usize)
}

/// Returns the original node identifier corresponding to the given GBWT node (handle).
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn node_id(id: usize) -> usize {
    id / 2
}

/// Returns the orientation of the original node corresponding to the given GBWT node (handle).
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn node_orientation(id: usize) -> Orientation {
    match id & 1 {
        0 => Orientation::Forward,
        _ => Orientation::Reverse,
    }
}

/// Decodes a GBWT node identifier (handle) as a node identifier and orientation in the original graph.
#[inline]
pub fn decode_node(id: usize) -> (usize, Orientation) {
    (node_id(id), node_orientation(id))
}

/// Returns the GBWT node identifier (handle) for the same original node in the other orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn flip_node(id: usize) -> usize {
    id ^ 1
}

/// Returns the encoded node side corresponding to the given original node and node side.
///
/// # Arguments
///
/// * `id`: Identifier of the original node.
/// * `side`: Side of the node.
///
/// # Panics
///
/// May panic if `id > usize::MAX / 2`.
#[inline]
pub fn encode_node_side(id: usize, side: NodeSide) -> usize {
    2 * id + (side as usize)
}

/// Returns the original node identifier corresponding to the given encoded node side.
#[inline]
pub fn node_side_id(id: usize) -> usize {
    id / 2
}

/// Returns the node side corresponding to the given encoded node side.
#[inline]
pub fn node_side(id: usize) -> NodeSide {
    match id & 1 {
        0 => NodeSide::Left,
        _ => NodeSide::Right,
    }
}

/// Decodes an encoded node side as a node identifier and node side in the original graph.
#[inline]
pub fn decode_node_side(id: usize) -> (usize, NodeSide) {
    (node_side_id(id), node_side(id))
}

/// Returns the encoded node side for the other side of the same original node.
#[inline]
pub fn flip_node_side(id: usize) -> usize {
    id ^ 1
}

/// Returns the entry side of the given orientation.
#[inline]
pub fn entry_side(orientation: Orientation) -> NodeSide {
    match orientation {
        Orientation::Forward => NodeSide::Left,
        Orientation::Reverse => NodeSide::Right,
    }
}

/// Returns the exit side of the given orientation.
#[inline]
pub fn exit_side(orientation: Orientation) -> NodeSide {
    match orientation {
        Orientation::Forward => NodeSide::Right,
        Orientation::Reverse => NodeSide::Left,
    }
}

/// Returns the orientation corresponding to the given entry side.
#[inline]
pub fn entry_orientation(side: NodeSide) -> Orientation {
    match side {
        NodeSide::Left => Orientation::Forward,
        NodeSide::Right => Orientation::Reverse,
    }
}

/// Returns the orientation corresponding to the given exit side.
#[inline]
pub fn exit_orientation(side: NodeSide) -> Orientation {
    match side {
        NodeSide::Left => Orientation::Reverse,
        NodeSide::Right => Orientation::Forward,
    }
}

//-----------------------------------------------------------------------------

/// Returns `true` if the given edge is in canonical orientation.
///
/// An edge is canonical, if the destination node (`to`) has a higher identifier than the source node (`from`).
/// A self-loop with at least one node in forward orientation is also canonical.
pub fn edge_is_canonical(from: (usize, Orientation), to: (usize, Orientation)) -> bool {
    if from.1 == Orientation::Forward {
        to.0 >= from.0
    } else {
        (to.0 > from.0) || (to.0 == from.0 && to.1 == Orientation::Forward)
    }
}

/// Returns `true` if the given edge is in canonical orientation.
///
/// This version takes the edge as a pair of GBWT node identifiers.
/// See [`edge_is_canonical`] for more information on canonicality and [`crate::GBWT`] on GBWT node identifiers.
pub fn encoded_edge_is_canonical(from: usize, to: usize) -> bool {
    edge_is_canonical(decode_node(from), decode_node(to))
}

//-----------------------------------------------------------------------------

/// Returns the sequence identifier corresponding to the given path and orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
///
/// # Arguments
///
/// * `id`: Identifier of the path.
/// * `orientation`: Orientation of the path.
///
/// # Panics
///
/// May panic if `id > usize::MAX / 2`.
#[inline]
pub fn encode_path(id: usize, orientation: Orientation) -> usize {
    2 * id + (orientation as usize)
}

/// Returns the path identifier corresponding to the given sequence.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn path_id(id: usize) -> usize {
    id / 2
}

/// Returns the orientation of the path corresponding to the given sequence.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn path_orientation(id: usize) -> Orientation {
    match id & 1 {
        0 => Orientation::Forward,
        _ => Orientation::Reverse,
    }
}

/// Decodes a sequence identifier as a path identifier and orientation in the original graph.
#[inline]
pub fn decode_path(id: usize) -> (usize, Orientation) {
    (path_id(id), path_orientation(id))
}

/// Returns the sequence identifier for the same path in the other orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn flip_path(id: usize) -> usize {
    id ^ 1
}

/// Returns `true` if the path is in canonical orientation.
///
/// A path is canonical, if:
///
/// 1. It is empty.
/// 2. Both the first and the last node are in forward orientation.
/// 3. The edge from the first node to the last node would be canonical.
pub fn path_is_canonical(path: &[(usize, Orientation)]) -> bool {
    if path.is_empty() {
        return true;
    }

    let first = path[0];
    let last = path[path.len() - 1];
    if first.1 == last.1 {
        return first.1 == Orientation::Forward;
    }

    edge_is_canonical(first, last)
}

/// Returns `true` if the path is in canonical orientation.
///
/// The path is assumed to be a sequence of GBWT node identifiers.
/// See [`path_is_canonical`] for more information on canonicality and [`crate::GBWT`] on GBWT node identifiers.
pub fn encoded_path_is_canonical(path: &[usize]) -> bool {
    if path.is_empty() {
        return true;
    }

    let first = decode_node(path[0]);
    let last = decode_node(path[path.len() - 1]);
    if first.1 == last.1 {
        return first.1 == Orientation::Forward;
    }

    edge_is_canonical(first, last)
}

/// Returns the reverse of the given path.
///
/// The path is assumed to be a sequence of GBWT node identifiers.
/// A reverse path visits the other orientation of each node in reverse order.
/// See also [`reverse_path_in_place`] and [`reverse_path_into`].
pub fn reverse_path(path: &[usize]) -> Vec<usize> {
    let mut result: Vec<usize> = path.iter().map(|x| flip_node(*x)).collect();
    result.reverse();
    result
}

/// Reverses the given path in place.
///
/// The path is assumed to be a sequence of GBWT node identifiers.
/// A reverse path visits the other orientation of each node in reverse order.
/// See also [`reverse_path`] and [`reverse_path_into`].
pub fn reverse_path_in_place(path: &mut [usize]) {
    path.reverse();
    for node in path.iter_mut() {
        *node = flip_node(*node);
    }
}

/// Inserts the reversed path into the given buffer.
///
/// The buffer uses 32-bit node identifiers to save space during GBWT construction.
/// The path is assumed to be a sequence of GBWT node identifiers.
/// A reverse path visits the other orientation of each node in reverse order.
/// See also [`reverse_path`] and [`reverse_path_in_place`].
pub fn reverse_path_into(path: &[usize], buffer: &mut Vec<u32>) {
    for node in path.iter().rev() {
        buffer.push(flip_node(*node) as u32);
    }
}

//-----------------------------------------------------------------------------

/// Returns the intersection of two ranges.
#[inline]
pub fn intersect(a: &Range<usize>, b: &Range<usize>) -> Range<usize> {
    cmp::max(a.start, b.start)..cmp::min(a.end, b.end)
}

//-----------------------------------------------------------------------------

/// Returns the full file name for a specific test file.
pub fn get_test_data(filename: &'static str) -> PathBuf {
    let mut buf = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    buf.push("test-data");
    buf.push(filename);
    buf
}

//-----------------------------------------------------------------------------

/// An immutable array of immutable strings.
///
/// The strings are concatenated and stored in a single byte vector.
/// This reduces the space overhead for the strings and the time overhead for serializing and loading them.
/// Serialization with [`Serialize`] further compresses the starting positions and compacts the alphabet in an attempt to use fewer than 8 bits per byte.
/// Compressed serialization using [`Self::compress`] and [`Self::decompress`] uses Zstandard compression for the concatenated strings.
///
/// `StringArray` can be built from a [`Vec`] or a slice of any type that can be converted to a string slice.
/// Construction from an iterator is not feasible, as `StringArray` needs to know the total length of the strings in advance.
///
/// Because the bytes may come from an untrusted source, `StringArray` does not assume that the bytes are valid UTF-8 strings.
///
/// # Examples
///
/// ```
/// use gbz::support::StringArray;
///
/// let source = vec!["first", "second", "third", "fourth"];
/// let array = StringArray::from(source.as_slice());
/// assert_eq!(array.len(), source.len());
/// for i in 0..array.len() {
///     assert_eq!(array.str(i).unwrap(), source[i]);
/// }
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StringArray {
    index: IntVector,
    strings: Vec<u8>,
}

impl StringArray {
    /// Default compression level for Zstandard.
    pub const DEFAULT_COMPRESSION_LEVEL: i32 = 3;

    /// Returns the number of strings in the array.
    #[inline]
    pub fn len(&self) -> usize {
        self.index.len() - 1
    }

    /// Returns `true` if the array is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the length of the `i`th string in bytes.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn str_len(&self, i: usize) -> usize {
        (self.index.get(i + 1) - self.index.get(i)) as usize
    }

    /// Returns the total length of a range of strings in bytes.
    ///
    /// # Panics
    ///
    /// May panic if `strings.start > strings.end` or `strings.end > self.len()`.
    pub fn range_len(&self, strings: Range<usize>) -> usize {
        (self.index.get(strings.end) - self.index.get(strings.start)) as usize
    }

    /// Returns a byte slice corresponding to the `i`th string.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn bytes(&self, i: usize) -> &[u8] {
        let start = self.index.get(i) as usize;
        let limit = self.index.get(i + 1) as usize;
        &self.strings[start..limit]
    }

    /// Returns a byte slice corresponding to a range of strings.
    ///
    /// # Panics
    ///
    /// May panic if `strings.start > strings.end` or `strings.end > self.len()`.
    pub fn range(&self, strings: Range<usize>) -> &[u8] {
        let start = self.index.get(strings.start) as usize;
        let limit = self.index.get(strings.end) as usize;
        &self.strings[start..limit]
    }

    /// Returns a string slice corresponding to the `i`th string or an error if the bytes are not valid UTF-8.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn str(&self, i: usize) -> Result<&str, Utf8Error> {
        std::str::from_utf8(self.bytes(i))
    }

    /// Returns a copy of the `i`th string or an error if the bytes are not valid UTF-8.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn string(&self, i: usize) -> Result<String, Utf8Error> {
        match self.str(i) {
            Ok(v) => Ok(v.to_string()),
            Err(e) => Err(e),
        }
    }

    /// Returns an iterator over the string array.
    pub fn iter(&self) -> StringIter<'_> {
        StringIter {
            parent: self,
            next: 0,
            limit: self.len(),
        }
    }

    // Builds an empty string array with capacity for `n` strings of total length `total_len`.
    fn with_capacity(n: usize, total_len: usize) -> StringArray {
        let mut index = IntVector::with_capacity(n + 1, bits::bit_len(total_len as u64)).unwrap();
        index.push(0);
        let strings: Vec<u8> = Vec::with_capacity(total_len);
        StringArray {
            index, strings,
        }
    }

    // Appends a new string to the array, assuming that there is space for it.
    fn append(&mut self, string: &str) {
        self.strings.extend(string.bytes());
        self.index.push(self.strings.len() as u64);
    }

    // Returns (bytes to packed, packed to bytes, packed character width).
    fn alphabet(data: &[u8]) -> (Vec<usize>, Vec<u8>, usize) {
        // Determine the byte values that are present.
        let mut bytes_to_packed: Vec<usize> = vec![0; 1 << 8];
        for byte in data {
            bytes_to_packed[*byte as usize] = 1;
        }

        // Determine alphabet size.
        let sigma = bytes_to_packed.iter().sum();
        let width = bits::bit_len(cmp::max(sigma, 1) as u64 - 1);

        // Build the alphabet mappings.
        let mut packed_to_bytes: Vec<u8> = vec![0; sigma];
        let mut rank = 0;
        for (index, value) in bytes_to_packed.iter_mut().enumerate() {
            if *value != 0 {
                *value = rank;
                packed_to_bytes[rank] = index as u8;
                rank += 1;
            }
        }

        (bytes_to_packed, packed_to_bytes, width)
    }

    /// Serializes the struct as a compressed string array to the given writer.
    ///
    /// Uses Zstandard compression for the sequences.
    /// If a compression level is not provided, [`Self::DEFAULT_COMPRESSION_LEVEL`] is used.
    /// See [`Self::serialize`] for a non-compressed serialization format.
    ///
    /// # Errors
    ///
    /// Passes through any I/O errors and compression errors.
    pub fn compress<T: io::Write>(&self, writer: &mut T, compression_level: Option<i32>) -> io::Result<()> {
        // Compress the index without the past-the-end sentinel.
        let sv = SparseVector::try_from_iter(self.index.iter().take(self.len()).map(|x| x as usize)).unwrap();
        sv.serialize(writer)?;
        drop(sv);

        // We want to know the uncompressed length of the concatenated strings.
        let total_len = self.strings.len();
        total_len.serialize(writer)?;

        // Compress the strings into a vector of bytes using zstd and serialize it.
        // We cannot write directly into the writer, as we need to know the compressed size in advance.
        let compression_level = compression_level.unwrap_or(Self::DEFAULT_COMPRESSION_LEVEL);
        let mut encoder = ZstdEncoder::new(Vec::new(), compression_level)?;
        encoder.write_all(&self.strings)?;
        let compressed = encoder.finish()?;
        compressed.serialize(writer)?;

        Ok(())
    }

    /// Deserializes a compressed string array from the given reader.
    ///
    /// Uses Zstandard compression for the sequences.
    /// See [`Self::load`] for a non-compressed serialization format.
    ///
    /// # Errors
    ///
    /// Passes through any I/O errors and decompression errors.
    /// Returns [`ErrorKind::InvalidData`] if the data is not internally consistent.
    pub fn decompress<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        // Load the compressed index and the total length of the concatenated strings.
        let sv = SparseVector::load(reader)?;
        let total_len = usize::load(reader)?;

        // Decompress the index.
        let mut index = IntVector::with_capacity(sv.count_ones() + 1, bits::bit_len(total_len as u64)).unwrap();
        index.extend(sv.one_iter().map(|(_, x)| x));
        index.push(total_len as u64);
        drop(sv);
        if index.get(0) != 0 {
            return Err(Error::new(ErrorKind::InvalidData, "StringArray: First string does not start at offset 0"));
        }

        // Load and decompress the strings.
        let compressed: Vec<u8> = Vec::load(reader)?;
        let mut decoder = ZstdDecoder::new(&compressed[..])?;
        let mut strings: Vec<u8> = Vec::with_capacity(total_len);
        decoder.read_to_end(&mut strings)?;
        if strings.len() != total_len {
            return Err(Error::new(ErrorKind::InvalidData, "StringArray: Decompressed string length does not match the expected length"));
        }

        Ok(StringArray {
            index, strings,
        })
    }

    /// Returns the size of the compressed struct in [`u64`] elements with the given compression level.
    ///
    /// If a compression level is not provided, [`Self::DEFAULT_COMPRESSION_LEVEL`] is used.
    /// See [`Self::size_in_elements`] for the size of the non-compressed struct.
    ///
    /// # Panics
    ///
    /// May panic due to compression errors.
    pub fn compressed_size_in_elements(&self, compression_level: Option<i32>) -> usize {
        let mut result = 0;

        let sv = SparseVector::try_from_iter(self.index.iter().take(self.len()).map(|x| x as usize)).unwrap();
        result += sv.size_in_elements();
        drop(sv);

        let total_len = self.strings.len();
        result += total_len.size_in_elements();

        let compression_level = compression_level.unwrap_or(Self::DEFAULT_COMPRESSION_LEVEL);
        let mut encoder = ZstdEncoder::new(Vec::new(), compression_level).unwrap();
        encoder.write_all(&self.strings).unwrap();
        let compressed = encoder.finish().unwrap();
        result += compressed.size_in_elements();

        result
    }
}

impl Serialize for StringArray {
    fn serialize_header<T: io::Write>(&self, _: &mut T) -> io::Result<()> {
        Ok(())
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        // Compress the index without the past-the-end sentinel.
        let sv = SparseVector::try_from_iter(self.index.iter().take(self.len()).map(|x| x as usize)).unwrap();
        sv.serialize(writer)?;
        drop(sv);

        // Determine and serialize the alphabet
        let (pack, alphabet, width) = Self::alphabet(&self.strings);
        alphabet.serialize(writer)?;

        // Pack and serialize the strings.
        let mut packed = IntVector::new(width).unwrap();
        packed.extend(self.strings.iter().map(|x| pack[*x as usize]));
        packed.serialize(writer)?;

        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        // Load the compressed index. We need the strings for the past-the-end sentinel.
        let sv = SparseVector::load(reader)?;

        // Load the alphabet.
        let alphabet = Vec::<u8>::load(reader)?;

        // Load and decompress the strings.
        let packed = IntVector::load(reader)?;
        let strings: Vec<u8> = packed.into_iter().map(|x| alphabet[x as usize]).collect();

        // Decompress the index.
        let mut index = IntVector::with_capacity(sv.count_ones() + 1, bits::bit_len(strings.len() as u64)).unwrap();
        index.extend(sv.one_iter().map(|(_, x)| x));
        index.push(strings.len() as u64);

        // Sanity checks.
        if index.get(0) != 0 {
            return Err(Error::new(ErrorKind::InvalidData, "StringArray: First string does not start at offset 0"));
        }
        Ok(StringArray {
            index, strings,
        })
    }

    fn size_in_elements(&self) -> usize {
        let sv = SparseVector::try_from_iter(self.index.iter().take(self.len()).map(|x| x as usize)).unwrap();
        let (_, alphabet, width) = Self::alphabet(&self.strings);

        sv.size_in_elements() + alphabet.size_in_elements() + IntVector::size_by_params(self.strings.len(), width)
    }
}

impl<T: AsRef<str>> From<&[T]> for StringArray {
    fn from(v: &[T]) -> Self {
        let total_len = v.iter().fold(0, |sum, item| sum + item.as_ref().len());
        let mut result = StringArray::with_capacity(v.len(), total_len);
        for string in v.iter() {
            result.append(string.as_ref());
        }
        result
    }
}
impl<T: AsRef<str>> From<Vec<T>> for StringArray {
    fn from(v: Vec<T>) -> Self {
        StringArray::from(v.as_slice())
    }
}

//-----------------------------------------------------------------------------

/// A read-only iterator over [`StringArray`].
///
/// The type of `Item` is `&[`[`u8`]`]`.
///
/// # Examples
///
/// ```
/// use gbz::support::StringArray;
/// use std::str;
///
/// let source = vec!["first", "second", "third"];
/// let array = StringArray::from(source.as_slice());
/// for (index, bytes) in array.iter().enumerate() {
///     assert_eq!(bytes, source[index].as_bytes());
/// }
/// ```
#[derive(Clone, Debug)]
pub struct StringIter<'a> {
    parent: &'a StringArray,
    // The first index we have not used.
    next: usize,
    // The first index we should not use.
    limit: usize,
}

impl<'a> Iterator for StringIter<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.limit {
            None
        } else {
            let result = Some(self.parent.bytes(self.next));
            self.next += 1;
            result
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.limit - self.next;
        (remaining, Some(remaining))
    }
}

impl<'a> DoubleEndedIterator for StringIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.next >= self.limit {
            None
        } else {
            self.limit -= 1;
            Some(self.parent.bytes(self.limit))
        }
    }
}

impl<'a> ExactSizeIterator for StringIter<'a> {}

impl<'a> FusedIterator for StringIter<'a> {}

//-----------------------------------------------------------------------------

/// An immutable set of immutable strings with integer identifiers.
///
/// The strings are stored in a [`StringArray`] and the identifiers are indexes into the array.
///
/// A `Dictionary` can be built from a [`StringArray`] or a [`Vec`] or a slice of any type that can be converted to a string slice.
/// The construction will fail if the source contains duplicate strings.
///
/// # Examples
///
/// ```
/// use gbz::support::Dictionary;
/// use std::convert::TryFrom;
///
/// let source = vec!["first", "second", "third", "fourth"];
/// let dict = Dictionary::try_from(source.as_slice()).unwrap();
/// for (index, value) in source.iter().enumerate() {
///     assert_eq!(dict.id(value), Some(index));
///     assert_eq!(dict.bytes(index), source[index].as_bytes());
/// }
/// assert_eq!(dict.id("fifth"), None);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Dictionary {
    strings: StringArray,
    sorted_ids: IntVector,
}

impl Dictionary {
    /// Returns the number of strings in the dictionary.
    #[inline]
    pub fn len(&self) -> usize {
        self.strings.len()
    }

    /// Returns `true` if the dictionary is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the identifier of the given string in the dictionary, or [`None`] if there is no such string.
    pub fn id<T: AsRef<[u8]>>(&self, string: T) -> Option<usize> {
        let mut low = 0;
        let mut high = self.len();
        while low < high {
            let mid = low + (high - low) / 2;
            let id = self.sorted_ids.get(mid) as usize;
            match string.as_ref().cmp(self.bytes(id)) {
                Ordering::Less => high = mid,
                Ordering::Equal => return Some(id),
                Ordering::Greater => low = mid + 1,
            }
        }
        None
    }

    /// Returns a byte slice corresponding to the string with identifier `i`.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn bytes(&self, i: usize) -> &[u8] {
        self.strings.bytes(i)
    }

    /// Returns a string slice corresponding to the string with identifier `i` or an error if the bytes are not valid UTF-8.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn str(&self, i: usize) -> Result<&str, Utf8Error> {
        self.strings.str(i)
    }

    /// Returns a copy of the string with identifier `i` or an error if the bytes are not valid UTF-8.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn string(&self, i: usize) -> Result<String, Utf8Error> {
        self.strings.string(i)
    }
}

impl Serialize for Dictionary {
    fn serialize_header<T: io::Write>(&self, _: &mut T) -> io::Result<()> {
        Ok(())
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.strings.serialize(writer)?;
        self.sorted_ids.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let strings = StringArray::load(reader)?;
        let sorted_ids = IntVector::load(reader)?;
        Ok(Dictionary {
            strings, sorted_ids,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.strings.size_in_elements() + self.sorted_ids.size_in_elements()
    }
}

impl TryFrom<StringArray> for Dictionary {
    type Error = String;

    fn try_from(source: StringArray) -> Result<Self, Self::Error> {
        // Sort the ids and check for duplicates.
        let mut sorted: Vec<usize> = Vec::with_capacity(source.len());
        for i in 0..source.len() {
            sorted.push(i);
        }
        sorted.sort_unstable_by(|a, b| source.bytes(*a).cmp(source.bytes(*b)));
        for i in 1..sorted.len() {
            if source.bytes(sorted[i - 1]) == source.bytes(sorted[i]) {
                return Err(String::from("Cannot build a dictionary from a source with duplicate strings"));
            }
        }

        // Compact the sorted ids.
        let width = if sorted.is_empty() { 1 } else { bits::bit_len(sorted.len() as u64 - 1) };
        let mut sorted_ids = IntVector::with_capacity(sorted.len(), width).unwrap();
        sorted_ids.extend(sorted);

        Ok(Dictionary {
            strings: source,
            sorted_ids,
        })
    }
}

impl<T: AsRef<str>> TryFrom<&[T]> for Dictionary {
    type Error = String;

    fn try_from(source: &[T]) -> Result<Self, Self::Error> {
        Self::try_from(StringArray::from(source))
    }
}

impl<T: AsRef<str>> TryFrom<Vec<T>> for Dictionary {
    type Error = String;

    fn try_from(source: Vec<T>) -> Result<Self, Self::Error> {
        Self::try_from(StringArray::from(source))
    }
}

impl AsRef<StringArray> for Dictionary {
    #[inline]
    fn as_ref(&self) -> &StringArray {
        &(self.strings)
    }
}

//-----------------------------------------------------------------------------

/// A key-value structure with strings as both keys and values.
///
/// The keys are case insensitive.
/// This structure is a simple wrapper over [`BTreeMap`]`<`[`String`]`, `[`String`]`>` that converts all keys to lower case.
///
/// # Examples
///
/// ```
/// use gbz::Tags;
///
/// let mut tags = Tags::new();
/// tags.insert("first-key", "first-value");
/// tags.insert("second-key", "second-value");
/// assert!(tags.contains_key("First-Key"));
/// assert_eq!(*tags.get("second-key").unwrap(), "second-value");
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Tags {
    tags: BTreeMap<String, String>,
}

impl Tags {
    /// Creates an empty `Tags` structure.
    pub fn new() -> Tags {
        Tags::default()
    }

    /// Returns the number of tags.
    pub fn len(&self) -> usize {
        self.tags.len()
    }

    /// Returns `true` if the structure is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the value corresponding to the key, or [`None`] no such tag exists.
    pub fn get(&self, key: &str) -> Option<&String> {
        let key = key.to_lowercase();
        self.tags.get(&key)
    }

    /// Returns `true` if there is a tag with the given key.
    pub fn contains_key(&self, key: &str) -> bool {
        let key = key.to_lowercase();
        self.tags.contains_key(&key)
    }

    /// Inserts a new tag, overwriting the possible old value associated with the same key.
    ///
    /// # Arguments
    ///
    /// * `key`: Key of the tag. The key is converted to lower case before it is inserted into the hash table.
    /// * `value`: Value of the tag.
    pub fn insert(&mut self, key: &str, value: &str) {
        let key = key.to_lowercase();
        let _ = self.tags.insert(key, value.to_string());
    }

    /// Removes a tag with the given key and returns its value, or [`None`] if there is no such tag.
    pub fn remove(&mut self, key: &str) -> Option<String> {
        let key = key.to_lowercase();
        self.tags.remove(&key)
    }

    /// Returns an iterator that visits all tags in sorted order by keys.
    ///
    /// The type of `Item` is `(&`[`String`]`, &`[`String`]`)`.
    pub fn iter(&self) -> TagIter<'_, String, String> {
        self.tags.iter()
    }

    // Returns the array of keys and values in serialized order.
    fn linearize(&self) -> StringArray {
        let mut linearized: Vec<&str> = Vec::with_capacity(2 * self.len());
        for (key, value) in self.iter() {
            linearized.push(key); linearized.push(value);
        }
        StringArray::from(linearized)
    }
}

impl Serialize for Tags {
    fn serialize_header<T: io::Write>(&self, _: &mut T) -> io::Result<()> {
        Ok(())
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        let linearized = self.linearize();
        linearized.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let linearized = StringArray::load(reader)?;
        if linearized.len() % 2 != 0 {
            return Err(Error::new(ErrorKind::InvalidData, "Tags: Key without a value"));
        }
        let mut result = Tags::new();
        for i in 0..linearized.len() / 2 {
            let key = linearized.str(2 * i).map_err(|_| Error::new(ErrorKind::InvalidData, "Tags: Invalid UTF-8 in a key"))?;
            let value = linearized.str(2 * i + 1).map_err(|_| Error::new(ErrorKind::InvalidData, "Tags: Invalid UTF-8 in a value"))?;
            result.insert(key, value);
        }
        if result.len() != linearized.len() / 2 {
            return Err(Error::new(ErrorKind::InvalidData, "Tags: Duplicate keys"));
        }
        Ok(result)
    }

    fn size_in_elements(&self) -> usize {
        let linearized = self.linearize();
        linearized.size_in_elements()
    }
}

impl AsRef<BTreeMap<String, String>> for Tags {
    #[inline]
    fn as_ref(&self) -> &BTreeMap<String, String> {
        &(self.tags)
    }
}

impl From<ByteCode> for Vec<u8> {
    fn from(source: ByteCode) -> Self {
        source.bytes
    }
}

//-----------------------------------------------------------------------------

/// A variable-length encoder for unsigned integers.
///
/// `ByteCode` encodes an integer as a sequence of bytes in little-endian order and stores it in the internal [`Vec`].
/// Each byte contains 7 bits of data, and the high bit indicates whether the encoding continues.
/// The bytes can be accessed with [`AsRef`] or extracted with [`From`], and [`ByteCodeIter`] can be used for decoding the integers.
/// Raw bytes can be appended to the encoding using [`ByteCode::write_byte`].
///
/// # Examples
///
/// ```
/// use gbz::support::ByteCode;
///
/// let mut encoder = ByteCode::new();
/// encoder.write(123); encoder.write(456); encoder.write(789);
/// let bytes = encoder.as_ref();
/// assert_eq!(*bytes, [123, 72 + 128, 3, 21 + 128, 6]);
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ByteCode {
    bytes: Vec<u8>
}

impl ByteCode {
    const MASK: u8 = 0x7F;
    const FLAG: u8 = 0x80;
    const SHIFT: usize = 7;

    /// Creates a new encoder.
    pub fn new() -> Self {
        ByteCode::default()
    }

    /// Encodes `value` and stores the encoding.
    pub fn write(&mut self, value: usize) {
        let mut value = value;
        while value > (Self::MASK as usize) {
            self.bytes.push(((value as u8) & Self::MASK) | Self::FLAG);
            value >>= Self::SHIFT;
        }
        self.bytes.push(value as u8);
    }

    /// Appends a byte to the encoding.
    pub fn write_byte(&mut self, byte: u8) {
        self.bytes.push(byte);
    }

    /// Returns the total number of bytes in the encoding.
    #[inline]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Returns `true` if the encoding is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl AsRef<[u8]> for ByteCode {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.bytes
    }
}

/// An iterator that decodes integers from a byte slice encoded by [`ByteCode`].
///
/// The type of `Item` is [`usize`].
/// Raw bytes can be read from the encoding using [`ByteCodeIter::byte`].
///
/// # Examples
///
/// ```
/// use gbz::support::{ByteCode, ByteCodeIter};
///
/// let mut source = ByteCode::new();
/// source.write(123); source.write(456); source.write(789);
///
/// let mut iter = ByteCodeIter::new(source.as_ref());
/// assert_eq!(iter.next(), Some(123));
/// assert_eq!(iter.next(), Some(456));
/// assert_eq!(iter.next(), Some(789));
/// assert_eq!(iter.next(), None);
/// ```
#[derive(Clone, Debug)]
pub struct ByteCodeIter<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> ByteCodeIter<'a> {
    /// Returns an iterator over the byte slice.
    pub fn new(bytes: &'a [u8]) -> Self {
        ByteCodeIter {
            bytes,
            offset: 0,
        }
    }

    /// Returns the next byte from the slice, or [`None`] if there are no more bytes left.
    pub fn byte(&mut self) -> Option<u8> {
        if self.offset >= self.bytes.len() {
            return None;
        }
        let result = Some(self.bytes[self.offset]);
        self.offset += 1;
        result
    }

    /// Returns the first unvisited offset in the byte slice.
    #[inline]
    pub fn offset(&self) -> usize {
        self.offset
    }
}

impl<'a> Iterator for ByteCodeIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let mut offset = 0;
        let mut result = 0;
        while self.offset < self.bytes.len() {
            let value = unsafe { *self.bytes.get_unchecked(self.offset) };
            self.offset += 1;
            result += ((value & ByteCode::MASK) as usize) << offset;
            offset += ByteCode::SHIFT;
            if value & ByteCode::FLAG == 0 {
                return Some(result);
            }
        }
        None
    }
}

impl<'a> FusedIterator for ByteCodeIter<'a> {}

//-----------------------------------------------------------------------------

/// A run-length encoder for non-empty runs of unsigned integers.
///
/// The exact encoding depends on alphabet size `sigma`.
/// If `sigma` is small, the encoder tries to encode short runs as a single byte.
/// For long runs, the remaining run length is encoded using [`ByteCode`].
/// For a large `sigma`, both the value and the run length are encoded using [`ByteCode`].
/// Alphabet size `sigma == 0` indicates a large alphabet of unknown size.
///
/// The bytes can be accessed with [`AsRef`] or extracted with [`From`], and [`RLEIter`] can be used for decoding the integers.
/// Raw bytes and [`ByteCode`]-encoded integers can be appended to the encoding using [`RLE::write_byte`] and [`RLE::write_int`].
/// The following functions can for creating a byte stream with various encodings:
///
/// # Examples
///
/// ```
/// use gbz::support::{Run, RLE};
///
/// let mut encoder = RLE::with_sigma(4);
/// encoder.write(Run::new(3, 12)); encoder.write(Run::new(2, 721)); encoder.write(Run::new(0, 34));
/// assert_eq!(*encoder.as_ref(), [3 + 4 * 11, 2 + 4 * 63, 17 + 128, 5, 0 + 4 * 33]);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RLE {
    bytes: ByteCode,
    sigma: usize,
    threshold: usize,
}

impl RLE {
    const THRESHOLD: usize = 255;
    const UNIVERSE: usize = 256;

    /// Creates a new encoder with alphabet size `0`.
    pub fn new() -> Self {
        RLE::default()
    }

    /// Creates a new encoder with the given alphabet size.
    pub fn with_sigma(sigma: usize) -> Self {
        let (sigma, threshold) = Self::sanitize(sigma);
        RLE {
            bytes: ByteCode::new(),
            sigma,
            threshold,
        }
    }

    /// Encodes and stores a run.
    ///
    /// Does nothing if `run.len == 0`.
    ///
    /// # Panics
    ///
    /// Panics if `run.value >= self.sigma()`.
    pub fn write(&mut self, run: Run) {
        if run.len == 0 {
            return;
        }
        assert!(run.value < self.sigma(), "RLE: Cannot encode value {} with alphabet size {}", run.value, self.sigma);
        unsafe { self.write_unchecked(run); }
    }

    /// Encodes and stores a run.
    ///
    /// # Safety
    ///
    /// Behavior is undefined if `run.len == 0` or `run.value >= self.sigma()`.
    pub unsafe fn write_unchecked(&mut self, run: Run) {
        if self.sigma >= Self::THRESHOLD {
            self.bytes.write(run.value);
            self.bytes.write(run.len - 1);
        } else if run.len < self.threshold {
            self.write_basic(run.value, run.len);
        } else {
            self.write_basic(run.value, self.threshold);
            self.bytes.write(run.len - self.threshold);
        }
    }

    /// Appends a byte to the encoding.
    pub fn write_byte(&mut self, byte: u8) {
        self.bytes.write_byte(byte);
    }

    /// Encodes `value` using [`ByteCode`] and stores the encoding.
    pub fn write_int(&mut self, value: usize) {
        self.bytes.write(value);
    }

    /// Returns the total number of bytes in the encoding.
    #[inline]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Returns `true` if the encoding is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the alphabet size.
    #[inline]
    pub fn sigma(&self) -> usize {
        self.sigma
    }

    /// Changes the alphabet size to `sigma`.
    pub fn set_sigma(&mut self, sigma: usize) {
        let (sigma, threshold) = Self::sanitize(sigma);
        self.sigma = sigma;
        self.threshold = threshold;
    }

    // Writes a single-byte run.
    fn write_basic(&mut self, value: usize, len: usize) {
        let code = value + self.sigma * (len - 1);
        self.bytes.write_byte(code as u8);
    }

    // Returns (effective sigma, threshold for short runs).
    fn sanitize(sigma: usize) -> (usize, usize) {
        let sigma = if sigma == 0 { usize::MAX } else { sigma };
        let threshold = if sigma < Self::THRESHOLD { Self::UNIVERSE / sigma } else { 0 };
        (sigma, threshold)
    }
}

impl Default for RLE {
    fn default() -> Self {
        let (sigma, threshold) = Self::sanitize(0);
        RLE {
            bytes: ByteCode::new(),
            sigma,
            threshold,
        }
    }
}

impl AsRef<[u8]> for RLE {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        self.bytes.as_ref()
    }
}

impl From<RLE> for Vec<u8> {
    fn from(source: RLE) -> Self {
        Self::from(source.bytes)
    }
}

//-----------------------------------------------------------------------------

/// An iterator that decodes runs from a byte slice encoded by [`RLE`].
///
/// The type of `Item` is [`Run`].
/// Alphabet size `sigma == 0` indicates a large alphabet of unknown size.
/// Raw bytes and [`ByteCode`]-encoded integers can be read from the encoding using [`RLEIter::byte`] and [`RLEIter::int`].
///
/// # Examples
///
/// ```
/// use gbz::support::{Run, RLE, RLEIter};
///
/// let mut source = RLE::with_sigma(4);
/// source.write(Run::new(3, 12)); source.write(Run::new(2, 721)); source.write(Run::new(0, 34));
///
/// let mut iter = RLEIter::with_sigma(source.as_ref(), 4);
/// assert_eq!(iter.next(), Some(Run::new(3, 12)));
/// assert_eq!(iter.next(), Some(Run::new(2, 721)));
/// assert_eq!(iter.next(), Some(Run::new(0, 34)));
/// assert_eq!(iter.next(), None);
/// ```
#[derive(Clone, Debug)]
pub struct RLEIter<'a> {
    source: ByteCodeIter<'a>,
    sigma: usize,
    threshold: usize,
}

impl<'a> RLEIter<'a> {
    /// Creates a new iterator over the byte slice with alphabet size `0`.
    pub fn new(bytes: &'a [u8]) -> Self {
        let (sigma, threshold) = RLE::sanitize(0);
        RLEIter {
            source: ByteCodeIter::new(bytes),
            sigma,
            threshold,
        }
    }

    /// Creates a new iterator.
    ///
    /// # Arguments
    ///
    /// * `bytes`: Byte slice.
    /// * `sigma`: Alphabet size.
    pub fn with_sigma(bytes: &'a [u8], sigma: usize) -> Self {
        let (sigma, threshold) = RLE::sanitize(sigma);
        RLEIter {
            source: ByteCodeIter::new(bytes),
            sigma,
            threshold,
        }
    }

    /// Returns the next byte from the slice, or [`None`] if there are no more bytes left.
    #[inline]
    pub fn byte(&mut self) -> Option<u8> {
        self.source.byte()
    }

    /// Returns the next [`ByteCode`]-encoded integer from the slice, or [`None`] if no more integers can be decoded.
    #[inline]
    pub fn int(&mut self) -> Option<usize> {
        self.source.next()
    }

    /// Returns the first unvisited offset in the byte slice.
    #[inline]
    pub fn offset(&self) -> usize {
        self.source.offset()
    }

    /// Returns the alphabet size.
    #[inline]
    pub fn sigma(&self) -> usize {
        self.sigma
    }

    /// Changes the alphabet size to `sigma`.
    pub fn set_sigma(&mut self, sigma: usize) {
        let (sigma, threshold) = RLE::sanitize(sigma);
        self.sigma = sigma;
        self.threshold = threshold;
    }
}

impl<'a> Iterator for RLEIter<'a> {
    type Item = Run;

    fn next(&mut self) -> Option<Self::Item> {
        let mut run = Run::default();
        if self.sigma >= RLE::THRESHOLD {
            run.value = self.source.next()?;
            run.len = self.source.next()? + 1;
        } else {
            let byte = self.source.byte()?;
            run.value = (byte as usize) % self.sigma;
            run.len = (byte as usize) / self.sigma + 1;
            if run.len == self.threshold {
                run.len += self.source.next()?;
            }
        }
        Some(run)
    }
}

impl<'a> FusedIterator for RLEIter<'a> {}

//-----------------------------------------------------------------------------

/// A quick and dirty disjoint sets implementation.
///
/// Uses path splitting and union by rank.
/// The implementation is for all values in range `offset..offset + len`.
/// Each value is initially in a separate set.
/// The root element of the set containing a value can be retrieved with [`DisjointSets::find`].
/// The sets containing two elements can be joined with [`DisjointSets::union`].
/// When the sets are extracted with [`DisjointSets::extract`], some values may be omitted.
///
/// # Examples
///
/// ```
/// use gbz::support::DisjointSets;
///
/// let mut sets = DisjointSets::new(7, 2);
/// assert_eq!(sets.len(), 7);
/// assert_eq!(sets.offset(), 2);
///
/// sets.union(3, 4);
/// sets.union(3, 5);
/// sets.union(5, 7);
/// assert_eq!(sets.find(7), 3 - sets.offset());
///
/// let sets = sets.extract(|value| value != 6);
/// assert_eq!(sets.len(), 3);
/// assert_eq!(sets[0], vec![2]);
/// assert_eq!(sets[1], vec![3, 4, 5, 7]);
/// assert_eq!(sets[2], vec![8]);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DisjointSets {
    // Offset of the parent.
    parents: Vec<usize>,
    // Rank is at most ~log(size).
    ranks: Vec<u8>,
    // Value `i` is stored at offset `i - offset`.
    offset: usize,
}

impl DisjointSets {
    /// Returns a new `DisjointSets` structure with the given length and starting offset.
    ///
    /// # Panics
    ///
    /// Panics if `len + offset > usize::MAX`.
    pub fn new(len: usize, offset: usize) -> Self {
        if len > usize::MAX - offset {
            panic!("DisjointSets: length {} + offset {} too large", len, offset);
        }
        DisjointSets {
            parents: (0..len).collect(),
            ranks: vec![0; len],
            offset,
        }
    }

    /// Returns the number of values in the structure.
    pub fn len(&self) -> usize {
        self.parents.len()
    }

    /// Returns the starting offset for the values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the starting offset for the values.
    pub fn offset(&self) -> usize {
        self.offset
    }

    /// Returns the root element for the set containing the given value.
    ///
    /// Replaces the parent of each element with the grandparent to speed up further queries.
    /// This is also known as path splitting.
    ///
    /// # Panics
    ///
    /// May panic if `value < self.offset()` or `value >= self.len() + self.offset`.
    pub fn find(&mut self, value: usize) -> usize {
        let mut value = value - self.offset;
        while self.parents[value] != value {
            let next = self.parents[value];
            self.parents[value] = self.parents[next];
            value = next;
        }
        value
    }

    /// Joins the sets containing values `a` and `b`.
    ///
    /// Uses union by rank.
    /// The rank of a set containing a single value is 0.
    /// The root of the set with the lower rank becomes the child of the root of the set with the higher rank.
    /// If the ranks are equal, `find(a)` becomes the root and the rank of the new set increases by 1.
    ///
    /// # Panics
    ///
    /// May panic if `value < self.offset()` or `value >= self.len() + self.offset`, for values `a` and `b`.
    pub fn union(&mut self, a: usize, b: usize) {
        let mut a = self.find(a);
        let mut b = self.find(b);
        if a == b {
            return;
        }
        if self.ranks[a] < self.ranks[b] {
            mem::swap(&mut a, &mut b);
        }
        self.parents[b] = a;
        if self.ranks[b] == self.ranks[a] {
            self.ranks[a] += 1;
        }
    }

    /// Returns the sets corresponding to this structure.
    ///
    /// Only includes values for which `include_value(value) == true`.
    /// The sets will be sorted by the minimum value, and each set will be in sorted order.
    pub fn extract<F: Fn(usize) -> bool>(&mut self, include_value: F) -> Vec<Vec<usize>> {
        let mut result: Vec<Vec<usize>> = Vec::new();
        let mut root_to_set: HashMap<usize, usize> = HashMap::new();

        for value in self.offset..self.len() + self.offset {
            if !include_value(value) {
                continue;
            }
            let root = self.find(value);
            match root_to_set.entry(root) {
                Entry::Occupied(e) => {
                    result[*e.get()].push(value);
                },
                Entry::Vacant(e) => {
                    e.insert(result.len());
                    result.push(vec![value]);
                },
            }
        }

        result
    }
}

//-----------------------------------------------------------------------------

/// A sorted list of edges stored as [`SmallPos`].
///
/// If there are at most [`EdgeList::SMALL_CAPACITY`] edges, they are stored inline in a sorted array.
/// Larger edge sets are stored as a map from node to offset.
///
/// # Examples
///
/// ```
/// use gbz::support::EdgeList;
/// use gbz::Pos;
///
/// let mut edges = EdgeList::new();
/// edges.increment(3, 10);
/// edges.increment(1, 5);
/// edges.increment(2, 7);
/// edges.increment(1, 6);
/// edges.increment(4, 12); // Triggers conversion to large representation.
///
/// assert_eq!(edges.len(), 4);
/// assert_eq!(edges.get(1), Some(11));
/// assert_eq!(edges.get(5), None);
///
/// let truth = vec![
///     Pos::new(1, 11),
///     Pos::new(2, 7),
///     Pos::new(3, 10),
///     Pos::new(4, 12)
/// ];
/// let edges: Vec<Pos> = edges.iter().collect();
/// assert_eq!(edges, truth);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EdgeList {
    // Up to three edges as (len, edges), with the edges in sorted order.
    Small(u8, [SmallPos; Self::SMALL_CAPACITY]),
    // More than three edges as a map from node to offset.
    Large(BTreeMap<u32, u32>),
}

impl EdgeList {
    /// Capacity of the small representation.
    pub const SMALL_CAPACITY: usize = 3;

    /// Creates an empty edge list.
    pub fn new() -> Self {
        EdgeList::Small(0, [SmallPos::default(); Self::SMALL_CAPACITY])
    }

    /// Returns the number of edges in the list.
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            EdgeList::Small(len, _) => *len as usize,
            EdgeList::Large(map) => map.len(),
        }
    }

    /// Returns `true` if the list is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Increments the offset in the given edge by the given amount.
    ///
    /// Returns the offset before the increment, or `0` if the edge did not exist.
    /// Inserts a new edge if it does not already exist.
    pub fn increment(&mut self, node: usize, amount: usize) -> u32 {
        match self {
            EdgeList::Small(len, edges) => {
                for edge in edges.iter_mut().take(*len as usize) {
                    if edge.node == node as u32 {
                        let val = edge.offset;
                        edge.offset += amount as u32;
                        return val;
                    }
                }
                if *len < Self::SMALL_CAPACITY as u8 {
                    edges[*len as usize] = SmallPos::new(node, amount);
                    // Bubble up the new edge to maintain sorted order.
                    for i in (1..=*len as usize).rev() {
                        if edges[i] < edges[i - 1] {
                            edges.swap(i, i - 1);
                        } else {
                            break;
                        }
                    }
                    *len += 1;
                } else {
                    // Convert to large representation.
                    let mut map = BTreeMap::new();
                    for edge in edges.iter().take(*len as usize) {
                        map.insert(edge.node, edge.offset);
                    }
                    map.insert(node as u32, amount as u32);
                    *self = EdgeList::Large(map);
                }
                0
            }
            EdgeList::Large(map) => {
                let entry = map.entry(node as u32).or_insert(0);
                let val = *entry;
                *entry += amount as u32;
                val
            }
        }
    }

    /// Replaces the offset in each edge with the rank of the node in the set of nodes.
    pub fn set_ranks(&mut self) {
        match self {
            EdgeList::Small(len, edges) => {
                for (i, edge) in edges.iter_mut().take(*len as usize).enumerate() {
                    edge.offset = i as u32;
                }
            }
            EdgeList::Large(map) => {
                for (i, offset) in map.values_mut().enumerate() {
                    *offset = i as u32;
                }
            }
        }
    }

    /// Sets the offset in each edge to `0`.
    pub fn clear_offsets(&mut self) {
        match self {
            EdgeList::Small(len, edges) => {
                for edge in edges.iter_mut().take(*len as usize) {
                    edge.offset = 0;
                }
            }
            EdgeList::Large(map) => {
                for offset in map.values_mut() {
                    *offset = 0;
                }
            }
        }
    }

    /// Returns the offset in the given edge, or [`None`] if it does not exist.
    pub fn get(&self, node: usize) -> Option<u32> {
        match self {
            EdgeList::Small(len, edges) => {
                for edge in edges.iter().take(*len as usize) {
                    if edge.node == node as u32 {
                        return Some(edge.offset);
                    }
                }
                None
            }
            EdgeList::Large(map) => map.get(&(node as u32)).copied(),
        }
    }

    /// Returns a mutable reference to the offset in the given edge, or [`None`] if it does not exist.
    pub fn get_mut(&mut self, node: usize) -> Option<&mut u32> {
        match self {
            EdgeList::Small(len, edges) => {
                for edge in edges.iter_mut().take(*len as usize) {
                    if edge.node == node as u32 {
                        return Some(&mut edge.offset);
                    }
                }
                None
            }
            EdgeList::Large(map) => map.get_mut(&(node as u32)),
        }
    }

    /// Returns an iterator over the edges in the list, in sorted order.
    pub fn iter(&self) -> EdgeListIter<'_> {
        match self {
            EdgeList::Small(_, _) => EdgeListIter {
                parent: self,
                index: 0,
                iter: std::collections::btree_map::Iter::<'_, u32, u32>::default(),
            },
            EdgeList::Large(map) => EdgeListIter {
                parent: self,
                index: 0,
                iter: map.iter(),
            },
        }
    }
}

impl Default for EdgeList {
    fn default() -> Self {
        Self::new()
    }
}

/// An iterator over the edges in an [`EdgeList`].
///
/// The value of `Item` is [`Pos`].
pub struct EdgeListIter<'a> {
    parent: &'a EdgeList,
    index: usize,
    iter: std::collections::btree_map::Iter<'a, u32, u32>,
}

impl<'a> Iterator for EdgeListIter<'a> {
    type Item = Pos;

    fn next(&mut self) -> Option<Self::Item> {
        match self.parent {
            EdgeList::Small(len, edges) => {
                if self.index < *len as usize {
                    let pos = edges[self.index];
                    self.index += 1;
                    Some(Pos::from(pos))
                } else {
                    None
                }
            }
            EdgeList::Large(map) => {
                if self.index < map.len() {
                    self.index += 1;
                }
                self.iter.next().map(|(node, offset)| Pos::new(*node as usize, *offset as usize))
            }
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.parent.len() - self.index;
        (len, Some(len))
    }
}

impl ExactSizeIterator for EdgeListIter<'_> {}

impl FusedIterator for EdgeListIter<'_> {}

//-----------------------------------------------------------------------------

// TODO: Add trivial chains and components to the serialized format.
/// A set of top-level chains represented as links between boundary nodes.
///
/// Top-level chains provide a linear high-level structure for each weakly connected component in the graph.
/// A chain is a sequence of nodes and snarls.
/// Boundary nodes bordering the snarls form a sketch of graph topology.
/// Given a pair of boundary nodes, the graph region between them is either a unary path or a snarl.
/// In both cases, no path can leave the region without visiting one of the boundary nodes.
///
/// This representation is based on storing links between successive boundary nodes.
/// Each link is stored twice, once in each orientation.
///
/// # Serialization
///
/// `Chains` implements Simple-SDS serialization, but the format may still change.
/// It is currently intended for extracting top-level chains from a vg distance index or snarl decomposition and using them in GBZ-base.
///
/// The header contains the number of chains as an [`usize`] element.
/// Each chain is stored as an [`IntVector`] storing the sequence of oriented boundary nodes as GBWT node identifiers / handles.
/// Each chain is expected to be in the canonical orientation and to have the minimal necessary bit width for the items.
/// That generally means that most nodes are in the forward orientation and the sequence of node identifiers is mostly increasing.
/// The chains are expected to be sorted in lexicographic order.
///
/// The serialized format does not store the number of trivial chains and components.
///
/// # Examples
///
/// ```
/// use gbz::support::{self, Chains, Orientation};
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("micb-kir3dl1.chains");
/// let chains = serialize::load_from(&filename);
/// assert!(chains.is_ok());
/// let chains: Chains = chains.unwrap();
///
/// assert_eq!(chains.len(), 2);
/// assert_eq!(chains.links(), 925);
/// let handle = support::encode_node(44, Orientation::Forward);
/// assert!(chains.has_handle(handle));
/// let next = support::encode_node(47, Orientation::Forward);
/// assert_eq!(chains.next(handle), Some(next));
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Chains {
    chains: usize,
    trivial_chains: Option<usize>,
    components: Option<usize>,
    next: BTreeMap<usize, usize>,
}

impl Chains {
    // Reads the serialized chains representation.
    fn read_data<R: Read>(reader: &mut R) -> io::Result<Vec<IntVector>> {
        let chains = usize::load(reader)?;
        let mut data: Vec<IntVector> = Vec::with_capacity(chains);
        for _ in 0..chains {
            let vec = IntVector::load(reader)?;
            data.push(vec);
        }
        Ok(data)
    }

    // Converts the chains to a bidirectional link map.
    fn link_map(data: Vec<IntVector>) -> io::Result<BTreeMap<usize, usize>> {
        let mut next = BTreeMap::new();
        for chain in data {
            for i in 1..chain.len() {
                let from = chain.get(i - 1) as usize;
                if next.contains_key(&from) {
                    let msg = format!("Duplicate link from {}", from);
                    return Err(Error::new(ErrorKind::InvalidData, msg));
                }
                let to = chain.get(i) as usize;
                next.insert(from, to);

                let rev_from = flip_node(to);
                if next.contains_key(&rev_from) {
                    let msg = format!("Duplicate link from {}", rev_from);
                    return Err(Error::new(ErrorKind::InvalidData, msg));
                }
                let rev_to = flip_node(from);
                next.insert(rev_from, rev_to);
            }
        }
        Ok(next)
    }

    // Returns the set of head/tail nodes in the orientation pointing inward.
    fn head_tail_nodes(&self) -> BTreeSet<usize> {
        let mut result = BTreeSet::new();
        for &handle in self.next.keys() {
            if !self.next.contains_key(&flip_node(handle)) {
                result.insert(handle);
            }
        }
        result
    }

    // Converts the link map back to chains.
    fn links_to_chains(&self) -> Vec<IntVector> {
        let mut active = self.head_tail_nodes();

        // Iterate over all chains, starting from the smallest remaining head/tail.
        let mut result = Vec::new();
        while let Some(head) = active.first().copied() {
            let mut path = vec![head];
            let mut curr = head;
            let mut max = head;
            while let Some(next) = self.next(curr) {
                path.push(next);
                curr = next;
                if next > max {
                    max = next;
                }
            }
            if let Some(head) = path.first() {
                active.remove(head);
            }
            if let Some(tail) = path.last() {
                active.remove(&flip_node(*tail));
            }

            // Encode the path in the canonical orientation.
            if !encoded_path_is_canonical(&path) {
                reverse_path_in_place(&mut path);
            }
            let width = simple_sds::bits::bit_len(max as u64);
            let mut packed = IntVector::with_capacity(path.len(), width).unwrap();
            packed.extend(path);
            result.push(packed);
        }

        result
    }

    /// Creates an empty set of chains.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a new bidirectional link between the two handles.
    ///
    /// That means a link from `from` to `to` and from `flip_node(to)` to `flip_node(from)`.
    /// Returns `true` if the links were added and `false` if there was already a link from `from` or `flip_node(to)`.
    ///
    /// New links may change the number of chains.
    /// After all links have been added, the number of chains should be determined using [`Self::count_chains`].
    ///
    /// # Examples
    ///
    /// ```
    /// use gbz::support::Chains;
    ///
    /// let mut chains = Chains::new();
    /// let _ = chains.add_link(2, 4);
    /// let _ = chains.add_link(4, 6);
    /// let _ = chains.add_link(12, 14);
    /// let _ = chains.add_link(14, 17);
    ///
    /// // This will fail, because link (2, 4) is also link (5, 3).
    /// let result = chains.add_link(5, 3);
    /// assert!(!result);
    ///
    /// chains.count_chains();
    /// assert_eq!(chains.len(), 2);
    /// assert_eq!(chains.links(), 4);
    /// // This is the other orientation of (14, 17).
    /// assert_eq!(chains.next(16), Some(15));
    /// ```
    pub fn add_link(&mut self, from: usize, to: usize) -> bool {
        if self.next.contains_key(&from) || self.next.contains_key(&flip_node(to)) {
            return false;
        }
        self.next.insert(from, to);
        self.next.insert(flip_node(to), flip_node(from));
        true
    }

    /// Determines the number of chains implied by the links.
    ///
    /// This should be called after using [`Self::add_link`].
    pub fn count_chains(&mut self) {
        let heads_tails = self.head_tail_nodes();
        self.chains = heads_tails.len() / 2;
    }

    /// Sets the number of trivial chains that do not contain any links.
    pub fn set_trivial_chains(&mut self, trivial_chains: Option<usize>) {
        self.trivial_chains = trivial_chains;
    }

    /// Sets the number of weakly connected components in the graph.
    pub fn set_components(&mut self, components: Option<usize>) {
        self.components = components;
    }

    /// Returns the number of chains.
    pub fn len(&self) -> usize {
        self.chains
    }

    /// Returns the number of trivial chains that do not contain any links.
    pub fn trivial_chains(&self) -> Option<usize> {
        self.trivial_chains
    }

    /// Returns the number of weakly connected components in the graph, if known.
    ///
    /// In the ideal case, this should be the same as the number of chains, including trivial ones.
    /// But in more complex graphs, there can be multiple chains in a component.
    pub fn components(&self) -> Option<usize> {
        self.components
    }

    /// Returns `true` if there are no chains.
    pub fn is_empty(&self) -> bool {
        self.chains == 0
    }

    /// Returns the total number of links in the chains.
    pub fn links(&self) -> usize {
        self.next.len() / 2
    }

    /// Returns the successor for the given handle in the chains, or [`None`] if there is no successor.
    pub fn next(&self, handle: usize) -> Option<usize> {
        self.next.get(&handle).copied()
    }

    /// Returns `true` if the given node is a boundary node in one of the chains.
    pub fn has_node(&self, node_id: usize) -> bool {
        let fw_handle = encode_node(node_id, Orientation::Forward);
        let rev_handle = encode_node(node_id, Orientation::Reverse);
        self.next.contains_key(&fw_handle) || self.next.contains_key(&rev_handle)
    }

    /// Returns `true` if the given handle refers to a boundary node.
    pub fn has_handle(&self, handle: usize) -> bool {
        let rev_handle = flip_node(handle);
        self.next.contains_key(&handle) || self.next.contains_key(&rev_handle)
    }

    /// Returns an iterator over the links, ordered by source handle.
    ///
    /// Filter using [`encoded_edge_is_canonical`] to visit each link in a single orientation.
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.next.iter().map(|(k, v)| (*k, *v))
    }
}

impl Serialize for Chains {
    fn serialize_header<T: Write>(&self, writer: &mut T) -> io::Result<()> {
        self.chains.serialize(writer)
    }

    fn serialize_body<T: Write>(&self, writer: &mut T) -> io::Result<()> {
        let data = self.links_to_chains();
        if data.len() != self.chains {
            let msg = format!("Expected {} chains, but found {}", self.chains, data.len());
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }
        for chain in data.iter() {
            chain.serialize(writer)?;
        }
        Ok(())
    }

    fn load<T: Read>(reader: &mut T) -> io::Result<Self> {
        let data = Self::read_data(reader)?;
        let chains = data.len();
        let next = Self::link_map(data)?;
        Ok(Self { chains, trivial_chains: None, components: None, next })
    }

    fn size_in_elements(&self) -> usize {
        let mut result = 1;
        let data = self.links_to_chains();
        for chain in data.iter() {
            result += chain.size_in_elements();
        }
        result
    }
}

//-----------------------------------------------------------------------------
