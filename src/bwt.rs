//! The BWT stored as an array of compressed node records.
//!
//! # Examples
//!
//! ```
//! use gbwt::Pos;
//! use gbwt::bwt::{BWT, BWTBuilder};
//! use gbwt::support::Run;
//!
//! // Encode the GBWT example from the paper.
//! let mut builder = BWTBuilder::new();
//! builder.append(&[Pos::new(1, 0)], &[Run::new(0, 3)]);
//! builder.append(&[Pos::new(2, 0), Pos::new(3, 0)], &[Run::new(0, 2), Run::new(1, 1)]);
//! builder.append(&[Pos::new(4, 0), Pos::new(5, 0)], &[Run::new(0, 1), Run::new(1, 1)]);
//! builder.append(&[Pos::new(4, 1)], &[Run::new(0, 1)]);
//! builder.append(&[Pos::new(5, 1), Pos::new(6, 0)], &[Run::new(1, 1), Run::new(0, 1)]);
//! builder.append(&[Pos::new(7, 0)], &[Run::new(0, 2)]);
//! builder.append(&[Pos::new(7, 2)], &[Run::new(0, 1)]);
//! builder.append(&[Pos::new(0, 0)], &[Run::new(0, 3)]);
//!
//! let bwt = BWT::from(builder);
//! assert_eq!(bwt.len(), 8);
//!
//! let node_2 = bwt.record(2).unwrap();
//! assert_eq!(node_2.id(), 2);
//! assert_eq!(node_2.outdegree(), 2);
//! assert_eq!(node_2.successor(1), 5);
//! assert_eq!(node_2.offset(1), 0);
//! assert_eq!(node_2.len(), 2);
//! assert_eq!(node_2.lf(1), Some(Pos::new(5, 0)));
//! assert_eq!(node_2.follow(0..2, 5), Some(0..1));
//!
//! // Determine the length of the BWT by iterating over the records.
//! let bwt_len = bwt.iter().fold(0, |len, record| len + record.len());
//! assert_eq!(bwt_len, 17);
//!
//! // Collect the record identifiers.
//! let ids: Vec<usize> = bwt.id_iter().collect();
//! assert_eq!(ids, vec![0, 1, 2, 3, 4, 5, 6, 7]);
//! ```

use crate::support::{ByteCodeIter, Run, RLE, RLEIter};
use crate::ENDMARKER;
use crate::support;

use simple_sds::sparse_vector::{SparseVector, SparseBuilder, OneIter};
use simple_sds::ops::{BitVec, Select};
use simple_sds::serialize::Serialize;

use std::cmp::Ordering;
use std::convert::TryFrom;
use std::io::{Error, ErrorKind};
use std::iter::FusedIterator;
use std::ops::Range;
use std::io;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// A GBWT position as a (node, offset) pair.
#[derive(Copy, Clone, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Pos {
    /// GBWT node identifier.
    pub node: usize,
    /// BWT offset within the node.
    pub offset: usize,
}

impl Pos {
    /// Creates a new position.
    #[inline]
    pub fn new(node: usize, offset: usize) -> Self {
        Pos {
            node, offset,
        }
    }
}

impl From<(usize, usize)> for Pos {
    #[inline]
    fn from(pos: (usize, usize)) -> Self {
        Self::new(pos.0, pos.1)
    }
}

//-----------------------------------------------------------------------------

/// The BWT encoded as a vector of bytes.
///
/// The encoding consists of `self.len()` concatenated node records.
/// Record identifiers are characters in the effective alphabet `0..self.len()`, but they are not necessarily the same as the node identifiers.
/// There may be empty records that do not correspond to any node in the graph.
/// See module-level documentation for an example.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BWT {
    index: SparseVector,
    data: Vec<u8>,
}

impl BWT {
    /// Returns the number of records in the BWT.
    #[inline]
    pub fn len(&self) -> usize {
        self.index.count_ones()
    }

    /// Returns `true` if the BWT is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    // Returns the byte slice corresponding to record `i`, assuming that it exists.
    fn record_bytes(&self, i: usize) -> &[u8] {
        let mut iter = self.index.select_iter(i);
        let (_, start) = iter.next().unwrap();
        let limit = if i + 1 < self.len() { iter.next().unwrap().1 } else { self.data.len() };
        &self.data[start..limit]
    }

    /// Returns the `i`th record, or [`None`] if there is no such node.
    pub fn record(&'_ self, i: usize) -> Option<Record<'_>> {
        if i >= self.len() {
            return None;
        }
        let bytes = self.record_bytes(i);
        Record::new(i, bytes)
    }

    /// Returns the compressed representation of the `i`th record, partitioned
    /// into edges and BWT, or [`None`] if there is no such record.
    pub fn compressed_record(&self, i: usize) -> Option<(&[u8], &[u8])> {
        if i >= self.len() {
            return None;
        }

        let bytes = self.record_bytes(i);
        let offset = Record::skip_edges(bytes)?;

        Some((&bytes[..offset], &bytes[offset..]))
    }

    /// Returns an iterator over the records in the BWT.
    ///
    /// Note that the iterator skips empty records.
    pub fn iter(&'_ self) -> RecordIter<'_> {
        RecordIter {
            parent: self,
            next: 0,
        }
    }

    /// Returns an iterator over the identifiers of non-empty records in the BWT.
    pub fn id_iter(&'_ self) -> IdIter<'_> {
        IdIter {
            parent: self,
            iter: self.index.one_iter(),
            next: 0,
        }
    }
}

impl Serialize for BWT {
    fn serialize_header<T: io::Write>(&self, _: &mut T) -> io::Result<()> {
        Ok(())
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.index.serialize(writer)?;
        self.data.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let index = SparseVector::load(reader)?;
        let data = Vec::<u8>::load(reader)?;
        if index.len() != data.len() {
            return Err(Error::new(ErrorKind::InvalidData, "BWT: Index / data length mismatch"));
        }
        Ok(BWT {
            index, data,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.index.size_in_elements() + self.data.size_in_elements()
    }
}

impl From<BWTBuilder> for BWT {
    fn from(source: BWTBuilder) -> Self {
        let mut builder = SparseBuilder::new(source.encoder.len(), source.offsets.len()).unwrap();
        for offset in source.offsets.iter() {
            unsafe { builder.set_unchecked(*offset); }
        }
        BWT {
            index: SparseVector::try_from(builder).unwrap(),
            data: Vec::<u8>::from(source.encoder),
        }
    }
}

//-----------------------------------------------------------------------------

/// A structure for building the BWT by appending node records.
///
/// This is mostly inteded for testing at the moment, as no BWT construction algorithms have been implemented.
/// See module-level documentation for an example.
#[derive(Clone, Debug, Default)]
pub struct BWTBuilder {
    offsets: Vec<usize>,
    encoder: RLE,
}

impl BWTBuilder {
    /// Creates a new builder.
    pub fn new() -> Self {
        BWTBuilder::default()
    }

    /// Returns the number of records.
    #[inline]
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Returns true if the builder is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Appends a new record to the BWT.
    ///
    /// The record consists of a list of edges and a list of runs.
    /// Each edge is a position in the successor node.
    /// Note that the successor node is given by its node id, which is not necessarily the same as its record id.
    /// Each run must have `run.value < edges.len()` and `run.len > 0`.
    pub fn append(&mut self, edges: &[Pos], runs: &[Run]) {
        self.offsets.push(self.encoder.len());
        self.encoder.write_int(edges.len());
        let mut prev = 0;
        for edge in edges {
            self.encoder.write_int(edge.node - prev); self.encoder.write_int(edge.offset);
            prev = edge.node;
        }
        self.encoder.set_sigma(edges.len());
        for run in runs {
            self.encoder.write(*run);
        }
    }
}

//-----------------------------------------------------------------------------

/// An iterator over the records in [`BWT`].
///
/// The type of `Item` is [`Record`].
/// Note that the iterator skips empty records.
/// See module-level documentation for an example.
#[derive(Clone, Debug)]
pub struct RecordIter<'a> {
    parent: &'a BWT,
    // The first index we have not visited.
    next: usize,
}

impl<'a> Iterator for RecordIter<'a> {
    type Item = Record<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.next < self.parent.len() {
            let result = self.parent.record(self.next);
            self.next += 1;
            if result.is_some() {
                return result;
            }
        }
        None
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.parent.len() - self.next))
    }
}

impl<'a> FusedIterator for RecordIter<'a> {}

//-----------------------------------------------------------------------------

/// An iterator over the identifiers of non-empty records in [`BWT`].
///
/// The type of `Item` is [`usize`].
/// See module-level documentation for an example.
#[derive(Clone, Debug)]
pub struct IdIter<'a> {
    parent: &'a BWT,
    iter: OneIter<'a>,
    // The first index we have not visited.
    next: usize,
}

impl<'a> Iterator for IdIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        for (node, offset) in self.iter.by_ref() {
            self.next = node + 1;
            if self.parent.data[offset] != 0 {
                return Some(node);
            }
        }
        None
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.parent.len() - self.next))
    }
}

impl<'a> FusedIterator for IdIter<'a> {}

//-----------------------------------------------------------------------------

/// A partially decompressed node record.
///
/// See module-level documentation for an example.
#[derive(Clone, Debug)]
pub struct Record<'a> {
    id: usize,
    edges: Vec<Pos>,
    bwt: &'a [u8],
}

impl<'a> Record<'a> {
    /// Returns a record corresponding to the byte slice, or [`None`] if the record would be empty.
    pub fn new(id: usize, bytes: &'a [u8]) -> Option<Self> {
        if bytes.is_empty() {
            return None;
        }
        let (edges, offset) = Self::decompress_edges(bytes)?;
        Some(Record {
            id,
            edges,
            bwt: &bytes[offset..],
        })
    }

    /// Returns a record composed of the given parts.
    ///
    /// # Safety
    ///
    /// The record will be invalid if `edges` is empty, the edges are not in sorted order,
    /// `bwt` does not encode a sequence of runs using [`RLE`] with alphabet size `edges.len()`,
    /// or the offsets implied by LF-mapping do not exist in the successor records.
    pub unsafe fn from_raw_parts(id: usize, edges: Vec<Pos>, bwt: &'a [u8]) -> Self {
        Record { id, edges, bwt }
    }

    /// Decomposes the record into raw parts.
    pub fn into_raw_parts(self) -> (usize, Vec<Pos>, &'a [u8]) {
        (self.id, self.edges, self.bwt)
    }

    /// Decompresses the adjacency list from a byte slice.
    ///
    /// Returns the list of edges and the slice offset after the adjacency
    /// list, or [`None`] if the list is empty.
    ///
    /// # Panics
    ///
    /// May panic if slice does not encode an adjacency list or if the slice
    /// ends early.
    pub fn decompress_edges(bytes: &[u8]) -> Option<(Vec<Pos>, usize)> {
        let mut iter = ByteCodeIter::new(bytes);
        let sigma = iter.next().unwrap();
        if sigma == 0 {
            return None;
        }

        let mut edges: Vec<Pos> = Vec::new();
        let mut prev = 0;
        for _ in 0..sigma {
            let node = iter.next().unwrap() + prev;
            prev = node;
            let offset = iter.next().unwrap();
            edges.push(Pos::new(node, offset));
        }

        Some((edges, iter.offset()))
    }

    // Skips the list of edges in the byte slice and returns the offset
    // after the edge list.
    fn skip_edges(bytes: &[u8]) -> Option<usize> {
        let mut iter = ByteCodeIter::new(bytes);
        let sigma = iter.next().unwrap();
        if sigma == 0 {
            return None;
        }

        for _ in 0..sigma {
            let _ = iter.next().unwrap();
            let _ = iter.next().unwrap();
        }

        Some(iter.offset())
    }

    /// Returns the identifier of the record.
    pub fn id(&self) -> usize {
        self.id
    }

    /// Returns the outdegree of the node.
    #[inline]
    pub fn outdegree(&self) -> usize {
        self.edges.len()
    }

    /// Returns the successor node of rank `i`.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.outdegree()`.
    #[inline]
    pub fn successor(&self, i: usize) -> usize {
        self.edges[i].node
    }

    /// Returns the BWT offset in the node of rank `i`.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.outdegree()`.
    #[inline]
    pub fn offset(&self, i: usize) -> usize {
        self.edges[i].offset
    }

    /// Returns the length of the offset range.
    ///
    /// This is somewhat slow, as it requires iterating over the run-lenght encoded BWT slice.
    /// Note that the length is always non-zero.
    pub fn len(&self) -> usize {
        let mut result = 0;
        for run in RLEIter::with_sigma(self.bwt, self.edges.len()) {
            result += run.len;
        }
        result
    }

    /// Returns `false`.
    ///
    /// Keeps Clippy happy.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Decompress the record as a vector of successor positions.
    pub fn decompress(&self) -> Vec<Pos> {
        let mut edges = self.edges.clone();
        let mut result: Vec<Pos> = Vec::new();
        for run in RLEIter::with_sigma(self.bwt, self.edges.len()) {
            for _ in 0..run.len {
                result.push(edges[run.value]);
                edges[run.value].offset += 1;
            }
        }
        result
    }

    /// Follows the sequence at offset `i` and returns the successor position.
    ///
    /// Returns [`None`] if the sequence ends or offset `i` does not exist.
    pub fn lf(&self, i: usize) -> Option<Pos> {
        let mut edges = self.edges.clone();
        let mut offset = 0;
        for run in RLEIter::with_sigma(self.bwt, self.edges.len()) {
            if offset + run.len > i {
                if self.successor(run.value) == ENDMARKER {
                    return None;
                } else {
                    edges[run.value].offset += i - offset;
                    return Some(edges[run.value]);
                }
            }
            edges[run.value].offset += run.len;
            offset += run.len;
        }
        None
    }

    /// Returns the predecessor node for the sequence at offset `i` in the other orientation of this node.
    ///
    /// This query assumes that the GBWT index is bidirectional.
    /// Returns [`None`] if the predecessor or the offset does not exist.
    pub fn predecessor_at(&self, i: usize) -> Option<usize> {
        // Determine the number of sequences going to each successor node.
        let mut edges: Vec<Pos> = Vec::with_capacity(self.edges.len());
        for rank in 0..self.edges.len() {
            edges.push(Pos::new(self.successor(rank), 0));
        }
        for run in RLEIter::with_sigma(self.bwt, self.edges.len()) {
            edges[run.value].offset += run.len;
        }

        // Flip the successor nodes to make them the predecessors of the other orientation of this node.
        for edge in &mut edges {
            if edge.node != ENDMARKER {
                edge.node = support::flip_node(edge.node);
            }
        }

        // Handle the special case where the predecessors are now in wrong order because they contain
        // both orientations of the same node.
        for rank in 1..edges.len() {
            if support::node_id(edges[rank - 1].node) == support::node_id(edges[rank].node) {
                edges.swap(rank - 1, rank);
            }
        }

        // Find the predecessor, if it exists.
        let mut offset = 0;
        for edge in edges {
            offset += edge.offset;
            if offset > i {
                if edge.node == ENDMARKER {
                    return None;
                }
                return Some(edge.node);
            }
        }

        None
    }

    // Returns the rank of the edge to the given node.
    fn edge_to(&self, node: usize) -> Option<usize> {
        let mut low = 0;
        let mut high = self.outdegree();
        while low < high {
            let mid = low + (high - low) / 2;
            match node.cmp(&self.edges[mid].node) {
                Ordering::Less => high = mid,
                Ordering::Equal => return Some(mid),
                Ordering::Greater => low = mid + 1,
            }
        }
        None
    }

    /// Returns the offset for which [`Self::lf`] would return `pos`, or [`None`] if no such offset exists.
    pub fn offset_to(&self, pos: Pos) -> Option<usize> {
        if pos.node == ENDMARKER {
            return None;
        }
        let outrank = self.edge_to(pos.node)?;

        // Rank of `pos.0` so far.
        let mut succ_rank = self.offset(outrank);
        if succ_rank > pos.offset {
            return None;
        }

        // Find the occurrence of `pos.0` of rank `pos.1 - succ_rank`.
        let mut offset = 0;
        for run in RLEIter::with_sigma(self.bwt, self.outdegree()) {
            offset += run.len;
            if run.value != outrank {
                continue;
            }
            succ_rank += run.len;
            if succ_rank > pos.offset {
                return Some(offset - (succ_rank - pos.offset));
            }
        }

        None
    }

    /// Follows all sequences in the offset range to the given node.
    ///
    /// Returns a semiopen offset range in the destination node, or [`None`] if no such sequences exist.
    /// See also [`Record::bd_follow`].
    ///
    /// # Arguments
    ///
    /// * `range`: Offset range in the record.
    /// * `node`: Destination node.
    pub fn follow(&self, range: Range<usize>, node: usize) -> Option<Range<usize>> {
        if range.is_empty() || node == ENDMARKER {
            return None;
        }
        let rank = self.edge_to(node)?;

        let mut result = self.offset(rank)..self.offset(rank);
        let mut offset = 0;
        for run in RLEIter::with_sigma(self.bwt, self.outdegree()) {
            if run.value == rank {
                let run_range = offset..offset + run.len;
                result.start += support::intersect(&run_range, &(0..range.start)).len();
                result.end += support::intersect(&run_range, &(0..range.end)).len();
            }
            offset += run.len;
            if offset >= range.end {
                break;
            }
        }

        if result.is_empty() { None } else { Some(result) }
    }

    /// Follows all sequences in the offset range to the given node.
    ///
    /// This query assumes that the GBWT index is bidirectional.
    /// Returns a semiopen offset range in the destination node, or [`None`] if no such sequences exist.
    /// The second return value is the number of occurrences of nodes `v` in the query range such that [`support::flip_node`]`(v) < `[`support::flip_node`]`(node)`.
    /// This information can be used for updating the reverse range in bidirectional search.
    /// See also [`Record::follow`].
    ///
    /// # Arguments
    ///
    /// * `range`: Offset range in the record.
    /// * `node`: Destination node.
    pub fn bd_follow(&self, range: Range<usize>, node: usize) -> Option<(Range<usize>, usize)> {
        if range.is_empty() || node == ENDMARKER {
            return None;
        }
        let rank = self.edge_to(node)?;
        let reverse = support::flip_node(node);

        let mut result = self.offset(rank)..self.offset(rank);
        let mut count = 0;
        let mut offset = 0;
        for run in RLEIter::with_sigma(self.bwt, self.outdegree()) {
            let run_range = offset..offset + run.len;
            if run.value == rank {
                result.start += support::intersect(&run_range, &(0..range.start)).len();
                result.end += support::intersect(&run_range, &(0..range.end)).len();
            }
            if support::flip_node(self.successor(run.value)) < reverse {
                count += support::intersect(&run_range, &range).len();
            }
            offset += run.len;
            if offset >= range.end {
                break;
            }
        }

        if result.is_empty() { None } else { Some((result, count)) }
    }
}

//-----------------------------------------------------------------------------
