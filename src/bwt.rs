//! The BWT stored as an array of compressed node records.
//!
//! # Examples
//!
//! ```
//! use gbwt::bwt::{BWT, BWTBuilder};
//!
//! // Encode the GBWT example from the paper.
//! let mut builder = BWTBuilder::new();
//! builder.append(&[(1, 0)], &[(0, 3)]);
//! builder.append(&[(2, 0), (3, 0)], &[(0, 2), (1, 1)]);
//! builder.append(&[(4, 0), (5, 0)], &[(0, 1), (1, 1)]);
//! builder.append(&[(4, 1)], &[(0, 1)]);
//! builder.append(&[(5, 1), (6, 0)], &[(1, 1), (0, 1)]);
//! builder.append(&[(7, 0)], &[(0, 2)]);
//! builder.append(&[(7, 2)], &[(0, 1)]);
//! builder.append(&[(0, 0)], &[(0, 3)]);
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
//! assert_eq!(node_2.lf(1), Some((5, 0)));
//! assert_eq!(node_2.follow(&(0..2), 5), Some(0..1));
//!
//! // Determine the length of the BWT by iterating over the records.
//! let bwt_len = bwt.iter().fold(0, |len, record| len + record.len());
//! assert_eq!(bwt_len, 17);
//! ```

use crate::support::{ByteCodeIter, RLE, RLEIter};
use crate::ENDMARKER;
use crate::support;

use simple_sds::sparse_vector::{SparseVector, SparseBuilder};
use simple_sds::ops::{BitVec, Select};
use simple_sds::serialize::Serialize;

use std::cmp::Ordering;
use std::convert::TryFrom;
use std::io::{Error, ErrorKind};
use std::iter::{FusedIterator};
use std::ops::Range;
use std::io;

#[cfg(test)]
mod tests;

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

    /// Returns the `i`th record, or `None` if there is no such node.
    pub fn record(&self, i: usize) -> Option<Record> {
        if i >= self.len() {
            return None;
        }
        let mut iter = self.index.select_iter(i);
        let (_, start) = iter.next().unwrap();
        let limit = if i + 1 < self.len() { iter.next().unwrap().1 } else { self.data.len() };
        Record::new(i, &self.data[start..limit])
    }

    /// Returns an iterator over the records in the BWT.
    ///
    /// Note that the iterator skips empty records.
    pub fn iter(&self) -> RecordIter {
        RecordIter {
            parent: self,
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
            index: index,
            data: data,
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
    /// Each edge is a pair (node id, BWT offset), where node id is not necessarily the same as record id.
    /// Each run is a pair `(rank, len)`, with `rank < edges.len()` and `len > 0`.
    pub fn append(&mut self, edges: &[(usize, usize)], runs: &[(usize, usize)]) {
        self.offsets.push(self.encoder.len());
        self.encoder.write_int(edges.len());
        let mut prev = 0;
        for (node, offset) in edges {
            self.encoder.write_int(*node - prev); self.encoder.write_int(*offset);
            prev = *node;
        }
        self.encoder.set_sigma(edges.len());
        for (rank, len) in runs {
            self.encoder.write(*rank, *len);
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

/// A partially decompressed node record.
///
/// See module-level documentation for an example.
#[derive(Clone, Debug)]
pub struct Record<'a> {
    id: usize,
    edges: Vec<(usize, usize)>,
    bwt: &'a [u8],
}

impl<'a> Record<'a> {
    /// Returns a record corresponding to the byte slice, or `None` if the record would be empty.
    pub fn new(id: usize, bytes: &'a [u8]) -> Option<Self> {
        if bytes.is_empty() {
            return None;
        }

        // Determine the outdegree.
        let mut iter = ByteCodeIter::new(bytes);
        let sigma = iter.next().unwrap();
        if sigma == 0 {
            return None;
        }

        // Decompress the edges.
        let mut edges: Vec<(usize, usize)> = Vec::new();
        let mut prev = 0;
        for _ in 0..sigma {
            let node = iter.next().unwrap() + prev;
            prev = node;
            let offset = iter.next().unwrap();
            edges.push((node, offset));
        }

        Some(Record {
            id: id,
            edges: edges,
            bwt: &bytes[iter.offset()..],
        })
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
        self.edges[i].0
    }

    /// Returns the BWT offset in the node of rank `i`.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.outdegree()`.
    #[inline]
    pub fn offset(&self, i: usize) -> usize {
        self.edges[i].1
    }

    /// Returns the length of the offset range.
    ///
    /// This is somewhat slow, as it requires iterating over the run-lenght encoded BWT slice.
    /// Note that the length is always non-zero.
    pub fn len(&self) -> usize {
        let mut result = 0;
        for (_, len) in RLEIter::with_sigma(self.bwt, self.edges.len()) {
            result += len;
        }
        result
    }

    /// Decompress the record as a vector of (successor node, offset in successor) pairs.
    pub fn decompress(&self) -> Vec<(usize, usize)> {
        let mut edges = self.edges.clone();
        let mut result: Vec<(usize, usize)> = Vec::new();
        for (rank, len) in RLEIter::with_sigma(self.bwt, self.edges.len()) {
            for _ in 0..len {
                result.push(edges[rank]);
                edges[rank].1 += 1;
            }
        }
        result
    }

    /// Follows the sequence at offset `i` and returns (successor node, offset in successor).
    ///
    /// Returns [`None`] if the sequence ends or offset `i` does not exist.
    pub fn lf(&self, i: usize) -> Option<(usize, usize)> {
        let mut edges = self.edges.clone();
        let mut offset = 0;
        for (rank, len) in RLEIter::with_sigma(self.bwt, self.edges.len()) {
            if offset + len > i {
                if self.successor(rank) == ENDMARKER {
                    return None;
                } else {
                    edges[rank].1 += i - offset;
                    return Some(edges[rank]);
                }
            }
            edges[rank].1 += len;
            offset += len;
        }
        None
    }

    /// Returns the predecessor node for the sequence at offset `i` in the other orientation of this node.
    ///
    /// This query assumes that the GBWT index is bidirectional.
    /// Returns [`None`] if the predecessor or the offset does not exist.
    pub fn predecessor_at(&self, i: usize) -> Option<usize> {
        // Determine the number of sequences going to each successor node.
        let mut edges: Vec<(usize, usize)> = Vec::with_capacity(self.edges.len());
        for rank in 0..self.edges.len() {
            edges.push((self.successor(rank), 0));
        }
        for (rank, len) in RLEIter::with_sigma(self.bwt, self.edges.len()) {
            edges[rank].1 += len;
        }

        // Flip the successor nodes to make them the predecessors of the other orientation of this node.
        for rank in 0..edges.len() {
            if edges[rank].0 != ENDMARKER {
                edges[rank].0 = support::flip_node(edges[rank].0);
            }
        }

        // Handle the special case where the predecessors are now in wrong order because they contain
        // both orientations of the same node.
        for rank in 1..edges.len() {
            if support::node_id(edges[rank - 1].0) == support::node_id(edges[rank].0) {
                edges.swap(rank - 1, rank);
            }
        }

        // Find the predecessor, if it exists.
        let mut offset = 0;
        for (id, count) in edges {
            offset += count;
            if offset > i {
                if id == ENDMARKER {
                    return None;
                }
                return Some(id);
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
            match node.cmp(&self.edges[mid].0) {
                Ordering::Less => high = mid,
                Ordering::Equal => return Some(mid),
                Ordering::Greater => low = mid + 1,
            }
        }
        None
    }

    /// Returns the offset for which [`Self::lf`] would return `pos`, or [`None`] if no such offset exists.
    pub fn offset_to(&self, pos: (usize, usize)) -> Option<usize> {
        if pos.0 == ENDMARKER {
            return None;
        }
        let outrank = self.edge_to(pos.0);
        if outrank == None {
            return None;
        }
        let outrank = outrank.unwrap();

        // Rank of `pos.0` so far.
        let mut succ_rank = self.offset(outrank);
        if succ_rank > pos.1 {
            return None;
        }

        // Find the occurrence of `pos.0` of rank `pos.1 - succ_rank`.
        let mut offset = 0;
        for (c, len) in RLEIter::with_sigma(&self.bwt, self.outdegree()) {
            offset += len;
            if c != outrank {
                continue;
            }
            succ_rank += len;
            if succ_rank > pos.1 {
                return Some(offset - (succ_rank - pos.1));
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
    pub fn follow(&self, range: &Range<usize>, node: usize) -> Option<Range<usize>> {
        if range.is_empty() || node == ENDMARKER {
            return None;
        }
        let rank = self.edge_to(node);
        if rank == None {
            return None;
        }
        let rank = rank.unwrap();

        let mut result = self.offset(rank)..self.offset(rank);
        let mut offset = 0;
        for (c, len) in RLEIter::with_sigma(&self.bwt, self.outdegree()) {
            if c == rank {
                let run = offset..offset + len;
                result.start += support::intersect(&run, &(0..range.start)).len();
                result.end += support::intersect(&run, &(0..range.end)).len();
            }
            offset += len;
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
    pub fn bd_follow(&self, range: &Range<usize>, node: usize) -> Option<(Range<usize>, usize)> {
        if range.is_empty() || node == ENDMARKER {
            return None;
        }
        let rank = self.edge_to(node);
        if rank == None {
            return None;
        }
        let rank = rank.unwrap();
        let reverse = support::flip_node(node);

        let mut result = self.offset(rank)..self.offset(rank);
        let mut count = 0;
        let mut offset = 0;
        for (c, len) in RLEIter::with_sigma(&self.bwt, self.outdegree()) {
            let run = offset..offset + len;
            if c == rank {
                result.start += support::intersect(&run, &(0..range.start)).len();
                result.end += support::intersect(&run, &(0..range.end)).len();
            }
            if support::flip_node(self.successor(c)) < reverse {
                count += support::intersect(&run, &range).len();
            }
            offset += len;
            if offset >= range.end {
                break;
            }
        }

        if result.is_empty() { None } else { Some((result, count)) }
    }
}

//-----------------------------------------------------------------------------
