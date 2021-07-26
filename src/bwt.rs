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
//! assert_eq!(node_2.outdegree(), 2);
//! assert_eq!(node_2.successor(1), 5);
//! assert_eq!(node_2.offset(1), 0);
//! assert_eq!(node_2.len(), 2);
//! assert_eq!(node_2.lf(1), Some((5, 0)));
//! assert_eq!(node_2.follow(0..2, 5), Some(0..1));
//! ```

use crate::support::{ByteCodeIter, RLE, RLEIter};
use crate::ENDMARKER;

use simple_sds::sparse_vector::{SparseVector, SparseBuilder};
use simple_sds::ops::{BitVec, Select};
use simple_sds::serialize::Serialize;

use std::cmp::Ordering;
use std::convert::TryFrom;
use std::io::{Error, ErrorKind};
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
        Record::new(&self.data[start..limit])
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

/// A partially decompressed node record.
///
/// See module-level documentation for an example.
#[derive(Clone, Debug)]
pub struct Record<'a> {
    edges: Vec<(usize, usize)>,
    bwt: &'a [u8],
}

impl<'a> Record<'a> {
    /// Returns a record corresponding to the byte slice, or `None` if the record would be empty.
    pub fn new(bytes: &'a [u8]) -> Option<Self> {
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
            edges: edges,
            bwt: &bytes[iter.offset()..],
        })
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

    /// Follows all sequences in the offset range to the given node.
    ///
    /// Returns a semiopen offset range in the destination node, or [`None`] if no such sequences exist.
    ///
    /// # Arguments
    ///
    /// * `range`: Offset range in the record.
    /// * `node`: Destination node.
    pub fn follow(&self, range: Range<usize>, node: usize) -> Option<Range<usize>> {
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
            if offset < range.end {
                if c == rank {
                    if offset < range.start {
                        result.start += if offset + len > range.start { range.start - offset } else { len };
                    }
                    result.end += if offset + len > range.end { range.end - offset } else { len };
                }
            }
            offset += len;
            if offset >= range.end {
                break;
            }
        }

        if result.is_empty() { None } else { Some(result) }
    }
}

//-----------------------------------------------------------------------------
