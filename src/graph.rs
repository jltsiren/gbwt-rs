//! GBWTGraph: Node sequences and node-to-segment translation.
//!
//! The [`Graph`] structure augments a GBWT index with node sequences and an optional node-to-segment translation.
//! This enables representing bidirected sequence graphs compatible with a subset of the [GFA format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).
//! Unlike in the [C++ implementation](https://github.com/jltsiren/gbwtgraph), the actual graph interface is provided by the GBZ structure.
//!
//! At the moment, this implementation only supports graphs built with other tools.
// FIXME document, link to GBWT, GBZ

use crate::headers::{Header, GraphPayload};
use crate::support::{StringArray, StringIter};

use simple_sds::ops::{BitVec, Select};
use simple_sds::serialize::Serialize;
use simple_sds::sparse_vector::{SparseVector, OneIter};

use std::io::{Error, ErrorKind};
use std::iter::FusedIterator;
use std::ops::Range;
use std::io;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

// FIXME need new test case for translation
/// Node sequences and an optional node-to-segment translation for a graph based on a GBWT index.
///
/// The graph stores a number of sequences (node labels) with consecutive identifiers starting from `0`.
/// Sequence identifier `id` corresponds to node identifier `id + min_node` in the original graph, where `min_node` is the smallest node identifier.
/// If some node identifiers are not used in the graph, the corresponding sequences are usually empty.
///
/// There may also be a translation between node identifiers and GFA segments.
/// Each segment has a string name, and it corresponds to the concatenation of nodes with consecutive identifiers.
/// If the translation is in use, node identifiers start from `1`.
///
/// # Examples
///
/// ```
/// use gbwt::Graph;
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gg");
/// let graph: Graph = serialize::load_from(&filename).unwrap();
///
/// assert_eq!(graph.nodes(), 12);
/// assert_eq!(graph.sequences(), 15);
/// assert_eq!(graph.sequence(13), "T".as_bytes());
/// assert_eq!(graph.sequence_len(3), 1);
///
/// let labels: Vec<&[u8]> = graph.iter().collect();
/// assert_eq!(labels.concat(), "GATTACAGATTA".as_bytes());
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Graph {
    header: Header<GraphPayload>,
    sequences: StringArray,
    segments: StringArray,
    mapping: SparseVector,
}

//-----------------------------------------------------------------------------

/// Sequences for nodes.
impl Graph {
    /// Returns the number of nodes in the original graph.
    ///
    /// This may be less than the number of sequences, if there are gaps in the node id space.
    #[inline]
    pub fn nodes(&self) -> usize {
        self.header.payload().nodes
    }

    /// Returns the number of sequences in the graph.
    ///
    /// Sequence ids are in the range `0..self.sequences()`.
    /// If there are gaps in the node id space of the graph, sequences corresponding to unused ids may be empty.
    #[inline]
    pub fn sequences(&self) -> usize {
        self.sequences.len()
    }

    /// Returns `true` if the graph is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Returns the sequence with the given identifier.
    ///
    /// # Panics
    ///
    /// May panic if `id >= self.sequences()`.
    #[inline]
    pub fn sequence(&self, id: usize) -> &[u8] {
        self.sequences.bytes(id)
    }

    /// Returns the length of the sequence with the given identifier.
    ///
    /// # Panics
    ///
    /// May panic if `id >= self.sequences()`.
    #[inline]
    pub fn sequence_len(&self, id: usize) -> usize {
        self.sequences.str_len(id)
    }

    /// Returns an iterator over the sequences in the graph.
    ///
    /// If there are gaps in the node id space of the graph, sequences corresponding to unused id may be empty.
    pub fn iter(&self) -> StringIter {
        self.sequences.iter()
    }
}

//-----------------------------------------------------------------------------

// FIXME tests
/// Segments.
impl Graph {
    /// Returns `true` if the graph contains a node-to-segment translation.
    #[inline]
    pub fn has_translation(&self) -> bool {
        self.header.is_set(GraphPayload::FLAG_TRANSLATION)
    }

    /// Returns the number of segments in the original graph.
    ///
    /// Note that some segments may be missing from the GBZ graph if the corresponding nodes are not used on any path.
    /// Returns `0` if the graph does not contain a node-to-segment translation.
    #[inline]
    pub fn segments(&self) -> usize {
        self.segments.len()
    }

    /// Returns the name of the segment with the given identifier.
    ///
    /// The name is assumed to be valid UTF-8.
    /// It is returned as a byte slice to avoid unnecessary checks.
    ///
    /// # Panics
    ///
    /// May panic if `id >= self.segments()`.
    #[inline]
    pub fn segment_name(&self, id: usize) -> &[u8] {
        self.segments.bytes(id)
    }

    /// Returns the range of node ids in the original graph corresponding to the given segment identifier.
    ///
    /// See also [`Graph::segment_range`].
    /// Note that random access to the node-to-segment translation is somewhat slow.
    ///
    /// # Panics
    ///
    /// May panic if `id >= self.segments()`.
    pub fn segment_nodes(&self, id: usize) -> Range<usize> {
        let mut iter = self.mapping.select_iter(id);
        let (_, start) = iter.next().unwrap();
        let end = if id + 1 < self.segments() { iter.next().unwrap().1 } else { self.sequences() };
        start..end
    }

    /// Returns the range of sequences corresponding to the given segment identifier.
    ///
    /// See also [`Graph::segment_nodes`].
    /// Note that random access to the node-to-segment translation is somewhat slow.
    ///
    /// # Panics
    ///
    /// May panic if `id >= self.segments()`.
    pub fn segment_range(&self, id: usize) -> Range<usize> {
        let range = self.segment_nodes(id);
        range.start - 1..range.end - 1
    }

    /// Returns the concatenation of sequences corresponding to the given segment identifier.
    ///
    /// Note that random access to the node-to-segment translation is somewhat slow.
    ///
    /// # Panics
    ///
    /// May panic if `id >= self.segments()`.
    pub fn segment_sequence(&self, id: usize) -> &[u8] {
        let range = self.segment_range(id);
        self.sequences.range(range)
    }

    /// Returns the total length of the sequences corresponding to the given segment identifier.
    ///
    /// Note that random access to the node-to-segment translation is somewhat slow.
    ///
    /// # Panics
    ///
    /// May panic if `id >= self.segments()`.
    pub fn segment_len(&self, id: usize) -> usize {
        let range = self.segment_range(id);
        self.sequences.range_len(range)
    }

    /// Returns an iterator over the segments in the graph.
    ///
    /// If there are gaps in the node id space of the graph, segments corresponding to unused ids may be empty.
    pub fn segment_iter(&self) -> SegmentIter {
        let mut iter = self.mapping.one_iter();
        let first_node = if let Some((_, value)) = iter.next() { value } else { 0 };
        let next = (0, first_node);
        let limit = (self.segments(), self.mapping.len());
        SegmentIter {
            parent: self,
            iter: iter,
            next: next,
            limit: limit,
        }
    }
}

//-----------------------------------------------------------------------------

impl Serialize for Graph {
    fn serialize_header<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.header.serialize(writer)
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.sequences.serialize(writer)?;
        self.segments.serialize(writer)?;
        self.mapping.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let header = Header::<GraphPayload>::load(reader)?;
        if let Err(msg) = header.validate() {
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }

        let sequences = StringArray::load(reader)?;

        let segments = StringArray::load(reader)?;
        if header.is_set(GraphPayload::FLAG_TRANSLATION) == segments.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Graph: Translation flag does not match the presence of segment names"));
        }

        let mapping = SparseVector::load(reader)?;
        if header.is_set(GraphPayload::FLAG_TRANSLATION) {
            if mapping.len() != header.payload().nodes + 1 {
                return Err(Error::new(ErrorKind::InvalidData, "Graph: Node-to-segment mapping does not match the number of nodes"));
            }
            if mapping.count_ones() != sequences.len() {
                return Err(Error::new(ErrorKind::InvalidData, "Graph: Node-to-segment mapping does not match the number of segments"));
            }
        } else if !segments.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Graph: Translation flag does not match the presence of node-to-segment mapping"));
        }

        Ok(Graph {
            header: header,
            sequences: sequences,
            segments: segments,
            mapping: mapping,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.header.size_in_elements() + self.sequences.size_in_elements() + self.segments.size_in_elements() + self.mapping.size_in_elements()
    }
}

//-----------------------------------------------------------------------------

// FIXME document, example, tests
// FIXME item is (name, node range, sequence)
#[derive(Clone, Debug)]
pub struct SegmentIter<'a> {
    parent: &'a Graph,
    // Iterator over `parent.mapping`.
    iter: OneIter<'a>,
    // The first (segment, node) identifier we have not used.
    next: (usize, usize),
    // The first (segment, node) identifier we should not use.
    limit: (usize, usize),
}

impl<'a> Iterator for SegmentIter<'a> {
    type Item = (&'a [u8], Range<usize>, &'a [u8]);

    fn next(&mut self) -> Option<Self::Item> {
        if self.next.0 >= self.limit.0 {
            None
        } else {
            let name = self.parent.segments.bytes(self.next.0);
            self.next.0 += 1;
            let prev = self.next.1;
            self.next.1 = if let Some((_, value)) = self.iter.next() { value } else { self.limit.1 };
            // Translate node range to sequence range.
            let sequence = self.parent.sequences.range(prev - 1..self.next.1 - 1);
            Some((name, prev..self.next.1, sequence))
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.limit.0 - self.next.0;
        (remaining, Some(remaining))
    }
}

impl<'a> DoubleEndedIterator for SegmentIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.next.0 >= self.limit.0 {
            None
        } else {
            self.limit.0 -= 1;
            let name = self.parent.segments.bytes(self.limit.0);
            let old_node_limit = self.limit.1;
            self.limit.1 = if let Some((_, value)) = self.iter.next_back() { value } else { self.next.1 };
            // Translate node range to sequence range.
            let sequence = self.parent.sequences.range(self.limit.1 - 1..old_node_limit - 1);
            Some((name, self.limit.1..old_node_limit, sequence))
        }
    }
}

impl<'a> ExactSizeIterator for SegmentIter<'a> {}

impl<'a> FusedIterator for SegmentIter<'a> {}

//-----------------------------------------------------------------------------
