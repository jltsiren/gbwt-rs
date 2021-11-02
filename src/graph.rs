//! GBWTGraph: Node sequences and node-to-segment translation.
//!
//! The [`Graph`] structure augments a [`crate::GBWT`] index with node sequences and an optional node-to-segment translation.
//! This enables representing bidirected sequence graphs compatible with a subset of the [GFA format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).
//! Unlike in the [C++ implementation](https://github.com/jltsiren/gbwtgraph), the actual graph interface is provided by the [`crate::GBZ`] structure.
//!
//! At the moment, this implementation only supports graphs built with other tools.

use crate::headers::{Header, GraphPayload};
use crate::support::{StringArray, StringIter};

use simple_sds::ops::{BitVec, Select, PredSucc};
use simple_sds::serialize::Serialize;
use simple_sds::sparse_vector::{SparseVector, OneIter};

use std::io::{Error, ErrorKind};
use std::iter::FusedIterator;
use std::ops::Range;
use std::io;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Node sequences and an optional node-to-segment translation for a graph based on a [`crate::GBWT`] index.
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
/// ## Sequences
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
///
/// ## Segments
///
/// ```
/// use gbwt::{Graph, Segment};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gg");
/// let graph: Graph = serialize::load_from(&filename).unwrap();
///
/// assert!(graph.has_translation());
/// assert_eq!(graph.segments(), 7);
///
/// let first = Segment::from_fields(0, "s11".as_bytes(), 1..3, "GAT".as_bytes());
/// assert_eq!(graph.segment(0), first);
///
/// let middle = Segment::from_fields(3, "s14".as_bytes(), 5..7, "CAG".as_bytes());
/// assert_eq!(graph.node_to_segment(6), middle);
///
/// let last = Segment::from_fields(6, "s17".as_bytes(), 9..10, "TA".as_bytes());
/// assert_eq!(graph.segment_name(6), last.name);
/// assert_eq!(graph.segment_nodes(6), last.nodes);
/// assert_eq!(graph.segment_sequence(6), last.sequence);
/// assert_eq!(graph.segment_len(6), last.sequence.len());
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

    /// Returns the segment with the given identifier.
    ///
    /// Note that random access to the node-to-segment translation is somewhat slow.
    ///
    /// # Panics
    ///
    /// May panic if `id >= self.segments()`.
    pub fn segment(&self, id: usize) -> Segment {
        let name = self.segment_name(id);
        let nodes = self.segment_nodes(id);
        let sequence = self.sequences.range(nodes.start - 1..nodes.end - 1);
        Segment::from_fields(id, name, nodes, sequence)
    }

    /// Returns the segment containing node `node_id` of the original graph.
    ///
    /// Note that random access to the node-to-segment translation is somewhat slow.
    ///
    /// # Panics
    ///
    /// May panic if `node_id == 0` or `node_id > self.sequences()` or if there is no node-to-segment translation.
    pub fn node_to_segment(&self, node_id: usize) -> Segment {
        let mut iter = self.mapping.predecessor(node_id);
        let (segment_id, start) = iter.next().unwrap();
        let end = if let Some((_, value)) = iter.next() { value } else { self.mapping.len() };
        let name = self.segment_name(segment_id);
        let sequence = self.sequences.range(start - 1..end - 1);
        Segment {
            id: segment_id,
            name: name,
            nodes: start..end,
            sequence: sequence,
        }
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
        let end = if let Some((_, value)) = iter.next() { value } else { self.mapping.len() };
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
            if mapping.len() != sequences.len() + 1 {
                return Err(Error::new(ErrorKind::InvalidData, "Graph: Node-to-segment mapping does not match the number of sequences"));
            }
            if mapping.count_ones() != segments.len() {
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

/// A segment in a GFA graph.
///
/// A segment has an integer identifier, a string name, and a stored sequence.
/// It segment corresponds to a concatenation of nodes with consecutive identifiers in the original graph.
/// The name and the sequence are assumed to be valid UTF-8, but this is not checked for performance reasons.
/// If the nodes are not used on any path in the original graph, the name and the sequence may be empty.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Segment<'a> {
    /// Identifier of the segment.
    pub id: usize,
    /// Name of the segment.
    pub name: &'a [u8],
    /// The corresponding range of node ids in the original graph.
    pub nodes: Range<usize>,
    /// Sequence corresponding to the segment.
    pub sequence: &'a [u8],
}

impl<'a> Segment<'a> {
    /// Builds a segment from the given fields.
    ///
    /// Arguments corresponds to the fields with the same name.
    pub fn from_fields(id: usize, name: &'a [u8], nodes: Range<usize>, sequence: &'a [u8]) -> Self {
        Segment {
            id: id,
            name: name,
            nodes: nodes,
            sequence: sequence,
        }
    }
}

//-----------------------------------------------------------------------------

/// A read-only iterator over the segments in a [`Graph`].
///
/// The type of `Item` is [`Segment`].
/// If there are gaps in the node id space of the graph, segments corresponding to unused ids may be empty.
///
/// # Examples
///
/// ```
/// use gbwt::{Graph, Segment};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gg");
/// let graph: Graph = serialize::load_from(&filename).unwrap();
///
/// assert!(graph.has_translation());
/// let mut iter = graph.segment_iter();
/// let first = Segment::from_fields(0, "s11".as_bytes(), 1..3, "GAT".as_bytes());
/// assert_eq!(iter.next(), Some(first));
/// let last = Segment::from_fields(6, "s17".as_bytes(), 9..10, "TA".as_bytes());
/// assert_eq!(iter.next_back(), Some(last));
/// ```
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
    type Item = Segment<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next.0 >= self.limit.0 {
            None
        } else {
            let id = self.next.0;
            self.next.0 += 1;
            let name = self.parent.segments.bytes(id);
            let prev = self.next.1;
            self.next.1 = if let Some((_, value)) = self.iter.next() { value } else { self.limit.1 };
            // Translate node range to sequence range.
            let sequence = self.parent.sequences.range(prev - 1..self.next.1 - 1);
            Some(Segment::from_fields(id, name, prev..self.next.1, sequence))
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
            let id = self.limit.0;
            let name = self.parent.segments.bytes(id);
            let old_node_limit = self.limit.1;
            self.limit.1 = if let Some((_, value)) = self.iter.next_back() { value } else { self.next.1 };
            // Translate node range to sequence range.
            let sequence = self.parent.sequences.range(self.limit.1 - 1..old_node_limit - 1);
            Some(Segment::from_fields(id, name, self.limit.1..old_node_limit, sequence))
        }
    }
}

impl<'a> ExactSizeIterator for SegmentIter<'a> {}

impl<'a> FusedIterator for SegmentIter<'a> {}

//-----------------------------------------------------------------------------
