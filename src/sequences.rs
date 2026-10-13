//! GBWTGraph: Node sequences and node-to-segment translation.
//!
//! The [`Sequences`] structure augments a [`crate::GBWT`] index with node sequences and an optional node-to-segment translation.
//! This enables representing bidirected sequence graphs compatible with a subset of the [GFA format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).
//! Unlike in the [C++ implementation](https://github.com/jltsiren/gbwtgraph), the actual graph interface is provided by the [`crate::GBZ`] structure.
//!
//! A [`Sequences`] structure can be built from the parts of a graph with
//! [`Sequences::from_graph`].

use crate::{ENDMARKER, GBWT, Orientation};
use crate::graph_builder::GraphBuilder;
use crate::headers::{Header, Payload, SequencesPayload};
use crate::support::{self, StringArray, StringIter};

use simple_sds::bit_vector::BitVector;
use simple_sds::ops::{BitVec, Select, PredSucc};
use simple_sds::raw_vector::{RawVector, AccessRaw};
use simple_sds::serialize::{Serialize, SerializeVersion};
use simple_sds::sparse_vector::{SparseVector, SparseBuilder, OneIter};

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
/// use gbz::sequences::Sequences;
/// use gbz::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gg");
/// let sequences: Sequences = serialize::load_from(&filename).unwrap();
///
/// assert_eq!(sequences.nodes(), 12);
/// assert_eq!(sequences.sequences(), 15);
/// assert_eq!(sequences.sequence(13), "T".as_bytes());
/// assert_eq!(sequences.sequence_len(3), 1);
///
/// let labels: Vec<&[u8]> = sequences.iter().collect();
/// assert_eq!(labels.concat(), "GATTACAGATTA".as_bytes());
/// ```
///
/// ## Segments
///
/// ```
/// use gbz::Segment;
/// use gbz::sequences::Sequences;
/// use gbz::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gg");
/// let sequences: Sequences = serialize::load_from(&filename).unwrap();
///
/// assert!(sequences.has_translation());
/// assert_eq!(sequences.segments(), 8);
///
/// let first = Segment::from_fields(0, "s11".as_bytes(), 1..3, "GAT".as_bytes());
/// assert_eq!(sequences.segment(0), first);
///
/// let middle = Segment::from_fields(3, "s14".as_bytes(), 5..7, "CAG".as_bytes());
/// assert_eq!(sequences.node_to_segment(6), middle);
///
/// let last = Segment::from_fields(7, "s17".as_bytes(), 11..12, "TA".as_bytes());
/// assert_eq!(sequences.segment_name(7), last.name);
/// assert_eq!(sequences.segment_nodes(7), last.nodes);
/// assert_eq!(sequences.segment_sequence(7), last.sequence);
/// assert_eq!(sequences.segment_len(7), last.sequence.len());
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequences {
    header: Header<SequencesPayload>,
    sequences: StringArray,
    segments: StringArray,
    mapping: SparseVector,
}

//-----------------------------------------------------------------------------

/// Sequences for nodes.
impl Sequences {
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
    pub fn iter(&'_ self) -> StringIter<'_> {
        self.sequences.iter()
    }
}

//-----------------------------------------------------------------------------

/// Segments.
impl Sequences {
    /// Returns `true` if the graph contains a node-to-segment translation.
    #[inline]
    pub fn has_translation(&self) -> bool {
        self.header.is_set(SequencesPayload::FLAG_TRANSLATION)
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
    pub fn segment(&'_ self, id: usize) -> Segment<'_> {
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
    pub fn node_to_segment(&'_ self, node_id: usize) -> Segment<'_> {
        let mut iter = self.mapping.predecessor(node_id);
        let (segment_id, start) = iter.next().unwrap();
        let end = if let Some((_, value)) = iter.next() { value } else { self.mapping.len() };
        let name = self.segment_name(segment_id);
        let sequence = self.sequences.range(start - 1..end - 1);
        Segment {
            id: segment_id,
            name,
            nodes: start..end,
            sequence,
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
    /// See also [`Sequences::segment_range`].
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
    /// See also [`Sequences::segment_nodes`].
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
    pub fn segment_iter(&'_ self) -> SegmentIter<'_> {
        let mut iter = self.mapping.one_iter();
        let first_node = if let Some((_, value)) = iter.next() { value } else { 0 };
        let next = (0, first_node);
        let limit = (self.segments(), self.mapping.len());
        SegmentIter {
            parent: self,
            iter,
            next,
            limit,
        }
    }
}

//-----------------------------------------------------------------------------

/// Construction.
impl Sequences {
    /// Builds a [`Sequences`] structure from a [`GraphBuilder`] and a
    /// bidirectional [`GBWT`] index over the same graph.
    ///
    /// The node sequences are stored for every potential original node in the
    /// interval determined by the GBWT alphabet. Nodes that are not visited by
    /// any path in the GBWT (and hence have an empty record) get an empty
    /// sequence. Nodes of the graph outside the interval are not stored.
    ///
    /// If the graph contains a node-to-segment translation, segments whose
    /// first node is not present in the GBWT get an empty name. Segments after
    /// the last node of the GBWT alphabet are not stored. As the serialized
    /// format assumes that node identifiers start from `1` when a translation
    /// is in use, the first node of the GBWT must be node `1` in that case.
    ///
    /// # Arguments
    ///
    /// * `gbwt`: A bidirectional GBWT index over the graph.
    /// * `graph`: The graph used for building the GBWT.
    ///
    /// # Errors
    ///
    /// Returns an error if the GBWT is not bidirectional, if a node visited by
    /// the GBWT is missing from the graph, or if the graph uses a translation
    /// and the first node of the GBWT is not node `1`.
    /// An empty GBWT produces an empty [`Sequences`] structure.
    pub fn from_graph(gbwt: &GBWT, graph: &GraphBuilder) -> Result<Sequences, String> {
        if !gbwt.is_bidirectional() {
            return Err(String::from("Sequences: The GBWT index must be bidirectional"));
        }
        if gbwt.is_empty() {
            // Like in the C++ implementation, an empty GBWT gives an empty graph.
            let mapping = SparseVector::try_from(SparseBuilder::new(0, 0).map_err(|e| format!("Sequences: {}", e))?)
                .map_err(|e| format!("Sequences: {}", e))?;
            return Ok(Sequences {
                header: Header::<SequencesPayload>::default(),
                sequences: StringArray::from_bytes(&[]),
                segments: StringArray::from_bytes(&[]),
                mapping,
            });
        }

        let min_node = support::node_id(gbwt.first_node());
        let max_node = support::node_id(gbwt.alphabet_size()) - 1;
        if graph.uses_translation() && min_node != 1 {
            return Err(String::from("Sequences: A graph with a translation must start from node 1 in the GBWT"));
        }

        // Collect the sequences of the present nodes without copying them.
        let present = present_nodes(gbwt);
        let mut header = Header::<SequencesPayload>::default();
        let mut node_sequences: Vec<&[u8]> = Vec::with_capacity(present.len());
        let mut nodes = 0;
        for sequence_id in 0..present.len() {
            if present.get(sequence_id) {
                let node_id = min_node + sequence_id;
                if !graph.has_node(node_id) {
                    return Err(format!("Sequences: Node {} is used in the GBWT but missing from the graph", node_id));
                }
                node_sequences.push(graph.node_sequence(node_id));
                nodes += 1;
            } else {
                node_sequences.push(&[]);
            }
        }
        header.payload_mut().nodes = nodes;

        // Build the node-to-segment translation.
        let mut segments: Vec<&[u8]> = Vec::new();
        let mut first_nodes: Vec<usize> = Vec::new();
        if graph.uses_translation() {
            header.set(SequencesPayload::FLAG_TRANSLATION);
            for (name, range) in graph.segment_iter() {
                // Segments are listed in ascending node id order, so we can stop
                // when we reach a segment outside the GBWT node id space.
                if range.start > max_node {
                    break;
                }
                first_nodes.push(range.start);
                if present.get(range.start - min_node) {
                    segments.push(name.as_bytes());
                } else {
                    segments.push(&[]);
                }
            }
        }

        // With a translation, `min_node == 1`, so the universe is `sequences + 1`.
        let universe = if graph.uses_translation() { max_node + 1 } else { 0 };
        let mut mapping_builder = SparseBuilder::new(universe, first_nodes.len())
            .map_err(|e| format!("Sequences: {}", e))?;
        for first in first_nodes {
            mapping_builder.set(first);
        }
        let mapping = SparseVector::try_from(mapping_builder)
            .map_err(|e| format!("Sequences: {}", e))?;

        Ok(Sequences {
            header,
            sequences: StringArray::from_bytes(&node_sequences),
            segments: StringArray::from_bytes(&segments),
            mapping,
        })
    }
}

// Returns a bitvector marking the original nodes visited by the GBWT.
//
// The bit at offset `i` corresponds to node `i + min_node`, where `min_node` is the
// first original node in the GBWT alphabet. The bitvector is empty if the GBWT is empty.
pub(crate) fn present_nodes(gbwt: &GBWT) -> BitVector {
    let potential_nodes = gbwt.alphabet_size().saturating_sub(gbwt.first_node()) / 2;
    let mut present = RawVector::with_len(potential_nodes, false);
    for record_id in gbwt.as_ref().id_iter() {
        if record_id == ENDMARKER {
            continue;
        }
        let gbwt_node = gbwt.record_to_node(record_id);
        if support::node_orientation(gbwt_node) == Orientation::Forward {
            present.set_bit((gbwt_node - gbwt.first_node()) / 2, true);
        }
    }
    BitVector::from(present)
}

//-----------------------------------------------------------------------------

impl Serialize for Sequences {
    fn serialize_header<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.serialize_header_version(writer, Self::DEFAULT_VERSION)
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.serialize_body_version(writer, Self::DEFAULT_VERSION)
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let mut header = Header::<SequencesPayload>::load(reader)?;
        if let Err(msg) = header.validate() {
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }

        // Decompress or deserialize the sequences.
        let sequences = if header.version() >= SequencesPayload::ZSTD_VERSION {
            StringArray::decompress(reader)?
        } else {
            StringArray::load(reader)?
        };

        let segments = StringArray::load(reader)?;
        if header.is_set(SequencesPayload::FLAG_TRANSLATION) == segments.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Sequences: Translation flag does not match the presence of segment names"));
        }

        let mapping = SparseVector::load(reader)?;
        if header.is_set(SequencesPayload::FLAG_TRANSLATION) {
            // If there are no gaps in the node id space, `mapping.len() == header.payload().nodes + 1`.
            // Unused nodes create gaps.
            if mapping.len() <= header.payload().nodes {
                return Err(Error::new(ErrorKind::InvalidData, "Sequences: Node-to-segment mapping does not match the number of nodes"));
            }
            if mapping.len() != sequences.len() + 1 {
                return Err(Error::new(ErrorKind::InvalidData, "Sequences: Node-to-segment mapping does not match the number of sequences"));
            }
            if mapping.count_ones() != segments.len() {
                return Err(Error::new(ErrorKind::InvalidData, "Sequences: Node-to-segment mapping does not match the number of segments"));
            }
        } else if !segments.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Sequences: Translation flag does not match the presence of node-to-segment mapping"));
        }

        // Update the header to the latest version after we have used the
        // serialized version for loading the correct data.
        header.update();

        Ok(Sequences {
            header, sequences, segments, mapping,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.header.size_in_elements() +
            self.sequences.compressed_size_in_elements(None) +
            self.segments.size_in_elements() +
            self.mapping.size_in_elements()
    }
}

impl SerializeVersion for Sequences {
    // This is currently the same as `SequencesPayload::MIN_VERSION`, but we could
    // plausibly be able to read versions we cannot serialize.
    const MIN_VERSION: usize = 3;

    // We should always be able to serialize the latest version.
    const MAX_VERSION: usize = SequencesPayload::VERSION as usize;

    // Unlike the C++ implementation, we have no reason to default to an older
    // version for backwards compatibility.
    const DEFAULT_VERSION: usize = SequencesPayload::VERSION as usize;

    fn serialize_header_version<T: io::Write>(&self, writer: &mut T, version: usize) -> io::Result<()> {
        Self::ensure_supported_version(version, "Sequences")?;
        let mut copy = self.header;
        copy.update_to_version(version as u32)?;
        copy.serialize(writer)
    }

    fn serialize_body_version<T: io::Write>(&self, writer: &mut T, version: usize) -> io::Result<()> {
        Self::ensure_supported_version(version, "Sequences")?;
        if version as u32 >= SequencesPayload::ZSTD_VERSION {
            self.sequences.compress(writer, None)?;
        } else {
            self.sequences.serialize(writer)?;
        }
        self.segments.serialize(writer)?;
        self.mapping.serialize(writer)?;
        Ok(())
    }

    fn determine_version<T: io::Read>(reader: &mut T) -> io::Result<usize> {
        let header = Header::<SequencesPayload>::load(reader)?;
        Ok(header.version() as usize)
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
            id, name, nodes, sequence,
        }
    }
}

//-----------------------------------------------------------------------------

/// A read-only iterator over the segments in a [`Sequences`].
///
/// The type of `Item` is [`Segment`].
/// If there are gaps in the node id space of the graph, segments corresponding to unused ids may be empty.
///
/// # Examples
///
/// ```
/// use gbz::Segment;
/// use gbz::sequences::Sequences;
/// use gbz::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gg");
/// let sequences: Sequences = serialize::load_from(&filename).unwrap();
///
/// assert!(sequences.has_translation());
/// let mut iter = sequences.segment_iter();
/// let first = Segment::from_fields(0, "s11".as_bytes(), 1..3, "GAT".as_bytes());
/// assert_eq!(iter.next(), Some(first));
/// let last = Segment::from_fields(7, "s17".as_bytes(), 11..12, "TA".as_bytes());
/// assert_eq!(iter.next_back(), Some(last));
/// ```
#[derive(Clone, Debug)]
pub struct SegmentIter<'a> {
    parent: &'a Sequences,
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
