//! GBZ: Space-efficient representation for a subset of GFA.
//!
//! The [`GBZ`] structure combines a [`GBWT`] index and the sequences in a [`Graph`] structure.
//! It is a space-efficient file format for storing a pangenome graph based on aligned haplotype sequences, with the haplotypes stored as paths.
//! The format is compatible with a subset of [GFA1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).
//!
//! See also the [C++ implementation](https://github.com/jltsiren/gbwtgraph) and the [file format specification](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md).

use crate::{ENDMARKER, SOURCE_KEY, SOURCE_VALUE};
use crate::{Graph, Segment, GBWT, Orientation};
use crate::bwt::Record;
use crate::gbwt::{SequenceIter, Metadata};
use crate::graph::SegmentIter as GraphSegmentIter;
use crate::headers::{Header, GBZPayload};
use crate::support::Tags;
use crate::support;

use simple_sds::bit_vector::{BitVector, OneIter, Identity};
use simple_sds::ops::{BitVec, Select};
use simple_sds::raw_vector::{RawVector, AccessRaw};
use simple_sds::serialize::Serialize;

use std::io::{Error, ErrorKind};
use std::iter::FusedIterator;
use std::ops::Range;
use std::io;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// GBZ file format for storing GFA graphs with many paths.
///
/// A GBZ graph combines a [`GBWT`] index and a [`Graph`].
/// It represents the subgraph of the original graph induced by the paths in the GBWT index.
/// `GBZ` methods used node identifiers in the original graph.
/// Functions [`support::encode_node`], [`support::node_id`], [`support::node_orientation`], [`support::decode_node`], and [`support::flip_node`] enable conversions between original and GBWT node ids.
/// While the methods in [`Graph`] may panic or return unpredictable results when the correspondin object does not exist, `GBZ` methods return [`None`] in such cases.
///
/// In addition to representing a bidirected sequence graph with integer node ids, a `GBZ` also contains an optional node-to-segment translation.
/// The translation enables representing a GFA graph, where each segment has a string name and corresponds to the concatenation of a range of nodes.
/// Nodes and edges refer to those in the original graph, while GFA graphs have segments and links.
///
/// # Examples
///
/// See also: [`NodeIter`], [`EdgeIter`], [`SegmentIter`], [`LinkIter`], [`PathIter`], [`SegmentPathIter`]
///
/// ## Nodes
///
/// ```
/// use gbwt::GBZ;
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// assert_eq!(gbz.nodes(), 12);
/// assert_eq!(gbz.min_node(), 11);
/// assert_eq!(gbz.max_node(), 25);
///
/// assert!(gbz.has_node(13));
/// assert_eq!(gbz.sequence(13).unwrap(), "T".as_bytes());
/// assert_eq!(gbz.sequence_len(13), Some(1));
/// ```
///
/// ## Segments
///
/// ```
/// use gbwt::{GBZ, Segment};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// assert!(gbz.has_translation());
///
/// let middle = Segment::from_fields(3, "s14".as_bytes(), 5..7, "CAG".as_bytes());
/// assert_eq!(gbz.node_to_segment(6), Some(middle));
/// ```
///
/// ## Paths
///
/// ```
/// use gbwt::{GBZ, Orientation};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// assert_eq!(gbz.paths(), 6);
/// let mut iter = gbz.path(3, Orientation::Forward).unwrap();
/// assert_eq!(iter.next(), Some((11, Orientation::Forward)));
/// assert_eq!(iter.nth(3), Some((17, Orientation::Forward)));
/// assert!(iter.next().is_none());
///
/// assert!(gbz.has_metadata());
/// let metadata = gbz.metadata().unwrap();
/// let path_name = metadata.path(3).unwrap();
/// assert_eq!(metadata.sample(path_name.sample()), Some("sample"));
/// assert_eq!(metadata.contig(path_name.contig()), Some("A"));
/// ```
///
/// # Notes
///
/// * Methods may panic if given path / node identifiers `id > usize::MAX / 2`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GBZ {
    header: Header<GBZPayload>,
    tags: Tags,
    index: GBWT,
    graph: Graph,
    real_nodes: BitVector,
}

//-----------------------------------------------------------------------------

// Private utilities.
impl GBZ {
    // Converts a node id in the original graph to a sequence id.
    #[inline]
    fn graph_node_to_sequence(&self, node_id: usize) -> usize {
        let gbwt_node = support::encode_node(node_id, Orientation::Forward);
        Self::gbwt_node_to_sequence(&self.index, gbwt_node)
    }

    // Converts a GBWT node id to a sequence id.
    #[inline]
    fn gbwt_node_to_sequence(index: &GBWT, gbwt_node: usize) -> usize {
        (gbwt_node - index.first_node()) / 2
    }

    // Converts a sequence id to a node id in the original graph.
    #[inline]
    fn sequence_to_graph_node(&self, sequence_id: usize) -> usize {
        support::node_id(sequence_id * 2 + self.index.first_node())
    }
}

/// Nodes and edges.
impl GBZ {
    /// Returns the number of nodes in the graph.
    #[inline]
    pub fn nodes(&self) -> usize {
        self.graph.nodes()
    }

    /// Returns the smallest node identifier in the original graph.
    #[inline]
    pub fn min_node(&self) -> usize {
        support::node_id(self.index.first_node())
    }

    /// Returns the largest node identifier in the original graph.
    #[inline]
    pub fn max_node(&self) -> usize {
        support::node_id(self.index.alphabet_size()) - 1
    }

    /// Returns `true` if the original graph contains a node with identifier `node_id`.
    #[inline]
    pub fn has_node(&self, node_id: usize) -> bool {
        let gbwt_node = support::encode_node(node_id, Orientation::Forward);
        self.index.has_node(gbwt_node) && self.real_nodes.get(Self::gbwt_node_to_sequence(&self.index, gbwt_node))
    }

    /// Returns the sequence for the node in the original graph, or [`None`] if there is no such node.
    pub fn sequence(&self, node_id: usize) -> Option<&[u8]> {
        if !self.has_node(node_id) {
            return None;
        }
        let sequence_id = self.graph_node_to_sequence(node_id);
        Some(self.graph.sequence(sequence_id))
    }

    /// Returns the length of the sequence for the node in the original graph, or [`None`] if there is no such node.
    pub fn sequence_len(&self, node_id: usize) -> Option<usize> {
        if !self.has_node(node_id) {
            return None;
        }
        let sequence_id = self.graph_node_to_sequence(node_id);
        Some(self.graph.sequence_len(sequence_id))
    }

    /// Returns an iterator over the node identifiers in the original graph.
    ///
    /// See [`NodeIter`] for an example.
    pub fn node_iter(&self) -> NodeIter {
        NodeIter {
            parent: &self,
            iter: self.real_nodes.one_iter()
        }
    }

    /// Returns an iterator over the successors of a node, or [`None`] if there is no such node.
    ///
    /// See [`EdgeIter`] for an example.
    ///
    /// # Arguments
    ///
    /// * `node_id`: Identifier of the node.
    /// * `orientation`: Orientation of the node.
    pub fn successors(&self, node_id: usize, orientation: Orientation) -> Option<EdgeIter> {
        if !self.has_node(node_id) {
            return None;
        }
        let gbwt_node = support::encode_node(node_id, orientation);
        let record_id = self.index.node_to_record(gbwt_node);
        let record = self.index.as_ref().record(record_id)?;
        let next = if record.outdegree() > 0 && record.successor(0) == ENDMARKER { 1 } else { 0 };
        let limit = record.outdegree();
        Some(EdgeIter {
            record: record,
            next: next,
            limit: limit,
            flip: false,
        })
    }

    /// Returns an iterator over the predecessors of a node, or [`None`] if there is no such node.
    ///
    /// See [`EdgeIter`] for an example.
    ///
    /// # Arguments
    ///
    /// * `node_id`: Identifier of the node.
    /// * `orientation`: Orientation of the node.
    pub fn predecessors(&self, node_id: usize, orientation: Orientation) -> Option<EdgeIter> {
        if !self.has_node(node_id) {
            return None;
        }
        let gbwt_node = support::encode_node(node_id, orientation.flip());
        let record_id = self.index.node_to_record(gbwt_node);
        let record = self.index.as_ref().record(record_id)?;
        let next = if record.outdegree() > 0 && record.successor(0) == ENDMARKER { 1 } else { 0 };
        let limit = record.outdegree();
        Some(EdgeIter {
            record: record,
            next: next,
            limit: limit,
            flip: true,
        })
    }
}

//-----------------------------------------------------------------------------

/// Segments and links.
impl GBZ {
    /// Returns `true` if the graph contains a node-to-segment translation.
    #[inline]
    pub fn has_translation(&self) -> bool {
        self.graph.has_translation()
    }

    /// Returns the segment containing node `node_id` of the original graph.
    ///
    /// Returns [`None`] if there is no such node or if the graph does not contain a node-to-segment translation.
    /// Note that random access to the translation is somewhat slow.
    pub fn node_to_segment(&self, node_id: usize) -> Option<Segment> {
        if self.has_translation() && self.has_node(node_id) {
            Some(self.graph.node_to_segment(node_id))
        } else {
            None
        }
    }

    /// Returns an iterator over the segments in the GFA graph, or [`None`] if there is no node-to-segment translation.
    ///
    /// See [`SegmentIter`] for an example.
    pub fn segment_iter(&self) -> Option<SegmentIter> {
        if self.has_translation() {
            Some(SegmentIter {
                parent: self,
                iter: self.graph.segment_iter(),
            })
        } else {
            None
        }
    }

    /// Returns an iterator over the successors of a segment.
    ///
    /// Returns [`None`] if there is no node-to-segment-translation or if the nodes corresponding to the segment do not exist.
    /// The result is unpredictable if `segment` is not a valid segment.
    /// See [`LinkIter`] for an example.
    ///
    /// # Arguments
    ///
    /// * `segment`: A valid segment in the GFA graph.
    /// * `orientation`: Orientation of the segment.
    pub fn segment_successors(&self, segment: &Segment, orientation: Orientation) -> Option<LinkIter> {
        if segment.nodes.is_empty() || !self.has_translation() {
            return None;
        }
        let node_id = match orientation {
            Orientation::Forward => segment.nodes.end - 1,
            Orientation::Reverse => segment.nodes.start,
        };
        let iter = self.successors(node_id, orientation)?;
        Some(LinkIter {
            parent: self,
            iter: iter,
        })
    }

    /// Returns an iterator over the predecessors of a segment.
    ///
    /// Returns [`None`] if there is no node-to-segment-translation or if the nodes corresponding to the segment do not exist.
    /// The result is unpredictable if `segment` is not a valid segment.
    /// See [`LinkIter`] for an example.
    ///
    /// # Arguments
    ///
    /// * `segment`: A valid segment in the GFA graph.
    /// * `orientation`: Orientation of the segment.
    pub fn segment_predecessors(&self, segment: &Segment, orientation: Orientation) -> Option<LinkIter> {
        if segment.nodes.is_empty() || !self.has_translation() {
            return None;
        }
        let node_id = match orientation {
            Orientation::Forward => segment.nodes.start,
            Orientation::Reverse => segment.nodes.end - 1,
        };
        let iter = self.predecessors(node_id, orientation)?;
        Some(LinkIter {
            parent: self,
            iter: iter,
        })
    }
}

//-----------------------------------------------------------------------------

// FIXME follow_forward, follow_backward
/// Paths.
impl GBZ {
    /// Returns the number of paths in the original graph.
    #[inline]
    pub fn paths(&self) -> usize {
        self.index.sequences() / 2
    }

    /// Returns an iterator over the given path, or [`None`] if there is no such path.
    ///
    /// See [`PathIter`] for an example.
    ///
    /// # Arguments
    ///
    /// * `path_id`: Path identifier in the original graph.
    /// * `orientation`: Orientation of the path.
    pub fn path(&self, path_id: usize, orientation: Orientation) -> Option<PathIter> {
        let iter = self.index.sequence(support::encode_path(path_id, orientation))?;
        Some(PathIter {
            iter: iter,
        })
    }

    /// Returns an iterator over the segments in the given path.
    ///
    /// The return value is [`None`] if there is no such path or there is no node-to-segment translation.
    /// See [`SegmentPathIter`] for an example.
    ///
    /// # Arguments
    ///
    /// * `path_id`: Path identifier in the original graph.
    /// * `orientation`: Orientation of the path.
    pub fn segment_path(&self, path_id: usize, orientation: Orientation) -> Option<SegmentPathIter> {
        if !self.has_translation() {
            return None;
        }
        let iter = self.index.sequence(support::encode_path(path_id, orientation))?;
        Some(SegmentPathIter {
            parent: self,
            iter: iter,
            next: None,
            nodes: 0..0,
            fail: false,
        })
    }

    /// Returns `true` if the GBWT index contains metadata.
    pub fn has_metadata(&self) -> bool {
        self.index.has_metadata()
    }

    /// Returns the metadata stored in the GBWT index, or [`None`] if there is no metadata.
    pub fn metadata(&self) -> Option<&Metadata> {
        self.index.metadata()
    }
}

//-----------------------------------------------------------------------------

impl Serialize for GBZ {
    fn serialize_header<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.header.serialize(writer)
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.tags.serialize(writer)?;
        self.index.serialize(writer)?;
        self.graph.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let header = Header::<GBZPayload>::load(reader)?;
        if let Err(msg) = header.validate() {
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }

        let mut tags = Tags::load(reader)?;
        tags.insert(SOURCE_KEY, SOURCE_VALUE);

        let index = GBWT::load(reader)?;
        if !index.is_bidirectional() {
            return Err(Error::new(ErrorKind::InvalidData, "GBZ: The GBWT index is not bidirectional"));
        }
        let potential_nodes = (index.alphabet_size() - index.first_node()) / 2;

        let graph = Graph::load(reader)?;
        if graph.sequences() != potential_nodes {
            return Err(Error::new(ErrorKind::InvalidData, "GBZ: Mismatch between GBWT alphabet size and Graph sequence count"));
        }

        // Cache real nodes.
        let mut real_nodes = RawVector::with_len(potential_nodes, false);
        for record_id in index.as_ref().id_iter() {
            if record_id == ENDMARKER {
                continue;
            }
            let gbwt_node = index.record_to_node(record_id);
            if support::node_orientation(gbwt_node) == Orientation::Forward {
                real_nodes.set_bit(Self::gbwt_node_to_sequence(&index, gbwt_node), true);
            }
        }

        Ok(GBZ {
            header: header,
            tags: tags,
            index: index,
            graph: graph,
            real_nodes: BitVector::from(real_nodes),
        })
    }

    fn size_in_elements(&self) -> usize {
        self.header.size_in_elements() + self.tags.size_in_elements() + self.index.size_in_elements() + self.graph.size_in_elements()
    }
}

//-----------------------------------------------------------------------------

impl AsRef<GBWT> for GBZ {
    fn as_ref(&self) -> &GBWT {
        &self.index
    }
}

impl AsRef<Graph> for GBZ {
    fn as_ref(&self) -> &Graph {
        &self.graph
    }
}

//-----------------------------------------------------------------------------

/// An iterator over the node identifiers in the original graph.
///
/// Because the GBZ stores a subgraph induced by the paths, nodes that are not visited by any path will be skipped.
///
/// # Examples
///
/// ```
/// use gbwt::GBZ;
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// let nodes: Vec<usize> = gbz.node_iter().collect();
/// assert_eq!(nodes, vec![11, 12, 13, 14, 15, 16, 17, 21, 22, 23, 24, 25]);
/// ```
#[derive(Clone, Debug)]
pub struct NodeIter<'a> {
    parent: &'a GBZ,
    iter: OneIter<'a, Identity>,
}

impl<'a> Iterator for NodeIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(_, id)| self.parent.sequence_to_graph_node(id))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<'a> DoubleEndedIterator for NodeIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back().map(|(_, id)| self.parent.sequence_to_graph_node(id))
    }
}

impl<'a> ExactSizeIterator for NodeIter<'a> {}

impl<'a> FusedIterator for NodeIter<'a> {}

//-----------------------------------------------------------------------------

/// An iterator over the predecessors or successors of a node.
///
/// The type of `Item` is `(`[`usize`]`, `[`Orientation`]`)`.
/// These values encode the identifier and the orientation of the node.
/// Successor nodes are always listed in sorted order.
/// Predecessor nodes are sorted by their identifiers, but the reverse orientation is listed before the forward orientation.
///
/// # Examples
///
/// ```
/// use gbwt::{GBZ, Orientation};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// // Successors of node 24 in forward orientation.
/// assert!(gbz.has_node(24));
/// let successors: Vec<(usize, Orientation)> =
///     gbz.successors(24, Orientation::Forward).unwrap().collect();
/// assert_eq!(successors, vec![(23, Orientation::Reverse), (25, Orientation::Forward)]);
///
/// // Predecessors of node 14 in reverse orientation.
/// assert!(gbz.has_node(14));
/// let predecessors: Vec<(usize, Orientation)> =
///     gbz.predecessors(14, Orientation::Reverse).unwrap().collect();
/// assert_eq!(predecessors, vec![(15, Orientation::Reverse), (16, Orientation::Reverse)]);
/// ```
#[derive(Clone, Debug)]
pub struct EdgeIter<'a> {
    record: Record<'a>,
    // The first outrank that has not been visited.
    next: usize,
    // The first outrank that we should not visit.
    limit: usize,
    // Flip the orientation in the iterated values.
    flip: bool,
}

impl<'a> Iterator for EdgeIter<'a> {
    type Item = (usize, Orientation);

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.limit {
            None
        } else {
            let gbwt_node = self.record.successor(self.next);
            self.next += 1;
            // We use the reverse record when we iterate over the predecessors. If we flip a successor
            // of the reverse record, we get a predecessor of the forward record.
            let mut orientation = support::node_orientation(gbwt_node);
            if self.flip {
                orientation = orientation.flip();
            }
            Some((support::node_id(gbwt_node), orientation))
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.limit - self.next;
        (remaining, Some(remaining))
    }
}

impl<'a> DoubleEndedIterator for EdgeIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.next >= self.limit {
            None
        } else {
            self.limit -= 1;
            let gbwt_node = self.record.successor(self.limit);
            // We use the reverse record when we iterate over the predecessors. If we flip a successor
            // of the reverse record, we get a predecessor of the forward record.
            let mut orientation = support::node_orientation(gbwt_node);
            if self.flip {
                orientation = orientation.flip();
            }
            Some((support::node_id(gbwt_node), orientation))
        }
    }
}

impl<'a> ExactSizeIterator for EdgeIter<'a> {}

impl<'a> FusedIterator for EdgeIter<'a> {}

//-----------------------------------------------------------------------------

/// A read-only iterator over the segments in the GFA graph.
///
/// The type of `Item` is [`Segment`].
/// Unlike [`crate::graph::SegmentIter`], this iterator will skip segments corresponding to unused node ids.
///
/// # Examples
///
/// ```
/// use gbwt::{GBZ, Segment};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// assert!(gbz.has_translation());
/// let mut iter = gbz.segment_iter().unwrap();
/// let first = Segment::from_fields(0, "s11".as_bytes(), 1..3, "GAT".as_bytes());
/// assert_eq!(iter.next(), Some(first));
/// let last = Segment::from_fields(7, "s17".as_bytes(), 11..12, "TA".as_bytes());
/// assert_eq!(iter.next_back(), Some(last));
/// ```
#[derive(Clone, Debug)]
pub struct SegmentIter<'a> {
    parent: &'a GBZ,
    iter: GraphSegmentIter<'a>,
}

impl<'a> Iterator for SegmentIter<'a> {
    type Item = Segment<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(segment) = self.iter.next() {
            if self.parent.has_node(segment.nodes.start) {
                return Some(segment);
            }
        }
        None
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.iter.len()))
    }
}

impl<'a> DoubleEndedIterator for SegmentIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        while let Some(segment) = self.iter.next_back() {
            if self.parent.has_node(segment.nodes.start) {
                return Some(segment);
            }
        }
        None
    }
}

impl<'a> FusedIterator for SegmentIter<'a> {}

//-----------------------------------------------------------------------------

/// An iterator over the predecessors or successors of a segment.
///
/// The type of `Item` is `(`[`Segment`]`, `[`Orientation`]`)`.
/// Values are listed in the same order as by [`EdgeIter`].
///
/// # Examples
///
/// ```
/// use gbwt::{GBZ, Segment, Orientation};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// assert!(gbz.has_translation());
///
/// // Successors of segment s14 (nodes `5..7`) in forward orientation.
/// let s14 = gbz.node_to_segment(5).unwrap();
/// let mut successors = gbz.segment_successors(&s14, Orientation::Forward).unwrap();
/// assert_eq!(successors.next().map(|(s, o)| (s.name, o)),
///     Some(("s15".as_bytes(), Orientation::Forward)));
/// assert_eq!(successors.next().map(|(s, o)| (s.name, o)),
///     Some(("s16".as_bytes(), Orientation::Forward)));
/// assert!(successors.next().is_none());
///
/// // Predecessors of segment s11 (nodes `1..3`) in reverse orientation.
/// let s11 = gbz.node_to_segment(1).unwrap();
/// let mut predecessors = gbz.segment_predecessors(&s11, Orientation::Reverse).unwrap();
/// assert_eq!(predecessors.next().map(|(s, o)| (s.name, o)),
///     Some(("s12".as_bytes(), Orientation::Reverse)));
/// assert_eq!(predecessors.next().map(|(s, o)| (s.name, o)),
///     Some(("s13".as_bytes(), Orientation::Reverse)));
/// assert!(predecessors.next().is_none());
/// ```
#[derive(Clone, Debug)]
pub struct LinkIter<'a> {
    parent: &'a GBZ,
    iter: EdgeIter<'a>,
}

impl<'a> Iterator for LinkIter<'a> {
    type Item = (Segment<'a>, Orientation);

    fn next(&mut self) -> Option<Self::Item> {
        let (node_id, orientation) = self.iter.next()?;
        self.parent.node_to_segment(node_id).map(|s| (s, orientation))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<'a> DoubleEndedIterator for LinkIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let (node_id, orientation) = self.iter.next_back()?;
        self.parent.node_to_segment(node_id).map(|s| (s, orientation))
    }
}

impl<'a> ExactSizeIterator for LinkIter<'a> {}

impl<'a> FusedIterator for LinkIter<'a> {}

//-----------------------------------------------------------------------------

/// An iterator over an oriented path in the original graph.
///
/// The type of `Item` is `(`[`usize`]`, `[`Orientation`]`)`.
/// This gives the identifier and the orientation of a node in the original graph.
/// `PathIter` is a thin wrapper over [`SequenceIter`].
///
/// # Examples
///
/// ```
/// use gbwt::{GBZ, Orientation};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// // Path 3 in reverse orientation.
/// if let Some(mut iter) = gbz.path(3, Orientation::Reverse) {
///     assert_eq!(iter.next(), Some((17, Orientation::Reverse)));
///     assert_eq!(iter.next(), Some((16, Orientation::Reverse)));
///     assert_eq!(iter.next(), Some((14, Orientation::Reverse)));
///     assert_eq!(iter.next(), Some((13, Orientation::Reverse)));
///     assert_eq!(iter.next(), Some((11, Orientation::Reverse)));
///     assert!(iter.next().is_none());
/// } else {
///     panic!("No iterator for path 3 (reverse)");
/// }
/// ```
#[derive(Clone, Debug)]
pub struct PathIter<'a> {
    iter: SequenceIter<'a>,
}

impl<'a> Iterator for PathIter<'a> {
    type Item = (usize, Orientation);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|gbwt_node| support::decode_node(gbwt_node))
    }
}

impl<'a> FusedIterator for PathIter<'a> {}

//-----------------------------------------------------------------------------

/// An iterator over a path as a concatenation of oriented segments.
///
/// The type of `Item` is `(`[`Segment`]`, `[`Orientation`]`)`.
/// If the path is not a valid concatenation of oriented segments, the iterator stops early.
///
/// # Examples
///
/// ```
/// use gbwt::{GBZ, Orientation};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// // Path 2 in reverse orientation.
/// if let Some(mut iter) = gbz.segment_path(2, Orientation::Reverse) {
///     assert_eq!(iter.next().map(|(s, o)| (s.name, o)),
///         Some(("s17".as_bytes(), Orientation::Reverse)));
///     assert_eq!(iter.next().map(|(s, o)| (s.name, o)),
///         Some(("s16".as_bytes(), Orientation::Reverse)));
///     assert_eq!(iter.next().map(|(s, o)| (s.name, o)),
///         Some(("s14".as_bytes(), Orientation::Reverse)));
///     assert_eq!(iter.next().map(|(s, o)| (s.name, o)),
///         Some(("s13".as_bytes(), Orientation::Reverse)));
///     assert_eq!(iter.next().map(|(s, o)| (s.name, o)),
///         Some(("s11".as_bytes(), Orientation::Reverse)));
///     assert!(iter.next().is_none());
/// } else {
///     panic!("No segment iterator for path 2 (reverse)");
/// }
/// ```
#[derive(Clone, Debug)]
pub struct SegmentPathIter<'a> {
    parent: &'a GBZ,
    iter: SequenceIter<'a>,
    // Next node visit in the current segment.
    next: Option<(usize, Orientation)>,
    // Node range for the last visited segment.
    nodes: Range<usize>,
    // The path is not a valid concatenation of segments.
    fail: bool,
}

impl<'a> SegmentPathIter<'a> {
    // Visit a new segment.
    fn visit(&mut self, nodes: Range<usize>, orientation: Orientation) {
        match orientation {
            Orientation::Forward => self.next = Some((nodes.start, orientation)),
            Orientation::Reverse => self.next = Some((nodes.end - 1, orientation)),
        }
        self.nodes = nodes;
        self.advance();
    }

    // Advance in the current segment and update `self.next`.
    fn advance(&mut self) {
        let (node_id, orientation) = self.next.unwrap();
        match orientation {
            Orientation::Forward => {
                if node_id + 1 < self.nodes.end {
                    self.next = Some((node_id + 1, orientation));
                } else {
                    self.next = None;
                }
            },
            Orientation::Reverse => {
                if node_id > self.nodes.start {
                    self.next = Some((node_id - 1, orientation));
                } else {
                    self.next = None;
                }
            },
        }
    }
}

impl<'a> Iterator for SegmentPathIter<'a> {
    type Item = (Segment<'a>, Orientation);

    fn next(&mut self) -> Option<Self::Item> {
        if self.fail {
            return None;
        }
        while let Some(gbwt_node) = self.iter.next() {
            let (node_id, orientation) = support::decode_node(gbwt_node);
            if let Some(expected) = self.next {
                if (node_id, orientation) != expected {
                    self.fail = true;
                    return None;
                }
                self.advance();
            } else {
                if let Some(segment) = self.parent.node_to_segment(node_id) {
                    self.visit(segment.nodes.clone(), orientation);
                    return Some((segment, orientation));
                } else {
                    self.fail = true;
                    return None;
                }
            }
        }
        None
    }
}

impl<'a> FusedIterator for SegmentPathIter<'a> {}

//-----------------------------------------------------------------------------
