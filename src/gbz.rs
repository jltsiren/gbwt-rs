//! GBZ: Space-efficient representation for a subset of GFA.
// FIXME document

use crate::{ENDMARKER, SOURCE_KEY, SOURCE_VALUE};
use crate::{Graph, GBWT, Orientation};
use crate::bwt::Record;
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

// FIXME
//#[cfg(test)]
//mod tests;

//-----------------------------------------------------------------------------

// FIXME note: node identifiers in the original graph vs in the GBWT; missing nodes
// FIXME document, example, tests
/// GBZ file format for storing GFA graphs with many paths.
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
/// assert_eq!(gbz.nodes(), 12);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GBZ {
    header: Header<GBZPayload>,
    tags: Tags,
    index: GBWT,
    graph: Graph,
    real_nodes: BitVector,
}

//-----------------------------------------------------------------------------

// FIXME example, tests
/// Statistics.
impl GBZ {
    /// Returns the number of nodes in the graph.
    #[inline]
    pub fn nodes(&self) -> usize {
        self.graph.nodes()
    }

    /// Returns `true` if the graph contains a node with identifier `node_id`.
    #[inline]
    pub fn has_node(&self, node_id: usize) -> bool {
        let gbwt_node = support::encode_node(node_id, Orientation::Forward);
        self.index.has_node(gbwt_node) && self.real_nodes.get(Self::gbwt_node_to_sequence(&self.index, gbwt_node))
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

// FIXME iterator similar to follow_paths() in the C++ version?
/// Nodes and edges.
impl GBZ {
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
        if node_id == ENDMARKER {
            return None;
        }
        let gbwt_node = support::encode_node(node_id, orientation);
        let record_id = self.index.node_to_record(gbwt_node);
        if let Some(record) = self.index.as_ref().record(record_id) {
            let next = if record.outdegree() > 0 && record.successor(0) == ENDMARKER { 1 } else { 0 };
            let limit = record.outdegree();
            return Some(EdgeIter {
                record: record,
                next: next,
                limit: limit,
                flip: false,
            });
        }
        return None;
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
        if node_id == ENDMARKER {
            return None;
        }
        let gbwt_node = support::encode_node(node_id, orientation.flip());
        let record_id = self.index.node_to_record(gbwt_node);
        if let Some(record) = self.index.as_ref().record(record_id) {
            let next = if record.outdegree() > 0 && record.successor(0) == ENDMARKER { 1 } else { 0 };
            let limit = record.outdegree();
            return Some(EdgeIter {
                record: record,
                next: next,
                limit: limit,
                flip: true,
            });
        }
        return None;
    }
}

//-----------------------------------------------------------------------------

// FIXME need some kind of mapping from node id to (segment id, node id range)
/// Segments and links.
impl GBZ {
    /// Returns `true` if the graph contains a node-to-segment translation.
    #[inline]
    pub fn has_translation(&self) -> bool {
        self.graph.has_translation()
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
}

//-----------------------------------------------------------------------------

// FIXME tests
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

// FIXME tests
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

// FIXME tests
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

// FIXME tests; need to test skipping over unused segments
/// A read-only iterator over the segments in the GFA graph.
///
/// The type of `Item` is `(&[`[`u8`]`], `[`Range`]`<`[`usize`]`>, &[`[`u8`]`])`.
/// This corresponds to (segment name, node id range, sequence).
/// Unlike [`crate::graph::SegmentIter`], this iterator will skip segments corresponding to unused node ids.
///
/// # Examples
///
/// ```
/// use gbwt::GBZ;
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// assert!(gbz.has_translation());
/// let mut iter = gbz.segment_iter().unwrap();
/// assert_eq!(iter.next(), Some(("s11".as_bytes(), 1..3, "GAT".as_bytes())));
/// assert_eq!(iter.next(), Some(("s12".as_bytes(), 3..4, "T".as_bytes())));
/// ```
#[derive(Clone, Debug)]
pub struct SegmentIter<'a> {
    parent: &'a GBZ,
    iter: GraphSegmentIter<'a>,
}

impl<'a> Iterator for SegmentIter<'a> {
    type Item = (&'a [u8], Range<usize>, &'a [u8]);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((name, range, sequence)) = self.iter.next() {
            if self.parent.has_node(range.start) {
                return Some((name, range, sequence));
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
        while let Some((name, range, sequence)) = self.iter.next_back() {
            if self.parent.has_node(range.start) {
                return Some((name, range, sequence));
            }
        }
        None
    }
}

impl<'a> FusedIterator for SegmentIter<'a> {}

//-----------------------------------------------------------------------------

// FIXME LinkIter?

//-----------------------------------------------------------------------------
