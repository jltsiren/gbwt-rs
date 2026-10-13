//! A helper structure for assembling the sequences and node-to-segment
//! translation needed for building a [`crate::Sequences`] structure.
//!
//! The [`GraphBuilder`] plays the same role as the `NaiveGraph` structure in the
//! C++ [gbwtgraph](https://github.com/jltsiren/gbwtgraph) implementation.
//! Unlike `NaiveGraph`, it does not store the edges of the graph, as the GBZ
//! derives its topology from the [`crate::GBWT`] index.
//!
//! The builder has two mutually exclusive modes:
//!
//! * [`GraphBuilder::add_node`]: Each node is added with an explicit integer
//!   identifier. This is used when the node identifiers are already valid node
//!   ids and no splitting is needed.
//! * [`GraphBuilder::translate_segment`]: Each segment is added with a string
//!   name and is split into consecutive nodes of at most `max_length` bases.
//!   This is used when the segment names are not valid node ids or when the
//!   segments are too long.
//!
//! # Examples
//!
//! ```
//! use gbz::{GraphBuilder, Orientation};
//! use gbz::support;
//!
//! let mut graph = GraphBuilder::new();
//!
//! // Translate two segments into nodes of at most 3 bp.
//! let a = graph.translate_segment("A", b"GATTACA", 3).unwrap();
//! let b = graph.translate_segment("B", b"CG", 3).unwrap();
//!
//! assert!(graph.uses_translation());
//! assert_eq!(a, 1..4);
//! assert_eq!(b, 4..5);
//! assert_eq!(graph.node_sequence(1), b"GAT");
//! assert_eq!(graph.node_sequence(3), b"A");
//! assert_eq!(graph.node_sequence(4), b"CG");
//!
//! // Expand a path over the segments into encoded node identifiers.
//! let mut path = Vec::new();
//! graph.expand_path(&[("A", Orientation::Forward), ("B", Orientation::Reverse)], &mut path).unwrap();
//! let expected = vec![
//!     support::encode_node(1, Orientation::Forward),
//!     support::encode_node(2, Orientation::Forward),
//!     support::encode_node(3, Orientation::Forward),
//!     support::encode_node(4, Orientation::Reverse),
//! ];
//! assert_eq!(path, expected);
//! ```

use crate::support::{self, Orientation};

use std::collections::BTreeMap;
use std::ops::Range;

//-----------------------------------------------------------------------------

/// A helper structure for assembling node sequences and node-to-segment
/// translation for [`crate::Sequences`] construction.
///
/// See the module documentation for further details.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GraphBuilder {
    // Are we translating segments into nodes?
    translation: bool,
    // Node identifier to sequence. The sequence is stored in upper case.
    nodes: BTreeMap<usize, Vec<u8>>,
    // Segment translations in insertion order, which is also ascending node id order.
    // Each element is (segment name, half-open range of node ids).
    segment_ranges: Vec<(String, Range<usize>)>,
    // Segment name to index in `segment_ranges`.
    segment_index: BTreeMap<String, usize>,
    // Smallest and largest node identifiers. `min_node` is `usize::MAX` and
    // `max_node` is `0` when the builder is empty.
    min_node: usize,
    max_node: usize,
}

impl GraphBuilder {
    /// Creates an empty builder.
    pub fn new() -> Self {
        GraphBuilder {
            translation: false,
            nodes: BTreeMap::new(),
            segment_ranges: Vec::new(),
            segment_index: BTreeMap::new(),
            min_node: usize::MAX,
            max_node: 0,
        }
    }

    /// Adds a new node with the given integer identifier and sequence.
    ///
    /// The sequence is normalized to upper case.
    ///
    /// # Errors
    ///
    /// Returns an error if the identifier is `0` (reserved for the GBWT
    /// endmarker) or too large to be encoded as a GBWT node, if a node with the
    /// same identifier already exists, if the sequence is empty, or if the
    /// builder is in translation mode (see [`Self::translate_segment`]).
    pub fn add_node(&mut self, id: usize, sequence: &[u8]) -> Result<(), String> {
        if self.translation {
            return Err(String::from("GraphBuilder: Cannot add nodes by id after translating segments"));
        }
        if id == 0 {
            return Err(String::from("GraphBuilder: Node identifier 0 is reserved"));
        }
        if id > (usize::MAX - 1) / 2 {
            return Err(format!("GraphBuilder: Node identifier {} is too large", id));
        }
        if sequence.is_empty() {
            return Err(format!("GraphBuilder: Empty sequence for node {}", id));
        }
        if self.nodes.contains_key(&id) {
            return Err(format!("GraphBuilder: Duplicate node {}", id));
        }

        self.nodes.insert(id, sequence.to_ascii_uppercase());
        if id < self.min_node {
            self.min_node = id;
        }
        if id > self.max_node {
            self.max_node = id;
        }

        Ok(())
    }

    /// Takes a segment with the given name and sequence and translates it into
    /// one or more nodes of at most `max_length` bases each.
    ///
    /// The nodes receive consecutive identifiers starting from the identifier
    /// after the current largest node identifier (or from `1` if the builder is
    /// empty). The sequence is normalized to upper case.
    ///
    /// Returns the half-open range of node identifiers corresponding to the segment.
    ///
    /// # Arguments
    ///
    /// * `name`: Segment name.
    /// * `sequence`: Segment sequence.
    /// * `max_length`: Maximum node length in bases. If `0`, the segment is not split.
    ///
    /// # Errors
    ///
    /// Returns an error if a segment with the same name already exists, if the
    /// sequence is empty, or if the builder already contains nodes added with
    /// [`Self::add_node`].
    pub fn translate_segment(&mut self, name: &str, sequence: &[u8], max_length: usize) -> Result<Range<usize>, String> {
        if !self.translation && !self.nodes.is_empty() {
            return Err(String::from("GraphBuilder: Cannot translate segments after adding nodes by id"));
        }
        if sequence.is_empty() {
            return Err(format!("GraphBuilder: Empty sequence for segment {}", name));
        }
        if self.segment_index.contains_key(name) {
            return Err(format!("GraphBuilder: Duplicate segment {}", name));
        }

        let chunk_size = if max_length == 0 { sequence.len() } else { max_length };
        let start = self.max_node + 1;
        let mut offset = 0;
        let mut id = start;
        while offset < sequence.len() {
            let length = std::cmp::min(chunk_size, sequence.len() - offset);
            self.nodes.insert(id, sequence[offset..offset + length].to_ascii_uppercase());
            offset += length;
            id += 1;
        }

        let range = start..id;
        self.translation = true;
        self.segment_index.insert(String::from(name), self.segment_ranges.len());
        self.segment_ranges.push((String::from(name), range.clone()));
        if start < self.min_node {
            self.min_node = start;
        }
        if id - 1 > self.max_node {
            self.max_node = id - 1;
        }

        Ok(range)
    }

    /// Returns the half-open range of node identifiers corresponding to the
    /// segment with the given name.
    ///
    /// Returns [`None`] if there is no such segment. In non-translation mode,
    /// the name is interpreted as an integer node identifier, and the method
    /// returns `id..id + 1` if a node with that identifier exists.
    pub fn translate(&self, name: &str) -> Option<Range<usize>> {
        if self.translation {
            self.segment_index.get(name).map(|&i| self.segment_ranges[i].1.clone())
        } else {
            match name.parse::<usize>() {
                Ok(id) if self.nodes.contains_key(&id) => Some(id..id + 1),
                _ => None,
            }
        }
    }

    /// Returns `true` if the builder uses node-to-segment translation.
    ///
    /// This is the case when nodes have been created with [`Self::translate_segment`].
    #[inline]
    pub fn uses_translation(&self) -> bool {
        self.translation
    }

    /// Returns `true` if the builder contains a node with the given identifier.
    #[inline]
    pub fn has_node(&self, id: usize) -> bool {
        self.nodes.contains_key(&id)
    }

    /// Returns the number of nodes in the builder.
    #[inline]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the smallest node identifier, or [`None`] if the builder is empty.
    #[inline]
    pub fn min_node(&self) -> Option<usize> {
        if self.nodes.is_empty() { None } else { Some(self.min_node) }
    }

    /// Returns the largest node identifier, or [`None`] if the builder is empty.
    #[inline]
    pub fn max_node(&self) -> Option<usize> {
        if self.nodes.is_empty() { None } else { Some(self.max_node) }
    }

    /// Returns the sequence of the node with the given identifier.
    ///
    /// Returns an empty slice if there is no such node.
    #[inline]
    pub fn node_sequence(&self, id: usize) -> &[u8] {
        match self.nodes.get(&id) {
            Some(sequence) => sequence.as_slice(),
            None => &[],
        }
    }

    /// Returns the number of segments in the translation.
    #[inline]
    pub fn segment_count(&self) -> usize {
        self.segment_ranges.len()
    }

    /// Returns an iterator over the segments in the translation.
    ///
    /// The type of `Item` is `(`[`&str`]`, `[`Range`]`<[`usize`]`>)`.
    /// The segments are listed in ascending order by their first node identifier.
    pub fn segment_iter(&self) -> impl Iterator<Item = (&str, Range<usize>)> {
        self.segment_ranges.iter().map(|(name, range)| (name.as_str(), range.clone()))
    }

    /// Expands a path over segments into a sequence of encoded node identifiers.
    ///
    /// Each step is a segment name and an orientation. A step in forward
    /// orientation expands to the segment nodes in ascending order, each encoded
    /// with [`support::encode_node`] in forward orientation. A step in reverse
    /// orientation expands to the segment nodes in descending order, each
    /// encoded in reverse orientation.
    ///
    /// # Arguments
    ///
    /// * `steps`: Path as a slice of (segment name, orientation) pairs.
    /// * `out`: Vector that the encoded node identifiers are appended to.
    ///
    /// # Errors
    ///
    /// Returns an error if a segment name does not exist.
    pub fn expand_path(&self, steps: &[(&str, Orientation)], out: &mut Vec<usize>) -> Result<(), String> {
        for &(name, orientation) in steps {
            let range = match self.translate(name) {
                Some(range) => range,
                None => return Err(format!("GraphBuilder: Invalid segment {}", name)),
            };
            match orientation {
                Orientation::Forward => {
                    for id in range {
                        out.push(support::encode_node(id, Orientation::Forward));
                    }
                }
                Orientation::Reverse => {
                    for id in range.rev() {
                        out.push(support::encode_node(id, Orientation::Reverse));
                    }
                }
            }
        }
        Ok(())
    }
}

//-----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_node() {
        let mut graph = GraphBuilder::new();
        assert!(!graph.uses_translation());
        assert_eq!(graph.node_count(), 0);
        assert_eq!(graph.min_node(), None);
        assert_eq!(graph.max_node(), None);

        graph.add_node(3, b"gatT").unwrap();
        graph.add_node(1, b"acg").unwrap();
        assert_eq!(graph.node_count(), 2);
        assert!(!graph.uses_translation());
        assert_eq!(graph.min_node(), Some(1));
        assert_eq!(graph.max_node(), Some(3));

        // Sequences are normalized to upper case.
        assert_eq!(graph.node_sequence(1), b"ACG");
        assert_eq!(graph.node_sequence(3), b"GATT");
        assert_eq!(graph.node_sequence(2), b"");

        // Non-translation translate().
        assert_eq!(graph.translate("3"), Some(3..4));
        assert_eq!(graph.translate("2"), None);
        assert_eq!(graph.translate("foo"), None);

        // Error cases.
        assert!(graph.add_node(3, b"A").is_err());
        assert!(graph.add_node(4, b"").is_err());
        assert!(graph.add_node(0, b"A").is_err());
        assert!(graph.add_node(usize::MAX, b"A").is_err());
        assert_eq!(graph.node_count(), 2);
    }

    #[test]
    fn add_node_then_translate_fails() {
        let mut graph = GraphBuilder::new();
        graph.add_node(1, b"A").unwrap();
        assert!(graph.translate_segment("S", b"AC", 1).is_err());
    }

    #[test]
    fn translate_segment_then_add_node_fails() {
        let mut graph = GraphBuilder::new();
        graph.translate_segment("S", b"AC", 1).unwrap();
        assert!(graph.add_node(1, b"A").is_err());
    }

    #[test]
    fn translate_segment_chunking() {
        let mut graph = GraphBuilder::new();
        graph.translate_segment("A", b"GATTACA", 3).unwrap();
        graph.translate_segment("B", b"CG", 3).unwrap();

        assert!(graph.uses_translation());
        assert_eq!(graph.node_count(), 4);
        assert_eq!(graph.min_node(), Some(1));
        assert_eq!(graph.max_node(), Some(4));
        assert_eq!(graph.segment_count(), 2);

        // A is split into 3 nodes of lengths 3, 3, 1.
        assert_eq!(graph.translate("A"), Some(1..4));
        assert_eq!(graph.node_sequence(1), b"GAT");
        assert_eq!(graph.node_sequence(2), b"TAC");
        assert_eq!(graph.node_sequence(3), b"A");
        // B is a single node.
        assert_eq!(graph.translate("B"), Some(4..5));
        assert_eq!(graph.node_sequence(4), b"CG");

        // Segments are listed in ascending node id order.
        let segments: Vec<(String, Range<usize>)> = graph.segment_iter().map(|(name, range)| (String::from(name), range)).collect();
        assert_eq!(segments, vec![(String::from("A"), 1..4), (String::from("B"), 4..5)]);

        assert_eq!(graph.translate("C"), None);
    }

    #[test]
    fn translate_segment_no_split() {
        let mut graph = GraphBuilder::new();
        // max_length 0 means no splitting.
        graph.translate_segment("A", b"GATTACA", 0).unwrap();
        assert_eq!(graph.translate("A"), Some(1..2));
        assert_eq!(graph.node_sequence(1), b"GATTACA");
    }

    #[test]
    fn translate_segment_errors() {
        let mut graph = GraphBuilder::new();
        assert!(graph.translate_segment("A", b"", 1).is_err());
        // A failed call must not switch the builder to translation mode.
        assert!(!graph.uses_translation());
        graph.add_node(1, b"A").unwrap();
        let mut graph = GraphBuilder::new();
        graph.translate_segment("A", b"AC", 1).unwrap();
        assert!(graph.translate_segment("A", b"GT", 1).is_err());
    }

    #[test]
    fn expand_path() {
        let mut graph = GraphBuilder::new();
        graph.translate_segment("A", b"GAT", 3).unwrap();
        graph.translate_segment("B", b"CG", 3).unwrap();

        let mut path = Vec::new();
        graph.expand_path(&[("A", Orientation::Forward), ("B", Orientation::Reverse)], &mut path).unwrap();
        let expected = vec![
            support::encode_node(1, Orientation::Forward),
            support::encode_node(2, Orientation::Reverse),
        ];
        assert_eq!(path, expected);

        let mut path = Vec::new();
        graph.expand_path(&[("B", Orientation::Reverse), ("A", Orientation::Forward)], &mut path).unwrap();
        let expected = vec![
            support::encode_node(2, Orientation::Reverse),
            support::encode_node(1, Orientation::Forward),
        ];
        assert_eq!(path, expected);

        assert!(graph.expand_path(&[("C", Orientation::Forward)], &mut Vec::new()).is_err());
    }

    #[test]
    fn expand_path_no_translation() {
        let mut graph = GraphBuilder::new();
        graph.add_node(2, b"GAT").unwrap();
        graph.add_node(5, b"CG").unwrap();

        let mut path = Vec::new();
        graph.expand_path(&[("2", Orientation::Forward), ("5", Orientation::Reverse)], &mut path).unwrap();
        let expected = vec![
            support::encode_node(2, Orientation::Forward),
            support::encode_node(5, Orientation::Reverse),
        ];
        assert_eq!(path, expected);
    }
}
