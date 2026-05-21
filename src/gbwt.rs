//! GBWT: A run-length encoded FM-index storing paths as sequences of node identifiers.
//!
//! The GBWT was originally described in:
//!
//! > Sirén, Garrison, Novak, Paten, Durbin: **Haplotype-aware graph indexes**.  
//! > Bioinformatics, 2020.
//! > DOI: [10.1093/bioinformatics/btz575](https://doi.org/10.1093/bioinformatics/btz575)
//!
//! This implementation includes basic GBWT construction using [`crate::GBWTBuilder`].
//! For other construction options, see the original [C++ implementation](https://github.com/jltsiren/gbwt).

use crate::{ENDMARKER, SOURCE_KEY, SOURCE_VALUE};
use crate::{Orientation, Pos, Metadata};
use crate::bwt::{BWT, Record};
use crate::headers::{Header, GBWTPayload};
use crate::support::Tags;
use crate::support;

use simple_sds::serialize::Serialize;

use std::io::{Error, ErrorKind};
use std::iter::FusedIterator;
use std::ops::Range;
use std::io;

pub mod builder;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// The GBWT index storing a collection of paths space-efficiently.
///
/// The GBWT stores integer sequences.
/// Each integer is assumed to be a node identifier, and each sequence is interpreted as a path in a graph.
/// If the index is not bidirectional, GBWT node and sequence identifiers correspond directly to node and path identifiers in the original graph.
///
/// In a bidirectional index, each node (path) in the original graph becomes two nodes (sequences) in the GBWT: one for the forward orientation and one for the reverse orientation.
/// A reverse path visits the other orientation of each node on the path in reverse order.
/// The following functions can be used for mapping between the identifiers used by the GBWT and the graph:
///
/// * [`support::encode_node`], [`support::node_id`], [`support::node_orientation`], [`support::decode_node`], and [`support::flip_node`] for node identifiers.
/// * [`support::encode_path`], [`support::path_id`], [`support::path_orientation`], [`support::decode_path`], and [`support::flip_path`] for sequence / path identifiers.
///
/// # Examples
///
/// ```
/// use gbz::{GBWT, SearchState, Orientation};
/// use gbz::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbwt");
/// let index: GBWT = serialize::load_from(&filename).unwrap();
///
/// // Statistics.
/// assert_eq!(index.len(), 68);
/// assert_eq!(index.sequences(), 12);
/// assert_eq!(index.alphabet_size(), 52);
/// assert!(index.is_bidirectional());
///
/// // Manually find the second-to-last node of path 2 in forward orientation.
/// let mut pos = index.start(support::encode_path(2, Orientation::Forward));
/// let mut last = None;
/// while pos.is_some() {
///     last = pos;
///     pos = index.forward(pos.unwrap());
/// }
/// let pos = index.backward(last.unwrap()).unwrap();
/// assert_eq!(pos.node, support::encode_node(15, Orientation::Forward));
///
/// // Search for subpath (12, forward), (14, forward), (15, forward).
/// let state = index.find(support::encode_node(12, Orientation::Forward)).unwrap();
/// let state = index.extend(&state, support::encode_node(14, Orientation::Forward)).unwrap();
/// let state = index.extend(&state, support::encode_node(15, Orientation::Forward)).unwrap();
/// assert_eq!(state.node, support::encode_node(15, Orientation::Forward));
/// assert_eq!(state.len(), 2);
///
/// // Bidirectional search for the same subpath.
/// let state = index.bd_find(support::encode_node(14, Orientation::Forward)).unwrap();
/// let state = index.extend_backward(&state, support::encode_node(12, Orientation::Forward)).unwrap();
/// let state = index.extend_forward(&state, support::encode_node(15, Orientation::Forward)).unwrap();
/// assert_eq!(state.forward.node, support::encode_node(15, Orientation::Forward));
/// assert_eq!(state.reverse.node, support::encode_node(12, Orientation::Reverse));
/// assert_eq!(state.len(), 2);
///
/// // Metadata and tags.
/// assert!(index.has_metadata());
/// let metadata = index.metadata().unwrap();
/// assert_eq!(metadata.paths(), 6);
/// assert_eq!(metadata.samples(), 2);
/// assert_eq!(metadata.contigs(), 2);
/// let tags = index.tags();
/// assert!(tags.contains_key("source"));
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GBWT {
    header: Header<GBWTPayload>,
    tags: Tags,
    bwt: BWT,
    endmarker: Vec<Pos>,
    da_samples: Vec<u64>, // We pass the data through but cannot interpret it.
    metadata: Option<Metadata>,
}

/// Index statistics.
impl GBWT {
    /// Returns the total length of the sequences in the index.
    #[inline]
    pub fn len(&self) -> usize {
        self.header.payload().size as usize
    }

    /// Returns `true` if the index is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of sequences in the index.
    #[inline]
    pub fn sequences(&self) -> usize {
        self.header.payload().sequences as usize
    }

    /// Returns the size of the alphabet.
    #[inline]
    pub fn alphabet_size(&self) -> usize {
        self.header.payload().alphabet_size as usize
    }

    /// Returns the alphabet offset for the effective alphabet.
    #[inline]
    pub fn alphabet_offset(&self) -> usize {
        self.header.payload().offset as usize
    }

    /// Returns the size of the effective alphabet.
    #[inline]
    pub fn effective_size(&self) -> usize {
        self.alphabet_size() - self.alphabet_offset()
    }

    /// Returns the smallest node identifier in the effective alphabet.
    #[inline]
    pub fn first_node(&self) -> usize {
        self.alphabet_offset() + 1
    }

    /// Converts node id to record id.
    ///
    /// The record id is valid if `self.has_node(node_id)`.
    #[inline]
    pub fn node_to_record(&self, node_id: usize) -> usize {
        node_id - self.alphabet_offset()
    }

    /// Converts record id to node id.
    ///
    /// The node id is valid if `record_id != ENDMARKER` and `record_id < self.effective_size()`.
    #[inline]
    pub fn record_to_node(&self, record_id: usize) -> usize {
        record_id + self.alphabet_offset()
    }

    /// Returns `true` if node identifier `id` is in the effective alphabet.
    #[inline]
    pub fn has_node(&self, id: usize) -> bool {
        id > self.alphabet_offset() && id < self.alphabet_size()
    }

    /// Returns `true` if the GBWT index is bidirectional.
    #[inline]
    pub fn is_bidirectional(&self) -> bool {
        self.header.is_set(GBWTPayload::FLAG_BIDIRECTIONAL)
    }
}

/// Metadata and tags.
impl GBWT {
    /// Returns `true` if the index contains metadata.
    pub fn has_metadata(&self) -> bool {
        self.header.is_set(GBWTPayload::FLAG_METADATA)
    }

    /// Returns a reference to the metadata, or [`None`] if there is no metadata.
    pub fn metadata(&self) -> Option<&Metadata> {
        self.metadata.as_ref()
    }

    /// Returns a reference to the tags.
    pub fn tags(&self) -> &Tags {
        &self.tags
    }

    /// Returns a mutable reference to the tags.
    pub fn tags_mut(&mut self) -> &mut Tags {
        &mut self.tags
    }
}

impl AsRef<BWT> for GBWT {
    fn as_ref(&self) -> &BWT {
        &self.bwt
    }
}

//-----------------------------------------------------------------------------

/// Sequence navigation.
impl GBWT {
    /// Returns the first position in sequence `id`.
    ///
    /// The return value is [`None`] if no such sequence exists or the sequence is empty.
    pub fn start(&self, id: usize) -> Option<Pos> {
        if id < self.endmarker.len() && self.endmarker[id].node != ENDMARKER {
            Some(self.endmarker[id])
        } else {
            None
        }
    }

    /// Follows the sequence forward and returns the next position, or [`None`] if no such position exists.
    pub fn forward(&self, pos: Pos) -> Option<Pos> {
        // This also catches the endmarker.
        if pos.node < self.first_node() {
            return None;
        }
        let record = self.bwt.record(self.node_to_record(pos.node))?;
        record.lf(pos.offset)
    }

    /// Follows the sequence backward and returns the previous position, or [`None`] if no such position exists.
    ///
    /// # Panics
    ///
    /// Panics if the index is not bidirectional.
    pub fn backward(&self, pos: Pos) -> Option<Pos> {
        assert!(self.is_bidirectional(), "Following sequences backward requires a bidirectional GBWT");
        // This also catches the endmarker.
        if pos.node <= self.first_node() {
            return None;
        }

        let reverse_id = self.node_to_record(support::flip_node(pos.node));
        let record = self.bwt.record(reverse_id)?;
        let predecessor = record.predecessor_at(pos.offset)?;
        let pred_record = self.bwt.record(self.node_to_record(predecessor))?;
        let offset = pred_record.offset_to(pos)?;

        Some(Pos::new(predecessor, offset))
    }

    /// Returns an iterator over sequence `id`, or [`None`] if there is no such sequence.
    pub fn sequence(&'_ self, id: usize) -> Option<SequenceIter<'_>> {
        if id >= self.sequences() {
            return None;
        }
        Some(SequenceIter {
            parent: self,
            next: self.start(id),
        })
    }
}

//-----------------------------------------------------------------------------

/// Subpath search.
impl GBWT {
    /// Returns a search state for all occurrences of the given node, or [`None`] if no such node exists.
    pub fn find(&self, node: usize) -> Option<SearchState> {
        // This also catches the endmarker.
        if node < self.first_node() {
            return None;
        }
        if let Some(record) = self.bwt.record(self.node_to_record(node)) {
            return Some(SearchState {
                node,
                range: 0..record.len(),
            });
        }
        None
    }

    /// Extends the search by the given node forward and returns the new search state, or [`None`] if no such extensions exist.
    ///
    /// Assume that the current search state corresponds to a set of substring occurrences ending with the same node.
    /// This method takes all of those substrings that continue with the given node, extends them with that node, and returns the new search state.
    ///
    /// # Arguments
    ///
    /// * `state`: A search state corresponding to a set of substring occurrences.
    /// * `node`: Node to extend the substrings with.
    pub fn extend(&self, state: &SearchState, node: usize) -> Option<SearchState> {
        // This also catches the endmarker.
        if node < self.first_node() {
            return None;
        }
        if let Some(record) = self.bwt.record(self.node_to_record(state.node))
            && let Some(range) = record.follow(state.range.clone(), node) {
            return Some(SearchState {
                node, range,
            })
        }
        None
    }

    /// Returns a bidirectional search state for all occurrences of the given node, or [`None`] if no such node exists.
    ///
    /// # Panics
    ///
    /// Will panic if the index is not bidirectional.
    pub fn bd_find(&self, node: usize) -> Option<BidirectionalState> {
        assert!(self.is_bidirectional(), "Bidirectional search requires a bidirectional GBWT");
        if let Some(state) = self.find(node) {
            let reverse = SearchState {
                node: support::flip_node(state.node),
                range: state.range.clone(),
            };
            return Some(BidirectionalState {
                forward: state,
                reverse,
            });
        }
        None
    }

    /// Extends the search by the given node forward and returns the new search state, or [`None`] if no such extensions exist.
    ///
    /// Assume that the current search state corresponds to a set of substring occurrences ending with the same node.
    /// This method takes all of those substrings that continue with the given node, extends them with that node, and returns the new search state.
    ///
    /// # Arguments
    ///
    /// * `state`: A bidirectional search state corresponding to a set of substring occurrences.
    /// * `node`: Node to extend the substrings with.
    ///
    /// # Panics
    ///
    /// Will panic if the index is not bidirectional.
    pub fn extend_forward(&self, state: &BidirectionalState, node: usize) -> Option<BidirectionalState> {
        assert!(self.is_bidirectional(), "Bidirectional search requires a bidirectional GBWT");
        // This also catches the endmarker.
        if node < self.first_node() {
            return None;
        }
        let record = self.bwt.record(self.node_to_record(state.forward.node))?;
        Self::bd_internal(&record, state, node)
    }

    /// Extends the search by the given node backward and returns the new search state, or [`None`] if no such extensions exist.
    ///
    /// Assume that the current search state corresponds to a set of substring occurrences starting with the same node.
    /// This method takes all of those substrings that are preceded by the given node, extends them with that node, and returns the new search state.
    ///
    /// # Arguments
    ///
    /// * `state`: A bidirectional search state corresponding to a set of substring occurrences.
    /// * `node`: Node to extend the substrings with.
    ///
    /// # Panics
    ///
    /// Will panic if the index is not bidirectional.
    pub fn extend_backward(&self, state: &BidirectionalState, node: usize) -> Option<BidirectionalState> {
        if let Some(result) = self.extend_forward(&state.flip(), support::flip_node(node)) {
            return Some(result.flip());
        }
        None
    }

    // Internal implementation of bidirectional search. Extends the state forward.
    #[doc(hidden)]
    pub fn bd_internal(record: &Record, state: &BidirectionalState, node: usize) -> Option<BidirectionalState> {
        let (range, offset) = record.bd_follow(state.forward.range.clone(), node)?;
        let forward = SearchState {
            node, range,
        };
        let pos = state.reverse.range.start + offset;
        let reverse = SearchState {
            node: state.reverse.node,
            range: pos..pos + forward.len(),
        };
        Some(BidirectionalState {
            forward, reverse,
        })
    }
}

//-----------------------------------------------------------------------------

impl Serialize for GBWT {
    fn serialize_header<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.header.serialize(writer)
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.tags.serialize(writer)?;
        self.bwt.serialize(writer)?;
        self.da_samples.serialize(writer)?; // Document array samples.
        self.metadata.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let mut header = Header::<GBWTPayload>::load(reader)?;
        if let Err(msg) = header.validate() {
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }

        let mut tags = Tags::load(reader)?;
        tags.insert(SOURCE_KEY, SOURCE_VALUE);

        let bwt = BWT::load(reader)?;

        // Decompress the endmarker, as the record can be poorly compressible.
        let endmarker = if bwt.is_empty() { Vec::new() } else { bwt.record(ENDMARKER).unwrap().decompress() };

        // We cannot interpret the document array samples from the C++ implementation.
        let da_samples = Vec::<u64>::load(reader)?;

        // Metadata.
        let metadata = Option::<Metadata>::load(reader)?;
        if header.is_set(GBWTPayload::FLAG_METADATA) != metadata.is_some() {
            return Err(Error::new(ErrorKind::InvalidData, "GBWT: Invalid metadata flag in the header"));
        }
        if let Some(meta) = metadata.as_ref() && meta.has_path_names() {
            let expected = if header.is_set(GBWTPayload::FLAG_BIDIRECTIONAL) { (header.payload().sequences / 2) as usize } else { header.payload().sequences as usize };
            if meta.paths() > 0 && meta.paths() != expected {
                return Err(Error::new(ErrorKind::InvalidData, "GBWT: Invalid path count in the metadata"));
            }
        }

        // Update the header to the latest version after we have used the
        // serialized version for loading the correct data.
        header.update();

        Ok(GBWT {
            header, tags, bwt, endmarker, da_samples, metadata,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.header.size_in_elements() + self.tags.size_in_elements() + self.bwt.size_in_elements() + self.da_samples.size_in_elements() + self.metadata.size_in_elements()
    }
}

//-----------------------------------------------------------------------------

/// A state of unidirectional search in [`GBWT`].
///
/// The state consists of the last matched GBWT node identifier and an offset range in that node.
/// This information is equivalent to a BWT range in a normal FM-index.
///
/// Because `SearchState` contains a [`Range`], it does not implement [`Copy`].
/// As search states are often reused, they are passed by reference instead of value.
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct SearchState {
    /// GBWT node identifier for the last matched node.
    pub node: usize,
    /// Offset range in the node.
    pub range: Range<usize>,
}

impl SearchState {
    /// Returns the number of matching substring occurrences (the length of the offset range).
    #[inline]
    pub fn len(&self) -> usize {
        self.range.len()
    }

    /// Returns `true` if there are no matching substring occurrences.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.range.is_empty()
    }
}

/// A state of bidirectional search in a bidirectional [`GBWT`].
///
/// The state consists of forward and reverse search states.
/// It usually corresponds to all occurrences of a substring `pattern`.
/// The forward state is then the search state for `pattern`, while the reverse state is for the reverse pattern obtained with [`support::reverse_path`].
///
/// Because `BidirectionalState` contains a [`Range`], it does not implement [`Copy`].
/// As search states are often reused, they are passed by reference instead of value.
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct BidirectionalState {
    /// GBWT search state for the forward pattern.
    pub forward: SearchState,
    /// GBWT search state for the reverse pattern.
    pub reverse: SearchState,
}

impl BidirectionalState {
    /// Returns the number of matching substring occurrences (the length of the offset range).
    #[inline]
    pub fn len(&self) -> usize {
        self.forward.len()
    }

    /// Returns `true` if there are no matching substring occurrences.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.forward.is_empty()
    }

    /// Returns a new search state with the forward and reverse states flipped.
    pub fn flip(&self) -> BidirectionalState {
        BidirectionalState {
            forward: self.reverse.clone(),
            reverse: self.forward.clone(),
        }
    }

    /// Returns the first node on the path corresponding to the search state.
    ///
    /// The return value consists of a node identifier in the original graph and the orientation of the node.
    #[inline]
    pub fn from(&self) -> (usize, Orientation) {
        support::decode_node(support::flip_node(self.reverse.node))
    }

    /// Returns the last node on the path corresponding to the search state.
    ///
    /// The return value consists of a node identifier in the original graph and the orientation of the node.
    #[inline]
    pub fn to(&self) -> (usize, Orientation) {
        support::decode_node(self.forward.node)
    }
}

//-----------------------------------------------------------------------------

/// An iterator over a sequence in [`GBWT`].
///
/// The type of `Item` is [`prim@usize`].
///
/// # Examples
///
/// ```
/// use gbz::{GBWT, Orientation};
/// use gbz::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbwt");
/// let index: GBWT = serialize::load_from(&filename).unwrap();
///
/// // Extract path 3 in reverse orientation.
/// let path: Vec<usize> = index.sequence(support::encode_path(3, Orientation::Reverse)).unwrap().collect();
/// assert_eq!(path, vec![35, 33, 29, 27, 23]);
/// ```
#[derive(Clone, Debug)]
pub struct SequenceIter<'a> {
    parent: &'a GBWT,
    // The next position.
    next: Option<Pos>,
}

impl<'a> Iterator for SequenceIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(pos) = self.next {
            self.next = self.parent.forward(pos);
            Some(pos.node)
        } else {
            None
        }
    }
}

impl<'a> FusedIterator for SequenceIter<'a> {}

//-----------------------------------------------------------------------------

