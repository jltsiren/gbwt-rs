//! GBWT: A run-length encoded FM-index storing paths as sequences of node identifiers.
//!
//! The GBWT was originally described in:
//!
//! > Sirén, Garrison, Novak, Paten, Durbin: **Haplotype-aware graph indexes**.  
//! > Bioinformatics, 2020.
//! > DOI: [10.1093/bioinformatics/btz575](https://doi.org/10.1093/bioinformatics/btz575)
//!
//! At the moment, this implementation only supports GBWT indexes built with other tools.
//! See also the original [C++ implementation](https://github.com/jltsiren/gbwt).

use crate::{ENDMARKER, SOURCE_KEY, SOURCE_VALUE};
use crate::{Orientation, Pos};
use crate::bwt::{BWT, Record};
use crate::headers::{Header, GBWTPayload, MetadataPayload};
use crate::support::{Dictionary, StringIter, Tags};
use crate::support;

use simple_sds::serialize::{Serialize, Serializable};
use simple_sds::serialize;

use std::io::{Error, ErrorKind};
use std::iter::FusedIterator;
use std::ops::Range;
use std::{io, slice};

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
/// use gbwt::{GBWT, SearchState, Orientation};
/// use gbwt::support;
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
    metadata: Option<Metadata>,
}

/// Index statistics.
impl GBWT {
    /// Returns the total length of the sequences in the index.
    #[inline]
    pub fn len(&self) -> usize {
        self.header.payload().size
    }

    /// Returns `true` if the index is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of sequences in the index.
    #[inline]
    pub fn sequences(&self) -> usize {
        self.header.payload().sequences
    }

    /// Returns the size of the alphabet.
    #[inline]
    pub fn alphabet_size(&self) -> usize {
        self.header.payload().alphabet_size
    }

    /// Returns the alphabet offset for the effective alphabet.
    #[inline]
    pub fn alphabet_offset(&self) -> usize {
        self.header.payload().offset
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
    pub fn sequence(&self, id: usize) -> Option<SequenceIter> {
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
        if let Some(record) = self.bwt.record(self.node_to_record(state.node)) {
            if let Some(range) = record.follow(state.range.clone(), node) {
                return Some(SearchState {
                    node, range,
                })
            }
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
        serialize::absent_option(writer)?; // Document array samples.
        self.metadata.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let header = Header::<GBWTPayload>::load(reader)?;
        if let Err(msg) = header.validate() {
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }

        let mut tags = Tags::load(reader)?;
        tags.insert(SOURCE_KEY, SOURCE_VALUE);

        let bwt = BWT::load(reader)?;

        // Decompress the endmarker, as the record can be poorly compressible.
        let endmarker = if bwt.is_empty() { Vec::new() } else { bwt.record(ENDMARKER).unwrap().decompress() };

        serialize::skip_option(reader)?; // Document array samples.

        // Metadata.
        let metadata = Option::<Metadata>::load(reader)?;
        if header.is_set(GBWTPayload::FLAG_METADATA) != metadata.is_some() {
            return Err(Error::new(ErrorKind::InvalidData, "GBWT: Invalid metadata flag in the header"));
        }
        if let Some(meta) = metadata.as_ref() {
            if meta.has_path_names() {
                let expected = if header.is_set(GBWTPayload::FLAG_BIDIRECTIONAL) { header.payload().sequences / 2 } else { header.payload().sequences };
                if meta.paths() > 0 && meta.paths() != expected {
                    return Err(Error::new(ErrorKind::InvalidData, "GBWT: Invalid path count in the metadata"));
                }
            }
        }

        Ok(GBWT {
            header, tags, bwt, endmarker, metadata,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.header.size_in_elements() + self.tags.size_in_elements() + self.bwt.size_in_elements() + serialize::absent_option_size() + self.metadata.size_in_elements()
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
/// The type of `Item` is [`usize`].
///
/// # Examples
///
/// ```
/// use gbwt::{GBWT, Orientation};
/// use gbwt::support;
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

/// Metadata for the paths in a GBWT index.
///
/// The metadata contains some basic statistics about the paths, and it may also contain names for paths, samples, and contigs.
/// Samples correspond to biological samples.
/// Contigs usually correspond to contigs in the reference sequence as well as to connected components in the graph.
/// The number of haplotypes is an estimate for the number of full-length paths in a connected component.
///
/// Path names are structured names associated with each path in the original graph.
/// If the GBWT index is bidirectional, it contains two sequences for each path name.
/// Each path name is a unique combination of four fields: sample, contig, phase, and fragment (see [`PathName`]).
/// The first two must be in the intervals `0..self.samples()` and `0..self.contigs()`, respectively.
///
/// # Examples
///
/// ```
/// use gbwt::{Metadata};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.meta");
/// let metadata: Metadata = serialize::load_from(&filename).unwrap();
///
/// assert!(metadata.has_path_names());
/// assert_eq!(metadata.paths(), 6);
/// assert_eq!(metadata.pan_sn_path(3), Some("sample#2#A".to_string()));
/// let path = metadata.path(3).unwrap();
///
/// assert!(metadata.has_sample_names());
/// assert_eq!(metadata.sample(path.sample()), Some("sample"));
///
/// assert!(metadata.has_contig_names());
/// assert_eq!(metadata.contig(path.contig()), Some("A"));
///
/// // Find the paths over contig B.
/// let id = metadata.contig_id("B").unwrap();
/// let mut paths = Vec::<usize>::new();
/// for (i, path) in metadata.path_iter().enumerate() {
///     if path.contig() == id {
///         paths.push(i);
///     }
/// }
/// assert_eq!(paths, vec![1, 4, 5]);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Metadata {
    header: Header<MetadataPayload>,
    path_names: Vec<PathName>,
    sample_names: Dictionary,
    contig_names: Dictionary,
}

/// Paths.
impl Metadata {
    /// Returns `true` if the metadata contains path names.
    pub fn has_path_names(&self) -> bool {
        self.header.is_set(MetadataPayload::FLAG_PATH_NAMES)
    }

    /// Returns the number path names in the metadata.
    ///
    /// If there are path names, each name corresponds to a path in the original graph.
    /// In a bidirectional GBWT, there are twice as many sequences as paths.
    pub fn paths(&self) -> usize {
        self.path_names.len()
    }

    /// Returns the name of the path with the given identifier, or [`None`] if there is no such name.
    ///
    /// Valid identifiers are in the interval `0..self.paths()`.
    pub fn path(&self, id: usize) -> Option<PathName> {
        if id < self.paths() {
            Some(self.path_names[id])
        } else {
            None
        }
    }

    /// Returns the name of the path with the given identifier in the [PanSN format](https://github.com/pangenome/PanSN-spec), or [`None`] if there is no such name.
    ///
    /// Valid identifiers are in the interval `0..self.paths()`.
    /// `'#'` will be used as the separator character.
    pub fn pan_sn_path(&self, id: usize) -> Option<String> {
        let path_name = self.path(id)?;
        let result = format!("{}#{}#{}", self.sample_name(path_name.sample()), path_name.phase(), self.contig_name(path_name.contig()));
        Some(result)
    }

    /// Returns an iterator over path names.
    pub fn path_iter(&self) -> slice::Iter<PathName> {
        self.path_names.iter()
    }
}

/// Samples.
impl Metadata {
    /// Returns `true` if the metadata contains sample names.
    pub fn has_sample_names(&self) -> bool {
        self.header.is_set(MetadataPayload::FLAG_SAMPLE_NAMES)
    }

    /// Returns the number samples.
    pub fn samples(&self) -> usize {
        self.header.payload().sample_count
    }

    /// Returns the number of haplotypes.
    ///
    /// This generally corresponds to the number of full-length paths in a graph component.
    pub fn haplotypes(&self) -> usize {
        self.header.payload().haplotype_count
    }

    /// Returns the name of the sample with the given identifier, or [`None`] if there is no such sample or name.
    ///
    /// Valid identifiers are in the interval `0..self.samples()`.
    /// Also returns [`None`] if the name exists but is not valid UTF-8.
    pub fn sample(&self, id: usize) -> Option<&str> {
        if self.has_sample_names() && id < self.samples() {
            self.sample_names.str(id).ok()
        } else {
            None
        }
    }

    /// Returns the name of the sample with the given identifier.
    ///
    /// Returns a string representation of the the sample identifier when [`Metadata::sample`] would return [`None`].
    pub fn sample_name(&self, id: usize) -> String {
        if let Some(name) = self.sample(id) {
            name.to_string()
        } else {
            id.to_string()
        }
    }
 
    /// Returns the sample idenfier corresponding to the given sample name, or [`None`] if there is no such sample.
    pub fn sample_id(&self, name: &str) -> Option<usize> {
        self.sample_names.id(name)
    }

    /// Returns an iterator over sample names.
    pub fn sample_iter(&self) -> StringIter {
        self.sample_names.as_ref().iter()
    }
}

/// Contigs.
impl Metadata {
    /// Returns `true` if the metadata contains contig names.
    pub fn has_contig_names(&self) -> bool {
        self.header.is_set(MetadataPayload::FLAG_CONTIG_NAMES)
    }

    /// Returns the number contigs.
    ///
    /// A contig usually corresponds to a graph component.
    pub fn contigs(&self) -> usize {
        self.header.payload().contig_count
    }

    /// Returns the name of the contig with the given identifier, or [`None`] if there is no such contig or name.
    ///
    /// Valid identifiers are in the interval `0..self.contigs()`.
    /// Also returns [`None`] if the name exists but is not valid UTF-8.
    pub fn contig(&self, id: usize) -> Option<&str> {
        if self.has_contig_names() && id < self.contigs() {
            self.contig_names.str(id).ok()
        } else {
            None
        }
    }

    /// Returns the name of the contig with the given identifier.
    ///
    /// Returns a string representation of the the contig identifier when [`Metadata::contig`] would return [`None`].
    pub fn contig_name(&self, id: usize) -> String {
        if let Some(name) = self.contig(id) {
            name.to_string()
        } else {
            id.to_string()
        }
    }

    /// Returns the contig idenfier corresponding to the given contig name, or [`None`] if there is no such contig.
    pub fn contig_id(&self, name: &str) -> Option<usize> {
        self.contig_names.id(name)
    }

    /// Returns an iterator over contig names.
    pub fn contig_iter(&self) -> StringIter {
        self.contig_names.as_ref().iter()
    }
}

impl Serialize for Metadata {
    fn serialize_header<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.header.serialize(writer)
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.path_names.serialize(writer)?;
        self.sample_names.serialize(writer)?;
        self.contig_names.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let header = Header::<MetadataPayload>::load(reader)?;
        if let Err(msg) = header.validate() {
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }

        let path_names = Vec::<PathName>::load(reader)?;
        if header.is_set(MetadataPayload::FLAG_PATH_NAMES) == path_names.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Metadata: Path name flag does not match the presence of path names"));
        }

        let sample_names = Dictionary::load(reader)?;
        if header.is_set(MetadataPayload::FLAG_SAMPLE_NAMES) {
            if header.payload().sample_count != sample_names.len() {
                return Err(Error::new(ErrorKind::InvalidData, "Metadata: Sample count does not match the number of sample names"));
            }
        } else if !sample_names.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Metadata: Sample names are present without the sample name flag"));
        }

        let contig_names = Dictionary::load(reader)?;
        if header.is_set(MetadataPayload::FLAG_CONTIG_NAMES) {
            if header.payload().contig_count != contig_names.len() {
                return Err(Error::new(ErrorKind::InvalidData, "Metadata: Contig count does not match the number of contig names"));
            }
        } else if !contig_names.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Metadata: Contig names are present without the contig name flag"));
        }

        Ok(Metadata {
            header, path_names, sample_names, contig_names,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.header.size_in_elements() + self.path_names.size_in_elements() + self.sample_names.size_in_elements() + self.contig_names.size_in_elements()
    }
}

//-----------------------------------------------------------------------------

/// A structured path name.
///
/// Each path name in the same metadata structure must be unique.
/// A path name has four components: sample, contig, phase (haplotype), and fragment (count).
///
/// * Samples correspond to sample identifiers in the metadata.
/// * Contigs correspond to contig identifiers in the metadata.
/// * Each (sample, phase) combination corresponds to a haplotype in the metadata.
/// * Fragment field can be used as a fragment index for fragments from the sequence identified by (sample, contig, phase).
///   It can also be used as the starting offset of the fragment in the corresponding sequence.
#[repr(C)]
#[derive(Copy, Clone, Default, Debug, Hash, PartialEq, Eq)]
pub struct PathName {
    /// Sample identifier.
    pub sample: u32,

    /// Contig identifier.
    pub contig: u32,

    /// Phase / haplotype identifier.
    pub phase: u32,

    /// Fragment identifier / running count.
    pub fragment: u32,
}

impl PathName {
    /// Returns a new path name with all components set to 0.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns a path name with the given values in each field.
    pub fn from_fields(sample: usize, contig: usize, phase: usize, fragment: usize) -> Self {
        PathName {
            sample: sample as u32,
            contig: contig as u32,
            phase: phase as u32,
            fragment: fragment as u32,
        }
    }

    /// Returns the sample identifier.
    pub fn sample(&self) -> usize {
        self.sample as usize
    }

    /// Returns the contig identifier.
    pub fn contig(&self) -> usize {
        self.contig as usize
    }

    /// Returns the phase / haplotype identifier.
    pub fn phase(&self) -> usize {
        self.phase as usize
    }

    /// Returns the fragment identifier / running count.
    pub fn fragment(&self) -> usize {
        self.fragment as usize
    }
}

impl Serializable for PathName {}

//-----------------------------------------------------------------------------
