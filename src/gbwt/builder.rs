//! [`GBWT`] construction.
//!
//! Like the C++ implementation, the construction algorithm uses 32-bit integers internally to save space.
//! This limits the maximum node identifier and the number of visits to a node to approximately [`u32::MAX`].
//! The interface uses [`usize`], as it is the semantically correct type and also used in the rest of the codebase.
//!
//! [`MutableGBWT`] provides a way of building [`GBWT`] incrementally by inserting batches of sequences.
//! [`GBWTBuilder`] is a more convenient wrapper for [`MutableGBWT`].
//! The user can choose whether to build a bidirectional GBWT and then insert paths one by one.
//! The builder buffers the inserted paths and does actual construction in a background thread.
//!
//! The construction algorithm is based on the BCR algorithm:
//!
//! > Markus J. Bauer, Anthony J. Cox, and Giovanna Rosone:\
//! > **Lightweight algorithms for constructing and inverting the BWT of string collections**.\
//! > Theoretical Computer Science 483:134–148, 2013.
//! > DOI: [10.1016/j.tcs.2012.02.002](https://doi.org/10.1016/j.tcs.2012.02.002)
//!
//! NOTE: The integrated tests using `cargo test` only cover small instances.
//! If the construction logic is modified, it should be tested with the `rebuild-gbwt` binary against existing GBWT indexes.

use crate::{ENDMARKER, SOURCE_KEY, SOURCE_VALUE};
use crate::{GBWT, Pos, FullPathName, Metadata, MetadataBuilder};
use crate::bwt::{BWT, BWTBuilder, SmallPos};
use crate::headers::{Header, GBWTPayload};
use crate::support::{self, Run, SmallRun, Tags, EdgeList};

use std::sync::mpsc::{SyncSender, Receiver};
use std::thread::JoinHandle;
use std::cmp;

//-----------------------------------------------------------------------------

/// A builder for constructing [`GBWT`] incrementally.
///
/// The user can insert paths one by one using [`Self::insert`].
/// Each path will be converted into one or two sequences, depending on whether the builder is bidirectional.
/// When building a GBWT with metadata, each path must be accompanied by a name.
///
/// The builder inserts the paths into buffers.
/// Once the buffers become full, they are sent to a background thread that inserts them into a [`MutableGBWT`].
/// Construction finishes when the user calls [`Self::build`], which also flushes any remaining data in the buffers.
///
/// # Examples
///
/// ```
/// use gbz::{GBWTBuilder, Orientation};
/// use gbz::support;
///
/// // Build a bidirectional GBWT without metadata, with a buffer capacity of 20 nodes.
/// // We ignore the error messages for simplicity.
/// let mut builder = GBWTBuilder::new(true, false, 20);
/// builder.insert(&[2, 6, 10], None).expect("Failed to insert path");
/// builder.insert(&[4, 8, 12], None).expect("Failed to insert path");
/// builder.insert(&[2, 8, 10], None).expect("Failed to insert path");
/// builder.insert(&[4, 6, 12], None).expect("Failed to insert path");
/// let mut gbwt = builder.build().expect("Failed to build GBWT");
///
/// assert!(gbwt.is_bidirectional());
/// assert_eq!(gbwt.len(), 32);
/// assert_eq!(gbwt.sequences(), 8);
///
/// // We get the reverse complement of the third path (id 2).
/// let iter = gbwt.sequence(support::encode_path(2, Orientation::Reverse));
/// assert!(iter.is_some());
/// let sequence: Vec<usize> = iter.unwrap().collect();
/// assert_eq!(sequence, &[11, 9, 3]);
/// ```
pub struct GBWTBuilder {
    // Are we building a bidirectional GBWT?
    bidirectional: bool,
    // `sequence_buffer` is not allowed to exceed this size, unless a single path needs more space.
    buffer_size: usize,
    // Buffer for concatenated sequences.
    seq_buffer: Vec<u32>,
    // Buffer for path names, if we are building a GBWT with path metadata.
    name_buffer: Option<Vec<FullPathName>>,
    // Construction thread that runs in the background.
    construction_thread: JoinHandle<Result<GBWT, String>>,
    // Channel for sending buffers to the construction thread.
    to_construction: SyncSender<(Vec<u32>, Option<Vec<FullPathName>>)>,
}

fn build_gbwt(bidirectional: bool, with_metadata: bool, receiver: Receiver<(Vec<u32>, Option<Vec<FullPathName>>)>)
    -> Result<GBWT, String>
{
    let mut builder = MutableGBWT::new(bidirectional, with_metadata);
    let mut error = None;
    while let Ok((seq_buffer, name_buffer)) = receiver.recv() {
        if seq_buffer.is_empty() {
            break;
        }
        if error.is_some() {
            continue;
        }
        // If we have a construction error, we store it and ignore subsequent batches
        // to avoid blocking the main thread.
        if let Err(e) = builder.insert(&seq_buffer, name_buffer.as_deref()) {
            error = Some(e);
        }
    }
    if let Some(e) = error {
        return Err(e);
    }
    Ok(GBWT::from(builder))
}

impl GBWTBuilder {
    /// Creates a new builder.
    ///
    /// # Arguments
    ///
    /// * `bidirectional`: Are we building a bidirectional GBWT?
    /// * `with_metadata`: Are we building a GBWT with path metadata?
    /// * `buffer_size`: Total length of a batch of sequences, including endmarkers, in nodes.
    ///   This will increase automatically if a single path exceeds the capacity.
    pub fn new(bidirectional: bool, with_metadata: bool, buffer_size: usize) -> Self {
        // We use up to three buffers: one in the main thread, one in the channel, and one in the construction thread.
        let (to_construction, receiver) = std::sync::mpsc::sync_channel(1);
        let construction_thread = std::thread::spawn(move || build_gbwt(bidirectional, with_metadata, receiver));
        Self {
            bidirectional,
            buffer_size,
            seq_buffer: Vec::with_capacity(buffer_size),
            name_buffer: if with_metadata { Some(Vec::new()) } else { None },
            construction_thread,
            to_construction,
        }
    }

    fn flush_buffers(&mut self) {
        if self.seq_buffer.is_empty() {
            return;
        }
        let seq_buffer = std::mem::replace(&mut self.seq_buffer, Vec::with_capacity(self.buffer_size));
        let name_buffer = if let Some(name_buffer) = self.name_buffer.as_mut() {
            Some(std::mem::take(name_buffer))
        } else {
            None
        };
        let _ = self.to_construction.send((seq_buffer, name_buffer));
    }

    /// Inserts a new path into the builder.
    ///
    /// # Errors
    ///
    /// Returns an error if the path contains an endmarker value ([`ENDMARKER`]).
    /// Returns an error if no path name was provided to a builder with metadata, or if a path name was provided to a builder without metadata.
    pub fn insert(&mut self, path: &[usize], name: Option<FullPathName>) -> Result<(), String> {
        if path.contains(&ENDMARKER) {
            return Err(String::from("GBWTBuilder: Path contains endmarker"));
        }
        if name.is_some() != self.name_buffer.is_some() {
            return Err(String::from("GBWTBuilder: Path names must be provided iff the builder is configured to have metadata"));
        }

        // Flush the buffers if we cannot fit the new path.
        let mut space_needed = path.len() + 1;
        if self.bidirectional {
            space_needed *= 2;
        }
        if self.seq_buffer.len() + space_needed > self.buffer_size {
            if space_needed > self.buffer_size {
                // Grow the buffer to fit the path.
                self.buffer_size = space_needed;
            }
            self.flush_buffers();
        }

        // Insert the path and the name into the buffers. If the buffers become
        // full, the next `insert` or `build` call will flush them. We do not
        // flush them here to avoid blocking the main thread if the construction
        // thread is busy.
        self.seq_buffer.extend(path.iter().map(|&node| node as u32));
        self.seq_buffer.push(ENDMARKER as u32);
        if self.bidirectional {
            support::reverse_path_into(path, &mut self.seq_buffer);
            self.seq_buffer.push(ENDMARKER as u32);
        }
        if let Some(name_buffer) = self.name_buffer.as_mut() {
            name_buffer.push(name.unwrap());
        }

        Ok(())
    }

    /// Finishes the construction and returns the built GBWT.
    ///
    /// # Errors
    ///
    /// Returns the first error the construction thread encountered, if any.
    pub fn build(self) -> Result<GBWT, String> {
        // Flush any remaining paths in the buffers, and send an empty buffer to signal the end of construction.
        let mut builder = self;
        builder.flush_buffers();
        let _ = builder.to_construction.send((Vec::new(), None));

        // Wait for the construction thread to finish and return the result.
        builder.construction_thread.join().unwrap_or_else(|_| Err(String::from("GBWTBuilder: Construction thread panicked")))
    }
}

//-----------------------------------------------------------------------------

/// A mutable node record corresponding to [`crate::bwt::Record`] used for GBWT construction.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MutableRecord {
    // Number of sequence visits to the node (BWT fragment length).
    len: usize,
    // Incoming edges as (source node id, number of visits from that node).
    incoming: EdgeList,
    // Outgoing edges as (target node id, BWT offset of the first visit to that node).
    outgoing: EdgeList,
    // Run-length encoded BWT fragment as a vector of (node id, run length) pairs.
    bwt: Vec<SmallRun>,
}

impl MutableRecord {
    /// Creates an empty record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the number of visits to the node (BWT fragment length).
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the record is empty (this should not happen).
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Sets the number of sequence visits from the given node.
    pub fn set_visits(&mut self, from: usize, count: usize) {
        if let Some(edge_count) = self.incoming.get_mut(from) {
            *edge_count = count as u32;
        } else {
            self.incoming.increment(from, count);
        }
    }

    /// Increases the number of sequence visits from the given node.
    pub fn add_visits(&mut self, from: usize, count: usize) {
        self.incoming.increment(from, count);
    }

    /// Adds an outgoing edge to the given node with an unknown offset.
    ///
    /// Does nothing if the edge already exists.
    pub fn add_edge(&mut self, to: usize) {
        self.outgoing.increment(to, 0);
    }

    /// Sets the offset of the outgoing edge to the given node.
    ///
    /// Does nothing if the edge does not exist.
    pub fn set_edge_offset(&mut self, to: usize, offset: usize) {
        if let Some(edge_offset) = self.outgoing.get_mut(to) {
            *edge_offset = offset as u32;
        }
    }

    /// Computes LF-mapping at the given offset.
    ///
    /// Returns the GBWT position the path at the given offset continues to, or [`None`] if the offset is out of bounds.
    pub fn follow_path(&self, offset: usize) -> Option<Pos> {
        if offset >= self.len {
            return None;
        }
        let mut bwt_offset = 0;
        let mut ranks = self.outgoing.clone();
        for run in self.bwt.iter() {
            bwt_offset += run.len as usize;
            ranks.increment(run.value as usize, run.len as usize);
            if bwt_offset > offset {
                let rank = ranks.get(run.value as usize).unwrap_or(0) as usize;
                let rank = rank - (bwt_offset - offset);
                return Some(Pos::new(run.value as usize, rank));
            }
        }
        None // This should be unreachable.
    }
}

//-----------------------------------------------------------------------------

// A sequence iterator used for GBWT construction.
// Sort order is `(next.node, pos)`, which is consistent with LF-mapping with `next.node`.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Sequence<'a> {
    // Sequence as a slice of node identifiers, excluding the endmarker.
    // We do not store the sequence offset, as it is shared between all sequences.
    nodes: &'a [u32],
    // GBWT position for the current sequence offset.
    curr: SmallPos,
    // GBWT position for the next sequence offset.
    // Offset in the position may be undetermined.
    next: SmallPos,
}

impl<'a> PartialOrd for Sequence<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> Ord for Sequence<'a> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.next.node, self.curr).cmp(&(other.next.node, other.curr))
    }
}

// A support structure for building a sequence of maximal runs incrementally.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct RunBuilder {
    source_len: usize,
    target_len: usize,
    runs: Vec<SmallRun>,
    ranks: EdgeList,
}

impl RunBuilder {
    fn new() -> Self {
        Self::default()
    }

    // Inserts a source run until the given offset, and returns the remaining run.
    fn insert_source_run(&mut self, source: SmallRun, offset: usize) -> SmallRun {
        let (to_insert, to_return) = if offset >= self.target_len + source.len as usize {
            (source, SmallRun::default())
        } else {
            let to_insert = SmallRun::new(source.value as usize, offset - self.target_len);
            let to_return = SmallRun::new(source.value as usize, self.target_len + source.len as usize - offset);
            (to_insert, to_return)
        };
        if let Some(last) = self.runs.last_mut() && last.value == to_insert.value {
            last.len += to_insert.len;
        } else {
            self.runs.push(to_insert);
        }
        self.source_len += to_insert.len as usize;
        self.target_len += to_insert.len as usize;
        self.ranks.increment(to_insert.value as usize, to_insert.len as usize);
        to_return
    }

    // Inserts a target value.
    fn insert_target_value(&mut self, value: usize) {
        if let Some(last) = self.runs.last_mut() && last.value == value as u32 {
            last.len += 1;
        } else {
            self.runs.push(SmallRun::new(value, 1));
        }
        self.target_len += 1;
        self.ranks.increment(value, 1);
    }

    // Returns the rank of the given value at the current target length.
    fn rank(&self, value: usize) -> u32 {
        self.ranks.get(value).unwrap_or(0)
    }
}

//-----------------------------------------------------------------------------

// TODO: from_gbwt (lossy conversion: drop partial metadata, DA samples)
/// A data structure for building [`GBWT`] indexes.
///
/// An empty index can be created with [`Self::new`].
/// Each [`Self::insert`] inserts a batch of sequences.
/// If the sequences are short and/or similar, inserting large batches can be much more efficient:
///
/// * Each step extends each sequence in the batch by one node.
///   If multiple sequences visit the same node in the same step, the node record is updated only once.
/// * If all sequences traverse the same part of the graph, the algorithm may benefit from memory locality.
///
/// # Examples
///
/// ```
/// use gbz::{GBWT, ENDMARKER};
/// use gbz::gbwt::builder::MutableGBWT;
///
/// // Unidirectional GBWT with no metadata.
/// let mut builder = MutableGBWT::new(false, false);
/// assert!(!builder.is_bidirectional());
/// assert!(!builder.has_metadata());
///
/// // Insert four sequences in two batches.
/// let buffer = vec![
///     1, 3, 5, ENDMARKER as u32,
///     2, 4, 6, ENDMARKER as u32,
/// ];
/// builder.insert(&buffer, None).expect("Failed to insert sequences");
/// let buffer = vec![
///     1, 4, 5, ENDMARKER as u32,
///     2, 3, 6, ENDMARKER as u32,
/// ];
/// builder.insert(&buffer, None).expect("Failed to insert sequences");
///
/// assert_eq!(builder.len(), 16);
/// assert_eq!(builder.sequences(), 4);
/// assert_eq!(builder.sequence(2), Some(vec![1, 4, 5]));
///
/// // Convert to GBWT.
/// let gbwt = GBWT::from(builder);
/// assert_eq!(gbwt.len(), 16);
/// assert_eq!(gbwt.sequences(), 4);
/// let iter = gbwt.sequence(2);
/// assert!(iter.is_some());
/// let sequence: Vec<usize> = iter.unwrap().collect();
/// assert_eq!(sequence, vec![1, 4, 5]);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MutableGBWT {
    // Total length of the sequences, including endmarkers.
    len: usize,
    // Is the GBWT bidirectional?
    bidirectional: bool,
    // GBWT tags.
    tags: Tags,
    // Endmarker record as an array of starting positions for each sequence.
    endmarker: Vec<Pos>,
    // Outgoing edges from the endmarker as (node id, number of sequences starting with that node).
    endmarker_edges: EdgeList,
    // Node identifier for record 0. This is `header.offset + 1`.
    first_node: usize,
    // Node records as a map from `node id - first_node` to mutable record.
    // We use `Option<Box<_>>` to minimize the overhead from resizing `records`.
    records: Vec<Option<Box<MutableRecord>>>,
    // Optional metadata.
    metadata: Option<MetadataBuilder>,
}

/// Queries.
impl MutableGBWT {
    /// Returns the total length of the sequences, including endmarkers.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the GBWT is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns `true` if the GBWT is bidirectional.
    pub fn is_bidirectional(&self) -> bool {
        self.bidirectional
    }

    /// Returns `true` if the GBWT has path metadata.
    pub fn has_metadata(&self) -> bool {
        self.metadata.is_some()
    }

    /// Returns the minimum node identifier, or [`None`] if there are no nodes.
    pub fn min_node(&self) -> Option<usize> {
        if self.records.is_empty() {
            None
        } else {
            Some(self.first_node)
        }
    }

    /// Returns the maximum node identifier, or [`None`] if there are no nodes.
    pub fn max_node(&self) -> Option<usize> {
        if self.records.is_empty() {
            None
        } else {
            Some(self.first_node + self.records.len() - 1)
        }
    }

    /// Returns the total number of sequences.
    pub fn sequences(&self) -> usize {
        self.endmarker.len()
    }

    /// Returns the total number of paths.
    pub fn paths(&self) -> usize {
        if self.bidirectional {
            self.endmarker.len() / 2
        } else {
            self.endmarker.len()
        }
    }

    /// Returns a reference to the GBWT tags.
    pub fn tags(&self) -> &Tags {
        &self.tags
    }

    /// Returns a mutable reference to the GBWT tags.
    pub fn tags_mut(&mut self) -> &mut Tags {
        &mut self.tags
    }

    /// Returns the given sequence, or [`None`] if the index is out of bounds.
    ///
    /// This is mostly intended for testing.
    pub fn sequence(&self, sequence_id: usize) -> Option<Vec<u32>> {
        if sequence_id >= self.endmarker.len() {
            return None;
        }

        let mut result = Vec::new();
        let mut pos = self.endmarker[sequence_id];
        while pos.node != ENDMARKER {
            result.push(pos.node as u32);
            let record = self.records.get(pos.node - self.first_node)?.as_ref()?;
            pos = record.follow_path(pos.offset)?;
        }

        Some(result)
    }

    /// Returns the name of the given path, or [`None`] if the index does not have metadata or the path id is out of bounds.
    ///
    /// This is mostly intended for testing.
    pub fn path_name(&self, path_id: usize) -> Option<FullPathName> {
        let metadata = self.metadata.as_ref()?;
        metadata.path_name(path_id)
    }
}

/// Construction.
impl MutableGBWT {
    /// Creates an empty GBWT.
    ///
    /// If `bidirectional` is `true`, the GBWT will be bidirectional.
    /// Each [`Self::insert`] call must then consist of an even number of sequences.
    /// The sequences at odd indices are assumed to be the reverse complements of the sequences at the preceding even indices.
    ///
    /// If `with_metadata` is `true`, the GBWT will have path metadata.
    /// Each [`Self::insert`] call must then be accompanied by path names for the inserted sequences, one per path.
    /// If the GBWT is bidirectional, each name corresponds to a pair of sequences.
    /// Otherwise each sequence is a separate path.
    ///
    /// The newly created index will contain the [`SOURCE_KEY`] tag but no other tags.
    pub fn new(bidirectional: bool, with_metadata: bool) -> Self {
        let mut tags = Tags::new();
        tags.insert(SOURCE_KEY, SOURCE_VALUE);
        Self {
            len: 0,
            bidirectional,
            tags,
            endmarker: Vec::new(),
            endmarker_edges: EdgeList::new(),
            first_node: 0,
            records: Vec::new(),
            metadata: if with_metadata { Some(MetadataBuilder::new()) } else { None },
        }
    }

    // Returns the sequences in the buffer without determining the GBWT positions.
    fn sequences_in_buffer<'a>(&'_ mut self, buffer: &'a [u32]) -> Result<(Vec<Sequence<'a>>, Option<(usize, usize)>), String> {
        let mut sequences: Vec<Sequence<'a>> = Vec::new();
        let mut start = 0;
        let mut node_range = None;
        for (i, &node) in buffer.iter().enumerate() {
            if node == ENDMARKER as u32 {
                sequences.push(Sequence {
                    nodes: &buffer[start..i],
                    curr: SmallPos::default(),
                    next: SmallPos::default(),
                });
                start = i + 1;
            } else {
                let node = node as usize;
                if let Some((min, max)) = node_range {
                    if node < min {
                        node_range = Some((node, max));
                    }
                    if node > max {
                        node_range = Some((min, node));
                    }
                } else {
                    node_range = Some((node, node));
                }
            }
        }
        if start < buffer.len() {
            return Err(String::from("MutableGBWT: Trailing nodes in the buffer"));
        }
        if self.bidirectional && !sequences.len().is_multiple_of(2) {
            return Err(String::from("MutableGBWT: Odd number of sequences in the buffer for a bidirectional GBWT"));
        }

        Ok((sequences, node_range))
    }

    // Appends the given path names to the metadata.
    fn append_metadata(&mut self, seq_count: usize, names: Option<&[FullPathName]>) -> Result<(), String> {
        match (&mut self.metadata, names) {
            (Some(metadata), Some(names)) => {
                let expected = if self.bidirectional { seq_count / 2 } else { seq_count };
                if names.len() != expected {
                    return Err(format!("MutableGBWT: Expected {} path names, got {}", expected, names.len()));
                }
                metadata.extend(names)?;
                Ok(())
            }
            (None, Some(_)) => Err(String::from("MutableGBWT: Path names provided for a GBWT without metadata")),
            (Some(_), None) => Err(String::from("MutableGBWT: No path names provided for a GBWT with metadata")),
            (None, None) => Ok(()),
        }
    }

    // Ensures that `self.records` covers node ids from `min_node` to `max_node`.
    fn expand_records(&mut self, min_node: usize, max_node: usize) {
        if self.records.is_empty() {
            self.first_node = min_node;
            self.records = vec![None; max_node - min_node + 1];
        } else if min_node < self.first_node || max_node >= self.first_node + self.records.len() {
            let new_first_node = cmp::min(min_node, self.first_node);
            let new_last_node = cmp::max(max_node, self.first_node + self.records.len() - 1);
            let mut new_records = Vec::with_capacity(new_last_node - new_first_node + 1);
            for _ in new_first_node..self.first_node {
                new_records.push(None);
            }
            for record in self.records.iter_mut() {
                new_records.push(record.take());
            }
            while new_records.len() < new_last_node - new_first_node + 1 {
                new_records.push(None);
            }
            self.first_node = new_first_node;
            self.records = new_records;
        }
    }

    // Recomputes the offsets in the outgoing edges from predecessors to the given node.
    fn recompute_outgoing_edges(
        endmarker_edges: &EdgeList,
        first_node: usize,
        records: &mut [Option<Box<MutableRecord>>],
        to: usize
    ) {
        let mut offset = 0;
        let successor = records[to - first_node].as_ref().unwrap();
        let incoming: Vec<Pos> = successor.incoming.iter().collect();
        for edge in incoming {
            if edge.node == ENDMARKER {
                // The offset is always 0 in an outgoing edge from the endmarker.
                offset += endmarker_edges.get(to).unwrap_or(0) as usize;
            } else {
                let predecessor = records[edge.node - first_node].as_mut().unwrap();
                predecessor.set_edge_offset(to, offset);
                offset += edge.offset;
            }
        }
    }

    // Special logic for iteration -1, with the current position at the endmarker
    // and the next position at the first node.
    fn handle_sequence_starts<'a>(&mut self, sequences: &'a mut [Sequence<'a>]) -> &'a mut [Sequence<'a>] {
        for (i, sequence) in sequences.iter_mut().enumerate() {
            sequence.curr = SmallPos::new(ENDMARKER, i);
            let start = sequence.nodes.first().copied().unwrap_or(ENDMARKER as u32) as usize;
            let offset = self.endmarker_edges.increment(start, 1) as usize;
            let pos = SmallPos::new(start, offset);
            sequence.next = pos;
            self.endmarker.push(Pos::from(pos));
        }

        let active_sequences = self.sort_sequences(sequences);
        let mut next = ENDMARKER as u32; // We never use the offsets in outgoing edges to the endmarker.
        for sequence in active_sequences.iter() {
            if sequence.next.node == next {
                continue;
            }
            next = sequence.next.node;
            let visits = self.endmarker_edges.get(next as usize).unwrap_or(0) as usize;
            self.records[next as usize - self.first_node].get_or_insert_default().set_visits(ENDMARKER, visits);
            Self::recompute_outgoing_edges(&self.endmarker_edges, self.first_node, &mut self.records, next as usize);
        }
        self.advance_sequences(active_sequences, 0);

        active_sequences
    }

    // Inserts `next.node` to record offset `curr.offset` for node `curr.node` in the GBWT.
    // Sets `next.offset` to the rank of `next.node` at `curr.offset`.
    fn update_bwts(&mut self, sequences: &mut [Sequence]) {
        // Update the BWT fragments.
        let mut i = 0;
        while i < sequences.len() {
            let curr = sequences[i].curr.node;
            let current = self.records[curr as usize - self.first_node].as_mut().unwrap();
            let mut run_iter = current.bwt.iter().copied();
            let mut remaining = SmallRun::default();
            let mut new_bwt = RunBuilder::new();
            while i < sequences.len() && sequences[i].curr.node == curr {
                let seq = &mut sequences[i];
                while new_bwt.target_len < seq.curr.offset as usize {
                    if remaining.len == 0 {
                        // This can only fail if MutableGBWT is in an inconsistent state.
                        let result = run_iter.next();
                        if result.is_none() {
                            panic!(
                                "MutableGBWT: Existing runs ended at position {} in node {}; inserted position is {}",
                                new_bwt.target_len, curr, seq.curr.offset
                            );
                        }
                        remaining = result.unwrap();
                    }
                    remaining = new_bwt.insert_source_run(remaining, seq.curr.offset as usize);
                }
                // We need to get the rank before inserting the new value.
                seq.next.offset = new_bwt.rank(seq.next.node as usize);
                new_bwt.insert_target_value(seq.next.node as usize);
                i += 1;
                current.len += 1;
            }
            if remaining.len > 0 {
                new_bwt.insert_source_run(remaining, usize::MAX);
            }
            for run in run_iter {
                new_bwt.insert_source_run(run, usize::MAX);
            }
            current.bwt = new_bwt.runs;
        }
    }

    // Adds the records and edges implied by the active sequences, if necessary.
    // Also increments the number of visits from `curr.node` in the record for `next.node`.
    fn add_records_and_edges(&mut self, sequences: &[Sequence]) {
        let mut node_pairs: Vec<(u32, u32)> = Vec::with_capacity(sequences.len());
        for sequence in sequences.iter() {
            node_pairs.push((sequence.curr.node, sequence.next.node));
        }

        // Ensure that the outgoing edges exist.
        node_pairs.sort_unstable();
        for (i, &(from, to)) in node_pairs.iter().enumerate() {
            if i > 0 && node_pairs[i - 1] == (from, to) {
                continue;
            }
            let predecessor = self.records[from as usize - self.first_node].as_mut().unwrap();
            predecessor.add_edge(to as usize);
        }

        // Ensure that successor records exist and increment the visits.
        node_pairs.sort_unstable_by_key(|&(_, to)| to);
        let mut start = 0;
        for (i, &(from, to)) in node_pairs.iter().enumerate() {
            if i + 1 >= node_pairs.len() || node_pairs[i + 1] != (from, to) {
                if to != ENDMARKER as u32 {
                    let successor = self.records[to as usize - self.first_node].get_or_insert_default();
                    successor.add_visits(from as usize, i - start + 1);
                }
                start = i + 1;
            }
        }
    }

    // Sorts the active sequences by `(next.node, curr)` and returns a new slice of active sequences.
    // Active sequences are those that have not reached the endmarker at `next`.
    fn sort_sequences<'a>(&mut self, sequences: &'a mut [Sequence<'a>]) -> &'a mut [Sequence<'a>] {
        sequences.sort_unstable();
        let mut i = 0;
        while i < sequences.len() && sequences[i].next.node == ENDMARKER as u32 {
            i += 1;
        }
        &mut sequences[i..]
    }

    // Updates the offsets of the outgoing edges from the current nodes of the active sequences.
    // Then determines the record offsets in the GBWT positions for the next nodes.
    fn determine_offsets(&mut self, sequences: &mut [Sequence]) {
        let mut next = 0; // We never use the offsets in outgoing edges to the endmarker.
        for sequence in sequences.iter() {
            if sequence.next.node == next {
                continue;
            }
            next = sequence.next.node;
            Self::recompute_outgoing_edges(&self.endmarker_edges, self.first_node, &mut self.records, next as usize);
        }

        for sequence in sequences.iter_mut() {
            let current = self.records[sequence.curr.node as usize - self.first_node].as_ref().unwrap();
            let start_offset = current.outgoing.get(sequence.next.node as usize).unwrap_or(0);
            // The offset was the rank of `next.node` at offset `curr.offset`.
            sequence.next.offset += start_offset;
        }
    }

    // Advances the sequences to the given offset, assuming that they are at the previous offset.
    // The sequences are now sorted by `curr`.
    // Determines the node but not the offset in the new successor GBWT position.
    fn advance_sequences(&mut self, sequences: &mut [Sequence], curr_offset: usize) {
        for sequence in sequences.iter_mut() {
            let next_node = sequence.nodes.get(curr_offset + 1).copied().unwrap_or(ENDMARKER as u32) as usize;
            sequence.curr = sequence.next;
            sequence.next = SmallPos::new(next_node, 0);
        }
    }

    /// Inserts a batch of sequences into the GBWT.
    ///
    /// The sequences are concatenated in the buffer, with each terminated by [`ENDMARKER`].
    /// If the GBWT is bidirectional, each sequence is assumed to be followed by its reverse complement.
    ///
    /// If the GBWT contains metadata, path names must be provided.
    /// Each name corresponds to a sequence (in unidirectional GBWT) or a pair of sequences (in bidirectional GBWT).
    ///
    /// The newly inserted sequences receive sequence and path identifiers starting from [`Self::sequences()`] and [`Self::paths()`], respectively.
    ///
    /// # Errors
    ///
    /// Returns an error, if:
    ///
    /// * There are trailing nodes after the last endmarker.
    /// * The number of sequences is odd in a bidirectional GBWT.
    /// * The number of path names is not as expected.
    ///
    /// Does not modify the GBWT if an error is returned.
    pub fn insert(&mut self, buffer: &[u32], names: Option<&[FullPathName]>) -> Result<(), String> {
        // Only these steps can fail, and the failure happens before we actually modify the metadata.
        let (mut sequences, node_range) = self.sequences_in_buffer(buffer)?;
        self.append_metadata(sequences.len(), names)?;

        if let Some((min_node, max_node)) = node_range {
            self.expand_records(min_node, max_node);
        }

        // Special logic for the endmarker at the start.
        let mut active_sequences = self.handle_sequence_starts(&mut sequences);
        self.len += buffer.len();

        // Invariants:
        // * GBWT position `curr` corresponds to `nodes[curr_offset]`.
        // * The sequences are sorted by `curr`.
        // * `nodes[0..=curr_offset]` have been inserted into the BWTs.
        // * `next.node` (`nodes[curr_offset + 1]`) should be inserted to BWT offset `curr.offset` in the record for `curr.node`.
        let mut curr_offset = 0;
        while !active_sequences.is_empty() {
            self.update_bwts(active_sequences);
            self.add_records_and_edges(active_sequences);
            active_sequences = self.sort_sequences(active_sequences);
            self.determine_offsets(active_sequences);
            curr_offset += 1;
            self.advance_sequences(active_sequences, curr_offset);
        }

        Ok(())
    }
}

impl From<MutableGBWT> for GBWT {
    fn from(builder: MutableGBWT) -> Self {
        // Header and tags.
        let mut header = Header::<GBWTPayload>::default();
        if builder.bidirectional {
            header.set(GBWTPayload::FLAG_BIDIRECTIONAL);
        }
        if builder.has_metadata() {
            header.set(GBWTPayload::FLAG_METADATA);
        }
        header.payload_mut().sequences = builder.sequences() as u64;
        header.payload_mut().size = builder.len() as u64;
        if let Some(min_node) = builder.min_node() {
            header.payload_mut().offset = (min_node - 1) as u64;
        }
        if let Some(max_node) = builder.max_node() {
            header.payload_mut().alphabet_size = (max_node + 1) as u64;
        } else if !builder.is_empty() {
            // Special case when we have only empty sequences.
            header.payload_mut().alphabet_size = 1;
        }

        // Encode the BWT, including the endmarker record.
        let mut bwt_builder = BWTBuilder::new(&builder.endmarker_edges, &builder.endmarker);
        for node_id in (header.payload().offset as usize + 1)..header.payload().alphabet_size as usize {
            if let Some(record) = builder.records[node_id - builder.first_node].as_ref() {
                bwt_builder.append(&record.outgoing, record.bwt.iter().map(|&run| Run::from(run)));
            } else {
                bwt_builder.append_empty();
            }
        }

        GBWT {
            header,
            tags: builder.tags,
            bwt: BWT::from(bwt_builder),
            endmarker: builder.endmarker,
            da_samples: Vec::new(),
            metadata: builder.metadata.map(Metadata::from),
        }
    }
}

//-----------------------------------------------------------------------------

/// Compares two GBWT indexes to ensure that they are bitwise identical.
///
/// The first index is assumed to be one built using this implementation.
/// The second is assumed to be a known correct index.
/// This does not compare document array samples, as only the C++ implementation currently builds them.
/// Does nothing if the indexes were built with different options (bidirectional or with metadata).
///
/// # Errors
///
/// Returns an error message if the indexes differ in any of the compared components.
pub fn compare_gbwts(index: &GBWT, truth: &GBWT, test_case: &str) -> Result<(), String> {
    if index.is_bidirectional() != truth.is_bidirectional() || index.has_metadata() != truth.has_metadata() {
        // The index was built with different options.
        return Ok(());
    }

    if index.header != truth.header {
        return Err(format!("GBWT headers do not match ({})", test_case));
    }
    if index.tags != truth.tags {
        return Err(format!("GBWT tags do not match ({})", test_case));
    }
    if index.bwt != truth.bwt {
        if index.bwt.len() != truth.bwt.len() {
            return Err(format!("GBWT BWT lengths do not match ({}): {} vs {}", test_case, index.bwt.len(), truth.bwt.len()));
        }
        let mut mismatching_records = 0;
        let mut first_mismatch = None;
        for (i, (record, truth_record)) in index.bwt.iter().zip(truth.bwt.iter()).enumerate() {
            if record != truth_record {
                mismatching_records += 1;
                if first_mismatch.is_none() {
                    first_mismatch = Some(i);
                }
            }
        }
        return Err(format!(
            "GBWT BWT records do not match ({}): {} mismatches out of {}, first mismatch at record {:?}",
            test_case, mismatching_records, index.bwt.len(), first_mismatch
        ));
    }
    // While the truth is derived from the BWT, the endmarker in a built index
    // was created independently.
    if index.endmarker != truth.endmarker {
        return Err(format!("GBWT endmarker records do not match ({})", test_case));
    }
    // We do not compare DA samples, as only the C++ implementation currently builds them.
    if index.metadata != truth.metadata {
        return Err(format!("GBWT metadata do not match ({})", test_case));
    }

    Ok(())
}

//-----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    use crate::Orientation;
    use crate::support;
    use simple_sds::serialize;
    use std::collections::HashSet;
    use rand::seq::SliceRandom;

    #[test]
    #[ignore]
    fn assumptions() {
        assert_eq!(std::mem::size_of::<EdgeList>(), 32, "EdgeList should be 32 bytes");
        assert_eq!(std::mem::size_of::<MutableRecord>(), 96, "MutableRecord should be 96 bytes");
        assert_eq!(std::mem::size_of::<Sequence>(), 32, "Sequence should be 32 bytes");
    }

//-----------------------------------------------------------------------------

    fn expand_runs(runs: &[SmallRun]) -> Vec<u32> {
        let mut result = Vec::new();
        for run in runs.iter() {
            result.extend(std::iter::repeat(run.value).take(run.len as usize));
        }
        result
    }

    // Target positions are (offset, value, expected rank).
    fn check_run_builder(source: &[SmallRun], target: &[(usize, usize, usize)]) {
        // Build the target runs.
        let mut builder = RunBuilder::new();
        let mut source_offset = 0;
        let mut remaining = SmallRun::default();
        for &(offset, value, rank) in target.iter() {
            while builder.target_len < offset {
                if remaining.len == 0 {
                    remaining = source[source_offset];
                    source_offset += 1;
                }
                remaining = builder.insert_source_run(remaining, offset);
            }
            assert_eq!(builder.rank(value), rank as u32, "Rank of value {} at offset {} should be {}", value, offset, rank);
            builder.insert_target_value(value);
        }
        if remaining.len > 0 {
            builder.insert_source_run(remaining, usize::MAX);
        }
        while source_offset < source.len() {
            builder.insert_source_run(source[source_offset], usize::MAX);
            source_offset += 1;
        }

        // Now determine the truth.
        let expanded_source = expand_runs(source);
        let mut truth = Vec::new();
        let mut source_offset = 0;
        let mut target_offset = 0;
        for &(offset, value, _) in target.iter() {
            while target_offset < offset {
                truth.push(expanded_source[source_offset]);
                source_offset += 1;
                target_offset += 1;
            }
            truth.push(value as u32);
            target_offset += 1;
        }
        while source_offset < expanded_source.len() {
            truth.push(expanded_source[source_offset]);
            source_offset += 1;
        }

        // Check the number of maximal runs.
        let mut prev = None;
        let mut runs = 0;
        for i in 0..truth.len() {
            if Some(truth[i]) != prev {
                prev = Some(truth[i]);
                runs += 1;
            }
        }
        assert_eq!(builder.runs.len(), runs, "Builder should have {} runs", runs);

        // Check the actual runs.
        let expanded_runs = expand_runs(&builder.runs);
        assert_eq!(expanded_runs.len(), truth.len(), "Expanded runs should have total length {}", truth.len());
        for i in 0..truth.len() {
            assert_eq!(expanded_runs[i], truth[i], "Expanded runs should match truth at offset {}", i);
        }
    }

    #[test]
    fn run_builder_empty() {
        let builder = RunBuilder::new();
        assert_eq!(builder.source_len, 0, "Run builder should start with source length 0");
        assert_eq!(builder.target_len, 0, "Run builder should start with target length 0");
        assert!(builder.runs.is_empty(), "Run builder should start with empty runs");
        assert!(builder.ranks.is_empty(), "Run builder should start with no ranks");
    }

    #[test]
    fn run_builder_insert_to_empty() {
        let source = Vec::new();
        let target = vec![
            (0, 1, 0),
            (1, 2, 0),
            (2, 1, 1),
            (3, 3, 0),
        ];
        check_run_builder(&source, &target);
    }

    #[test]
    fn run_builder_merge_runs() {
        let source = vec![
            SmallRun::new(1, 2),
            SmallRun::new(2, 3),
            SmallRun::new(2, 2),
            SmallRun::new(3, 1),
        ];
        let target = Vec::new();
        check_run_builder(&source, &target);
    }

    #[test]
    fn run_builder_insert_middle() {
        let source = vec![
            SmallRun::new(1, 4),
            SmallRun::new(2, 3),
            SmallRun::new(3, 4),
        ];
        let target = vec![
            (2, 1, 2), // Same value in the middle.
            (6, 1, 5), // Different value in the middle.
        ];
        check_run_builder(&source, &target);
    }

    #[test]
    fn run_builder_insert_same_ends() {
        let source = vec![
            SmallRun::new(1, 4),
            SmallRun::new(2, 3),
            SmallRun::new(3, 4),
        ];
        let target = vec![
            (0, 1, 0),  // Start of vector and run.
            (5, 1, 5),  // End of run.
            (9, 3, 0),  // Start of run.
            (14, 3, 5), // End of vector and run.
        ];
        check_run_builder(&source, &target);
    }

    #[test]
    fn run_builder_insert_different_ends() {
        let source = vec![
            SmallRun::new(1, 4),
            SmallRun::new(2, 3),
            SmallRun::new(3, 4),
        ];
        let target = vec![
            (0, 3, 0),  // Start of vector and run.
            (5, 3, 1),  // End of run.
            (9, 1, 4),  // Start of run.
            (14, 1, 5), // End of vector and run.
        ];
        check_run_builder(&source, &target);
    }

//-----------------------------------------------------------------------------

    #[test]
    fn mutable_record_empty() {
        let record = MutableRecord::new();
        assert_eq!(record.len(), 0, "Mutable record should start with length 0");
        assert!(record.is_empty(), "Mutable record should be empty");
        assert!(record.incoming.is_empty(), "Mutable record should start with no incoming edges");
        assert!(record.outgoing.is_empty(), "Mutable record should start with no outgoing edges");
        assert!(record.bwt.is_empty(), "Mutable record should start with empty BWT fragment");
        assert!(record.follow_path(0).is_none(), "Following path in an empty record should return None");
    }

    #[test]
    fn mutable_record_follow_path() {
        // We do not need incoming edges.
        let mut record = MutableRecord::new();
        record.add_edge(1);
        record.set_edge_offset(1, 10);
        record.add_edge(2);
        record.set_edge_offset(2, 20);
        record.add_edge(3);
        record.set_edge_offset(3, 30);
        record.bwt.push(SmallRun::new(1, 3));
        record.bwt.push(SmallRun::new(3, 4));
        record.bwt.push(SmallRun::new(2, 2));
        record.bwt.push(SmallRun::new(1, 1));
        record.len = 10;

        let expanded_runs = expand_runs(&record.bwt);
        for offset in 0..record.len() {
            let expected_rank = expanded_runs.iter().take(offset).filter(|&&value| value == expanded_runs[offset]).count();
            let expected_offset = expected_rank + 10 * expanded_runs[offset] as usize;
            let expected_pos = Pos::new(expanded_runs[offset] as usize, expected_offset);
            assert_eq!(
                record.follow_path(offset), Some(expected_pos),
                "Following path at offset {} should return ({}, {})", offset, expected_pos.node, expected_pos.offset
            );
        }
        assert!(record.follow_path(record.len).is_none(), "Following path at offset equal to length should return None");
    }

//-----------------------------------------------------------------------------

    fn expected_sequence(path: &[usize], reverse: bool) -> Vec<u32> {
        if reverse {
            let reversed = support::reverse_path(path);
            reversed.iter().copied().map(|node| node as u32).collect()
        } else {
            path.iter().copied().map(|node| node as u32).collect()
        }
    }

    fn total_length(paths: &[Vec<usize>], bidirectional: bool) -> usize {
        let mut len = paths.iter().map(|path| path.len() + 1).sum();
        if bidirectional {
            len *= 2;
        }
        len
    }

    // This assumes that the flags in the builder are correct.
    fn check_mutable_gbwt(builder: &MutableGBWT, paths: &[Vec<usize>], path_names: &[FullPathName], test_case: &str) {
        let len = total_length(paths, builder.is_bidirectional());
        assert_eq!(builder.len(), len, "MutableGBWT should have correct length ({})", test_case);

        let mut nodes: HashSet<usize> = HashSet::new();
        for path in paths.iter() {
            for &node in path.iter() {
                nodes.insert(node);
                if builder.is_bidirectional() {
                    nodes.insert(support::flip_node(node));
                }
            }
        }

        let sequence_count = if builder.is_bidirectional() { paths.len() * 2 } else { paths.len() };
        assert_eq!(builder.sequences(), sequence_count, "MutableGBWT should have correct sequence count ({})", test_case);
        assert_eq!(builder.paths(), paths.len(), "MutableGBWT should have correct path count ({})", test_case);

        // Check sequences and path names.
        for (i, path) in paths.iter().enumerate() {
            if builder.is_bidirectional() {
                let expected_forward = expected_sequence(path, false);
                assert_eq!(builder.sequence(2 * i), Some(expected_forward), "MutableGBWT should return correct forward sequence for path {} ({})", i, test_case);
                let expected_reverse = expected_sequence(path, true);
                assert_eq!(builder.sequence(2 * i + 1), Some(expected_reverse), "MutableGBWT should return correct reverse sequence for path {} ({})", i, test_case);
            } else {
                let expected = expected_sequence(path, false);
                assert_eq!(builder.sequence(i), Some(expected), "MutableGBWT should return correct sequence for path {} ({})", i, test_case);
            }

            if builder.has_metadata() {
                let expected_name = path_names[i].clone();
                assert_eq!(builder.path_name(i), Some(expected_name), "MutableGBWT should return correct path name for path {} ({})", i, test_case);
            } else {
                assert_eq!(builder.path_name(i), None, "MutableGBWT should not have path names ({})", test_case);
            }
        }

        // Check out-of-bounds queries.
        assert_eq!(builder.sequence(sequence_count), None, "Extracting sequence with out-of-bounds id should return None ({})", test_case);
        assert_eq!(builder.path_name(paths.len()), None, "Getting path name with out-of-bounds id should return None ({})", test_case);
    }

    fn check_built_gbwt(
        index: &GBWT,
        paths: &[Vec<usize>], path_names: &[FullPathName],
        bidirectional: bool, with_metadata: bool,test_case: &str
    ) {
        // We need to check separately that the flags survived the conversion.
        assert_eq!(index.is_bidirectional(), bidirectional, "GBWT should have correct bidirectional flag ({})", test_case);
        assert_eq!(index.has_metadata(), with_metadata, "GBWT should have correct metadata flag ({})", test_case);
        let metadata = index.metadata();

        let len = total_length(paths, index.is_bidirectional());
        assert_eq!(index.len(), len, "GBWT should have correct length ({})", test_case);

        let sequence_count = if index.is_bidirectional() { paths.len() * 2 } else { paths.len() };
        assert_eq!(index.sequences(), sequence_count, "GBWT should have correct sequence count ({})", test_case);

        // Check sequences and path names.
        for (i, path) in paths.iter().enumerate() {
            if index.is_bidirectional() {
                let fw_iter = index.sequence(support::encode_path(i, Orientation::Forward));
                assert!(fw_iter.is_some(), "GBWT should return Some for forward sequence of path {} ({})", i, test_case);
                let fw_sequence: Vec<usize> = fw_iter.unwrap().collect();
                assert_eq!(&fw_sequence, path, "GBWT should return correct forward sequence for path {} ({})", i, test_case);
                let rev_iter = index.sequence(support::encode_path(i, Orientation::Reverse));
                assert!(rev_iter.is_some(), "GBWT should return Some for reverse sequence of path {} ({})", i, test_case);
                let rev_sequence: Vec<usize> = rev_iter.unwrap().collect();
                let expected_reverse = support::reverse_path(path);
                assert_eq!(rev_sequence, expected_reverse, "GBWT should return correct reverse sequence for path {} ({})", i, test_case);
            } else {
                let fw_iter = index.sequence(i);
                assert!(fw_iter.is_some(), "GBWT should return Some for sequence of path {} ({})", i, test_case);
                let fw_sequence: Vec<usize> = fw_iter.unwrap().collect();
                assert_eq!(&fw_sequence, path, "GBWT should return correct sequence for path {} ({})", i, test_case);
            }

            if let Some(metadata) = metadata {
                let path_name = FullPathName::from_metadata(metadata, i);
                let expected_name = path_names[i].clone();
                assert_eq!(path_name, Some(expected_name), "GBWT should return correct path name for path {} ({})", i, test_case);
            }
        }

        // Check out-of-bounds queries.
        assert!(index.sequence(sequence_count).is_none(), "Extracting sequence with out-of-bounds id should return None ({})", test_case);
        if let Some(metadata) = metadata {
            assert!(FullPathName::from_metadata(metadata, paths.len()).is_none(), "Getting path name with out-of-bounds id should return None ({})", test_case);
        }
    }

    fn extract_test_case(gbwt_name: &'static str) -> (GBWT, Vec<Vec<usize>>, Vec<FullPathName>) {
        let filename = support::get_test_data(gbwt_name);
        let index: GBWT = serialize::load_from(&filename).unwrap();

        let mut paths = Vec::new();
        let mut names = Vec::new();
        let metadata = index.metadata().unwrap();
        for path_id in 0..metadata.paths() {
            let path: Vec<usize> = index.sequence(support::encode_path(path_id, Orientation::Forward)).unwrap().collect();
            paths.push(path);
            names.push(FullPathName::from_metadata(metadata, path_id).unwrap());
        }

        (index, paths, names)
    }

    // Generates all subpaths of the given path, in random order.
    fn subpath_test_case(path: &[usize]) -> Vec<Vec<usize>> {
        let mut result = Vec::new();
        for start in 0..path.len() {
            for end in start..=path.len() {
                result.push(path[start..end].to_vec());
            }
        }
        let mut rng = rand::rng();
        result.shuffle(&mut rng);
        result
    }

    fn generate_batches<'a>(
        paths: &[Vec<usize>], path_names: &'a [FullPathName],
        two_batches: bool, bidirectional: bool, with_metadata: bool
    ) -> Vec<(Vec<u32>, Option<&'a [FullPathName]>)> {
        let path_sources = if two_batches {
            let mid = paths.len() / 2;
            vec![&paths[..mid], &paths[mid..]]
        } else {
            vec![paths]
        };
        let name_sources = if with_metadata {
            if two_batches {
                let mid = path_names.len() / 2;
                vec![Some(&path_names[..mid]), Some(&path_names[mid..])]
            } else {
                vec![Some(path_names)]
            }
        } else {
            vec![None; path_sources.len()]
        };

        let mut result = Vec::new();
        for (&path_source, &name_source) in path_sources.iter().zip(name_sources.iter()) {
            let mut buffer: Vec<u32> = Vec::new();
            for path in path_source.iter() {
                buffer.extend(path.iter().copied().map(|node| node as u32));
                buffer.push(ENDMARKER as u32);
                if bidirectional {
                    let reversed = support::reverse_path(path);
                    buffer.extend(reversed.iter().copied().map(|node| node as u32));
                    buffer.push(ENDMARKER as u32);
                }
            }
            result.push((buffer, name_source));
        }

        result
    }

//-----------------------------------------------------------------------------

    #[test]
    fn mutable_gbwt_empty() {
        for (bidirectional, with_metadata) in [(false, false), (false, true), (true, false), (true, true)] {
            let test_case = format!(
                "bidirectional={}, with_metadata={}",
                bidirectional, with_metadata
            );
            let builder = MutableGBWT::new(bidirectional, with_metadata);

            // Check that flags and tags get initialized correctly.
            assert_eq!(builder.is_bidirectional(), bidirectional, "MutableGBWT should have correct bidirectional flag ({})", test_case);
            assert_eq!(builder.has_metadata(), with_metadata, "MutableGBWT should have correct metadata flag ({})", test_case);
            assert_eq!(builder.tags().len(), 1, "MutableGBWT should start with one tag ({})", test_case);
            let expected_value = String::from(SOURCE_VALUE);
            assert_eq!(builder.tags().get(SOURCE_KEY), Some(&expected_value), "MutableGBWT should have the correct source tag ({})", test_case);

            check_mutable_gbwt(&builder, &[], &[], &test_case);
            let index = GBWT::from(builder);
            check_built_gbwt(&index, &[], &[], bidirectional, with_metadata, &test_case);
        }
    }

    #[test]
    fn mutable_gbwt_insert() {
        let (original, paths, names) = extract_test_case("example.gbwt");
        for (two_batches, bidirectional, with_metadata) in [
            (false, false, false),
            (false, false, true),
            (false, true, false),
            (false, true, true),
            (true, false, false),
            (true, false, true),
            (true, true, false),
            (true, true, true),
        ] {
            let test_case = format!(
                "two_batches={}, bidirectional={}, with_metadata={}",
                two_batches, bidirectional, with_metadata
            );

            // Build the GBWT and validate it.
            let batches = generate_batches(&paths, &names, two_batches, bidirectional, with_metadata);
            let mut builder = MutableGBWT::new(bidirectional, with_metadata);
            for (i, (buffer, path_names)) in batches.iter().enumerate() {
                let result = builder.insert(buffer, *path_names);
                assert!(result.is_ok(), "Builder failed with batch {} ({})", i, test_case);
            }
            check_mutable_gbwt(&builder, &paths, &names, &test_case);
            let index = GBWT::from(builder);
            check_built_gbwt(&index, &paths, &names, bidirectional, with_metadata, &test_case);
            let result = compare_gbwts(&index, &original, &test_case);
            assert!(result.is_ok(), "{}", result.unwrap_err());
        }
    }

    #[test]
    fn mutable_gbwt_subpaths() {
        let paths = subpath_test_case(&[2, 4, 6, 8, 10]);
        for (two_batches, bidirectional) in [(false, false), (false, true), (true, false), (true, true)] {
            let test_case = format!(
                "two_batches={}, bidirectional={}",
                two_batches, bidirectional
            );

            // Build the GBWT and validate it.
            let batches = generate_batches(&paths, &[], two_batches, bidirectional, false);
            let mut builder = MutableGBWT::new(bidirectional, false);
            for (i, (buffer, _)) in batches.iter().enumerate() {
                let result = builder.insert(buffer, None);
                assert!(result.is_ok(), "Builder failed with batch {} ({})", i, test_case);
            }
            check_mutable_gbwt(&builder, &paths, &[], &test_case);
            let index = GBWT::from(builder);
            check_built_gbwt(&index, &paths, &[], bidirectional, false, &test_case);
        }
    }

//-----------------------------------------------------------------------------

    #[test]
    fn gbwt_builder_empty() {
        for (bidirectional, with_metadata) in [(false, false), (false, true), (true, false), (true, true)] {
            let test_case = format!(
                "bidirectional={}, with_metadata={}",
                bidirectional, with_metadata
            );
            let builder = GBWTBuilder::new(bidirectional, with_metadata, 1000);
            let result = builder.build();
            assert!(result.is_ok(), "Empty builder failed ({}): {}", test_case, result.unwrap_err());
            let index = result.unwrap();
            check_built_gbwt(&index, &[], &[], bidirectional, with_metadata, &test_case);
        }
    }

    #[test]
    fn gbwt_builder_insert() {
        let (original, paths, names) = extract_test_case("example.gbwt");

        // Three buffer sizes: 0 (automatic rezize), 24 (multiple batches), 1000 (single batch).
        for (bidirectional, with_metadata, buffer_size) in [
            (false, false, 0), (false, false, 24), (false, false, 1000),
            (false, true, 0), (false, true, 24), (false, true, 1000),
            (true, false, 0), (true, false, 24), (true, false, 1000),
            (true, true, 0), (true, true, 24), (true, true, 1000),
        ] {
            let test_case = format!(
                "bidirectional={}, with_metadata={}, buffer_size={}",
                bidirectional, with_metadata, buffer_size
            );
            let mut builder = GBWTBuilder::new(bidirectional, with_metadata, buffer_size);
            for (i, (path, name)) in paths.iter().zip(names.iter()).enumerate() {
                let name = if with_metadata { Some(name.clone()) } else { None };
                let result = builder.insert(path, name);
                assert!(result.is_ok(), "Builder failed with path {} ({}): {}", i, test_case, result.unwrap_err());
            }

            let result = builder.build();
            assert!(result.is_ok(), "Build failed ({}): {}", test_case, result.unwrap_err());
            let index = result.unwrap();
            check_built_gbwt(&index, &paths, &names, bidirectional, with_metadata, &test_case);
            let result = compare_gbwts(&index, &original, &test_case);
            assert!(result.is_ok(), "{}", result.unwrap_err());
        }
    }

    #[test]
    fn gbwt_builder_insert_empty() {
        let mut builder = GBWTBuilder::new(false, false, 1000);
        let result = builder.insert(&[], None);
        assert!(result.is_ok(), "Builder failed with empty path: {}", result.unwrap_err());
        let result = builder.build();
        assert!(result.is_ok(), "Build failed with empty path: {}", result.unwrap_err());
        let index = result.unwrap();
        check_built_gbwt(&index, &[Vec::new()], &[], false, false, "empty path");
    }

    #[test]
    fn gbwt_builder_subpaths() {
        let paths = subpath_test_case(&[2, 4, 6, 8, 10]);
        for (bidirectional, buffer_size) in [(false, 0), (false, 24), (false, 1000), (true, 0), (true, 24), (true, 1000)] {
            let test_case = format!(
                "bidirectional={}, buffer_size={}",
                bidirectional, buffer_size
            );
            let mut builder = GBWTBuilder::new(bidirectional, false, buffer_size);
            for (i, path) in paths.iter().enumerate() {
                let result = builder.insert(path, None);
                assert!(result.is_ok(), "Builder failed with path {} ({}): {}", i, test_case, result.unwrap_err());
            }

            let result = builder.build();
            assert!(result.is_ok(), "Build failed ({}): {}", test_case, result.unwrap_err());
            let index = result.unwrap();
            check_built_gbwt(&index, &paths, &[], bidirectional, false, &test_case);
        }
    }
}

//-----------------------------------------------------------------------------
