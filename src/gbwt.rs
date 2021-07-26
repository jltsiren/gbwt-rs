// FIXME top-level documentation, example

use crate::{ENDMARKER, SOURCE_KEY, SOURCE_VALUE};
use crate::bwt::BWT;
use crate::headers::{Header, GBWTPayload};
use crate::support::Tags;

use simple_sds::serialize::Serialize;

use std::io::{Error, ErrorKind};
use std::io;

//-----------------------------------------------------------------------------

// FIXME document, tests
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GBWT {
    header: Header<GBWTPayload>,
    tags: Tags,
    bwt: BWT,
}

// FIXME tests
// Statistics.
impl GBWT {
    /// Returns the total length of the sequences in the index.
    pub fn len(&self) -> usize {
        self.header.payload().size
    }

    /// Returns `true` if the index is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of sequences in the index.
    pub fn sequences(&self) -> usize {
        self.header.payload().sequences
    }

    /// Returns the size of the alphabet.
    pub fn alphabet_size(&self) -> usize {
        self.header.payload().alphabet_size
    }

    /// Returns the alphabet offset for the effective alphabet.
    pub fn alphabet_offset(&self) -> usize {
        self.header.payload().offset
    }

    /// Returns the size of the effective alphabet.
    pub fn effective_size(&self) -> usize {
        self.alphabet_size() - self.alphabet_offset()
    }

    /// Returns the smallest node identifier in the effective alphabet.
    pub fn first_node(&self) -> usize {
        self.alphabet_offset() + 1
    }

    /// Returns `true` if node identifier `id` is in the effective alphabet.
    pub fn has_node(&self, id: usize) -> bool {
        id > self.alphabet_offset() && id < self.alphabet_size()
    }

    /// Returns `true` if the GBWT index is bidirectional.
    pub fn is_bidirectional(&self) -> bool {
        self.header.is_set(GBWTPayload::FLAG_BIDIRECTIONAL)
    }
}

// FIXME default? new?

//-----------------------------------------------------------------------------

// FIXME tests
// Sequence navigation.
impl GBWT {
    /// Returns the first position in sequence `sequence`, or [`None`] if no such sequence exists.
    ///
    /// The return value is a pair (node identifier, offset in node).
    pub fn start(&self, sequence: usize) -> Option<(usize, usize)> {
        if let Some(record) = self.bwt.record(ENDMARKER) {
            return record.lf(sequence);
        }
        None
    }

    /// Follows the sequence forward and returns the next position, or [`None`] if no such position exists.
    ///
    /// The argument and the return value are pairs (node identifier, offset in node).
    pub fn forward(&self, pos: (usize, usize)) -> Option<(usize, usize)> {
        // This also catches the endmarker.
        if pos.0 <= self.first_node() {
            return None;
        }
        if let Some(record) = self.bwt.record(pos.0 - self.alphabet_offset()) {
            return record.lf(pos.1);
        }
        None
    }

    /// Follows the sequence backward and returns the previous position, or [`None`] if no such position exists.
    ///
    /// The argument and the return value are pairs (node identifier, offset in node).
    ///
    /// # Panics
    ///
    /// Panics if the index is not bidirectional.
    pub fn backward(&self, pos: (usize, usize)) -> Option<(usize, usize)> {
        if !self.is_bidirectional() {
            panic!("Following sequences backward is only possible in a bidirectional GBWT");
        }
        // This also catches the endmarker.
        if pos.0 <= self.first_node() {
            return None;
        }
        if let Some(record) = self.bwt.record(pos.0 - self.alphabet_offset()) {
            // FIXME count the number of occurrences of each successor of the reverse node
            // FIXME this allows us to find the predecessor of the current node
            // FIXME note the special case if the right successor is (v, forward) but there are also edges to (v, reverse)
            // FIXME because orientation flipping changes the order
            // FIXME then go to the predecessor and determine which position goes to the given offset
        }
        None
    }

    // FIXME iter(sequence)
}

//-----------------------------------------------------------------------------

// FIXME impl: find, extend, bd_find, extend_forward, extend_backward

//-----------------------------------------------------------------------------

// FIXME tests
impl Serialize for GBWT {
    fn serialize_header<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.header.serialize(writer)
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.tags.serialize(writer)?;
        self.bwt.serialize(writer)?;
        // FIXME da samples, metadata
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
        // FIXME we should decompress the endmarker
        // FIXME skip da samples, metadata
        Ok(GBWT {
            header: header,
            tags: tags,
            bwt: bwt,
        })
    }

    fn size_in_elements(&self) -> usize {
        // FIXME da samples, metadata
        self.header.size_in_elements() + self.tags.size_in_elements() + self.bwt.size_in_elements()
    }
}

//-----------------------------------------------------------------------------

// FIXME SearchState, BDSearchState

//-----------------------------------------------------------------------------

// FIXME Iter

//-----------------------------------------------------------------------------
