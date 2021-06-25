//! The BWT stored as an array of compressed node records.

use crate::support::ByteCodeIter;

use simple_sds::sparse_vector::SparseVector;
use simple_sds::ops::{BitVec, Select};
use simple_sds::serialize::Serialize;

use std::io::{Error, ErrorKind};
use std::io;

//-----------------------------------------------------------------------------

// FIXME document, example, tests
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BWT {
    index: SparseVector,
    data: Vec<u8>,
}

impl BWT {
    /// Returns the number of records in the BWT.
    #[inline]
    pub fn len(&self) -> usize {
        self.index.count_ones()
    }

    /// Returns `true` if the BWT is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the `i`th record or `None` if there is no such node.
    pub fn record(&self, i: usize) -> Option<Record> {
        if i >= self.len() {
            return None;
        }
        let mut iter = self.index.select_iter(i);
        let (_, start) = iter.next().unwrap();
        let limit = if i + 1 >= self.len() { iter.next().unwrap().1 } else { self.data.len() };
        Record::new(&self.data[start..limit])
    }
}

impl Serialize for BWT {
    fn serialize_header<T: io::Write>(&self, _: &mut T) -> io::Result<()> {
        Ok(())
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.index.serialize(writer)?;
        self.data.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let index = SparseVector::load(reader)?;
        let data = Vec::<u8>::load(reader)?;
        if index.len() != data.len() {
            return Err(Error::new(ErrorKind::InvalidData, "BWT: Index / data length mismatch"));
        }
        Ok(BWT {
            index: index,
            data: data,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.index.size_in_elements() + self.data.size_in_elements()
    }
}

//-----------------------------------------------------------------------------

// FIXME document, example, tests
#[derive(Clone, Debug)]
pub struct Record<'a> {
    edges: Vec<(usize, usize)>,
    bwt: &'a [u8],
}

impl<'a> Record<'a> {
    /// Returns the outdegree of the node.
    #[inline]
    pub fn outdegree(&self) -> usize {
        self.edges.len()
    }

    /// Returns a record corresponding to the byte slice, or `None` if the record would be empty.
    pub fn new(bytes: &'a [u8]) -> Option<Self> {
        if bytes.is_empty() {
            return None;
        }

        // Determine the outdegree.
        let mut iter = ByteCodeIter::new(bytes);
        let sigma = iter.next().unwrap();
        if sigma == 0 {
            return None;
        }

        // Decompress the edges.
        let mut edges: Vec<(usize, usize)> = Vec::new();
        let mut prev = 0;
        for _ in 0..sigma {
            let node = iter.next().unwrap() + prev;
            prev = node;
            let offset = iter.next().unwrap();
            edges.push((node, offset));
        }

        Some(Record {
            edges: edges,
            bwt: &bytes[iter.offset()..],
        })
    }
}

//-----------------------------------------------------------------------------
