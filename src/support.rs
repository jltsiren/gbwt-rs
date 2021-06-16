//! Support structures for GBWT and GBZ.

use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Vector, Access, Push, BitVec, Select};
use simple_sds::serialize::Serialize;
use simple_sds::sparse_vector::SparseVector;
use simple_sds::bits;

use std::io::{Error, ErrorKind};
use std::iter::FusedIterator;
use std::str::Utf8Error;
use std::{cmp, io};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// An immutable array of immutable strings.
///
/// The strings are concatenated and stored in a single byte vector.
/// This reduces the space overhead for the strings and the time overhead for serializing and loading them.
/// The serialization format further compresses the starting positions and compacts the alphabet in an attempt to use fewer than 8 bits per byte.
///
/// `StringArray` can be built from a [`Vec`] or a slice of any type that can be converted to a string slice.
/// Construction from an iterator is not feasible, as `StringArray` needs to know the total length of the strings in advance.
///
/// Because the bytes may come from an untrusted source, `StringArray` does not assume that the bytes are valid UTF-8 strings.
///
/// # Examples
///
/// ```
/// use gbwt::support::StringArray;
///
/// let source = vec!["first", "second", "third", "fourth"];
/// let array = StringArray::from(source.as_slice());
/// assert_eq!(array.len(), source.len());
/// for i in 0..array.len() {
///     assert_eq!(array.str(i).unwrap(), source[i]);
/// }
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StringArray {
    index: IntVector,
    strings: Vec<u8>,
}

impl StringArray {
    /// Returns the number of strings in the array.
    #[inline]
    pub fn len(&self) -> usize {
        self.index.len() - 1
    }

    /// Returns `true` if the array is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the length of the `i`th string in bytes.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn str_len(&self, i: usize) -> usize {
        (self.index.get(i + 1) - self.index.get(i)) as usize
    }

    /// Returns a byte slice corresponding to the `i`th string.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn bytes(&self, i: usize) -> &[u8] {
        let start = self.index.get(i) as usize;
        let limit = self.index.get(i + 1) as usize;
        &self.strings[start..limit]
    }

    /// Returns a string slice corresponding to the `i`th string or an error if the bytes are not valid UTF-8.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn str(&self, i: usize) -> Result<&str, Utf8Error> {
        std::str::from_utf8(self.bytes(i))
    }

    /// Returns a copy of the `i`th string or an error if the bytes are not valid UTF-8.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn string(&self, i: usize) -> Result<String, Utf8Error> {
        match self.str(i) {
            Ok(v) => Ok(v.to_string()),
            Err(e) => Err(e),
        }
    }

    /// Returns an iterator over the string array.
    pub fn iter(&self) -> StringIter<'_> {
        StringIter {
            parent: self,
            next: 0,
            limit: self.len(),
        }
    }

    // Builds an empty string array with capacity for `n` strings of total length `total_len`.
    fn with_capacity(n: usize, total_len: usize) -> StringArray {
        let mut index = IntVector::with_capacity(n + 1, bits::bit_len(total_len as u64)).unwrap();
        index.push(0);
        let strings: Vec<u8> = Vec::with_capacity(total_len);
        StringArray {
            index: index,
            strings: strings,
        }
    }

    // Appends a new string to the array, assuming that there is space for it.
    fn append(&mut self, string: &str) {
        self.strings.extend(string.bytes());
        self.index.push(self.strings.len() as u64);
    }

    // Returns (bytes to packed, packed to bytes).
    fn alphabet(data: &[u8]) -> (Vec<usize>, Vec<u8>) {
        // Determine the byte values that are present.
        let mut bytes_to_packed: Vec<usize> = vec![0; 1 << 8];
        for byte in data {
            bytes_to_packed[*byte as usize] = 1;
        }

        // Determine alphabet size.
        let sigma = cmp::max(bytes_to_packed.iter().sum(), 1);

        // Build the alphabet mappings.
        let mut packed_to_bytes: Vec<u8> = vec![0; sigma];
        let mut rank = 0;
        for i in 0..bytes_to_packed.len() {
            if bytes_to_packed[i] != 0 {
                bytes_to_packed[i] = rank;
                packed_to_bytes[rank] = i as u8;
                rank += 1;
            }
        }

        (bytes_to_packed, packed_to_bytes)
    }
}

impl Serialize for StringArray {
    fn serialize_header<T: io::Write>(&self, _: &mut T) -> io::Result<()> {
        Ok(())
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        // Compress the index without the past-the-end sentinel.
        let sv = SparseVector::try_from_iter(self.index.iter().take(self.len()).map(|x| x as usize)).unwrap();
        sv.serialize(writer)?;
        drop(sv);

        // Determine and serialize the alphabet
        let (pack, alphabet) = Self::alphabet(&self.strings);
        alphabet.serialize(writer)?;

        // Pack and serialize the strings.
        let mut packed = IntVector::new(bits::bit_len((alphabet.len() - 1) as u64)).unwrap();
        packed.extend(self.strings.iter().map(|x| pack[*x as usize]));
        packed.serialize(writer)?;

        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        // Load the compressed index. We need the strings for the past-the-end sentinel.
        let sv = SparseVector::load(reader)?;

        // Load the alphabet.
        let alphabet = Vec::<u8>::load(reader)?;

        // Load and decompress the strings.
        let packed = IntVector::load(reader)?;
        let strings: Vec<u8> = packed.into_iter().map(|x| alphabet[x as usize]).collect();

        // Decompress the index.
        let mut index = IntVector::with_capacity(sv.count_ones() + 1, bits::bit_len(strings.len() as u64)).unwrap();
        index.extend(sv.one_iter().map(|(_, x)| x));
        index.push(strings.len() as u64);

        // Sanity checks.
        if index.get(0) != 0 {
            return Err(Error::new(ErrorKind::InvalidData, "First string does not start at offset 0"));
        }
        Ok(StringArray {
            index: index,
            strings: strings,
        })
    }

    fn size_in_elements(&self) -> usize {
        let sv = SparseVector::try_from_iter(self.index.iter().take(self.len()).map(|x| x as usize)).unwrap();
        let (_, alphabet) = Self::alphabet(&self.strings);

        sv.size_in_elements() + alphabet.size_in_elements() + IntVector::size_by_params(self.strings.len(), bits::bit_len(self.strings.len() as u64))
    }
}

impl<T: AsRef<str>> From<&[T]> for StringArray {
    fn from(v: &[T]) -> Self {
        let total_len = v.iter().fold(0, |sum, item| sum + item.as_ref().len());
        let mut result = StringArray::with_capacity(v.len(), total_len);
        for string in v.iter() {
            result.append(string.as_ref());
        }
        result
    }
}

impl<T: AsRef<str>> From<Vec<T>> for StringArray {
    fn from(v: Vec<T>) -> Self {
        let total_len = v.iter().fold(0, |sum, item| sum + item.as_ref().len());
        let mut result = StringArray::with_capacity(v.len(), total_len);
        for string in v.iter() {
            result.append(string.as_ref());
        }
        result
    }
}

//-----------------------------------------------------------------------------

/// A read-only iterator over [`StringArray`].
///
/// The type of `Item` is `&[`[`u8`]`]`.
///
/// # Examples
///
/// ```
/// use gbwt::support::StringArray;
/// use std::str;
///
/// let source = vec!["first", "second", "third"];
/// let array = StringArray::from(source.as_slice());
/// for (index, bytes) in array.iter().enumerate() {
///     assert_eq!(bytes, source[index].as_bytes());
/// }
/// ```
#[derive(Clone, Debug)]
pub struct StringIter<'a> {
    parent: &'a StringArray,
    // The first index we have not used.
    next: usize,
    // The first index we should not use.
    limit: usize,
}

impl<'a> Iterator for StringIter<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.limit {
            None
        } else {
            let result = Some(self.parent.bytes(self.next));
            self.next += 1;
            result
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.limit - self.next;
        (remaining, Some(remaining))
    }
}

impl<'a> DoubleEndedIterator for StringIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.next >= self.limit {
            None
        } else {
            self.limit -= 1;
            Some(self.parent.bytes(self.limit))
        }
    }
}

impl<'a> ExactSizeIterator for StringIter<'a> {}

impl<'a> FusedIterator for StringIter<'a> {}

//-----------------------------------------------------------------------------

// FIXME Dictionary, Tags, ByteCode, Run

//-----------------------------------------------------------------------------
