//! Support structures for GBWT and GBZ.

use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Vector, Access, Push, BitVec, Select};
use simple_sds::serialize::Serialize;
use simple_sds::sparse_vector::SparseVector;
use simple_sds::bits;

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::collections::btree_map::Iter as TagIter;
use std::convert::TryFrom;
use std::io::{Error, ErrorKind};
use std::iter::FusedIterator;
use std::path::PathBuf;
use std::str::Utf8Error;
use std::{cmp, io};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Returns the GBWT node identifier corresponding to the given original node and orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
///
/// # Arguments
///
/// * `id`: Identifier of the original node.
/// * `is_reverse`: Is the original node in reverse orientation.
///
/// # Panics
///
/// May panic if `id > usize::MAX / 2`.
#[inline]
pub fn encode_node(id: usize, is_reverse: bool) -> usize {
    2 * id + (is_reverse as usize)
}

/// Returns the original node identifier corresponding to the given GBWT node.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn node_id(id: usize) -> usize {
    id / 2
}

/// Returns `true` if the given GBWT node corresponds to an original node in reverse orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn node_is_reverse(id: usize) -> bool {
    id & 1 != 0
}

/// Returns the GBWT node identifier for the same original node in the other orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn flip_node(id: usize) -> usize {
    id ^ 1
}

/// Returns the sequence identifier corresponding to the given path and orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
///
/// # Arguments
///
/// * `id`: Identifier of the path.
/// * `is_reverse`: Is the path in reverse orientation.
///
/// # Panics
///
/// May panic if `id > usize::MAX / 2`.
#[inline]
pub fn encode_path(id: usize, is_reverse: bool) -> usize {
    2 * id + (is_reverse as usize)
}

/// Returns the path identifier corresponding to the given sequence.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn path_id(id: usize) -> usize {
    id / 2
}

/// Returns `true` if the given sequence corresponds to a path in reverse orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn path_is_reverse(id: usize) -> bool {
    id & 1 != 0
}

/// Returns the sequence identifier for the same path in the other orientation.
///
/// This encoding is used in bidirectional GBWT indexes.
#[inline]
pub fn flip_path(id: usize) -> usize {
    id ^ 1
}

//-----------------------------------------------------------------------------

/// Returns the full file name for a specific test file.
pub fn get_test_data(filename: &'static str) -> PathBuf {
    let mut buf = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    buf.push("test-data");
    buf.push(filename);
    buf
}

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

    // Returns (bytes to packed, packed to bytes, packed character width).
    fn alphabet(data: &[u8]) -> (Vec<usize>, Vec<u8>, usize) {
        // Determine the byte values that are present.
        let mut bytes_to_packed: Vec<usize> = vec![0; 1 << 8];
        for byte in data {
            bytes_to_packed[*byte as usize] = 1;
        }

        // Determine alphabet size.
        let sigma = bytes_to_packed.iter().sum();
        let width = bits::bit_len(cmp::max(sigma, 1) as u64 - 1);

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

        (bytes_to_packed, packed_to_bytes, width)
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
        let (pack, alphabet, width) = Self::alphabet(&self.strings);
        alphabet.serialize(writer)?;

        // Pack and serialize the strings.
        let mut packed = IntVector::new(width).unwrap();
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
            return Err(Error::new(ErrorKind::InvalidData, "StringArray: First string does not start at offset 0"));
        }
        Ok(StringArray {
            index: index,
            strings: strings,
        })
    }

    fn size_in_elements(&self) -> usize {
        let sv = SparseVector::try_from_iter(self.index.iter().take(self.len()).map(|x| x as usize)).unwrap();
        let (_, alphabet, width) = Self::alphabet(&self.strings);

        sv.size_in_elements() + alphabet.size_in_elements() + IntVector::size_by_params(self.strings.len(), width)
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
        StringArray::from(v.as_slice())
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

/// An immutable set of immutable strings with integer identifiers.
///
/// The strings are stored in a [`StringArray`] and the identifiers are indexes into the array.
///
/// A `Dictionary` can be built from a [`StringArray`] or a [`Vec`] or a slice of any type that can be converted to a string slice.
/// The construction will fail if the source contains duplicate strings.
///
/// # Examples
///
/// ```
/// use gbwt::support::Dictionary;
/// use std::convert::TryFrom;
///
/// let source = vec!["first", "second", "third", "fourth"];
/// let dict = Dictionary::try_from(source.as_slice()).unwrap();
/// for (index, value) in source.iter().enumerate() {
///     assert_eq!(dict.id(value), Some(index));
///     assert_eq!(dict.bytes(index), source[index].as_bytes());
/// }
/// assert_eq!(dict.id("fifth"), None);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Dictionary {
    strings: StringArray,
    sorted_ids: IntVector,
}

impl Dictionary {
    /// Returns the number of strings in the dictionary.
    #[inline]
    pub fn len(&self) -> usize {
        self.strings.len()
    }

    /// Returns `true` if the dictionary is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the identifier of the given string in the dictionary, or [`None`] if there is no such string.
    pub fn id<T: AsRef<[u8]>>(&self, string: T) -> Option<usize> {
        let mut low = 0;
        let mut high = self.len();
        while low < high {
            let mid = low + (high - low) / 2;
            let id = self.sorted_ids.get(mid) as usize;
            match string.as_ref().cmp(self.bytes(id)) {
                Ordering::Less => high = mid,
                Ordering::Equal => return Some(id),
                Ordering::Greater => low = mid + 1,
            }
        }
        None
    }

    /// Returns a byte slice corresponding to the string with identifier `i`.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn bytes(&self, i: usize) -> &[u8] {
        self.strings.bytes(i)
    }

    /// Returns a string slice corresponding to the string with identifier `i` or an error if the bytes are not valid UTF-8.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn str(&self, i: usize) -> Result<&str, Utf8Error> {
        self.strings.str(i)
    }

    /// Returns a copy of the string with identifier `i` or an error if the bytes are not valid UTF-8.
    ///
    /// # Panics
    ///
    /// May panic if `i >= self.len()`.
    pub fn string(&self, i: usize) -> Result<String, Utf8Error> {
        self.strings.string(i)
    }
}

impl Serialize for Dictionary {
    fn serialize_header<T: io::Write>(&self, _: &mut T) -> io::Result<()> {
        Ok(())
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.strings.serialize(writer)?;
        self.sorted_ids.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let strings = StringArray::load(reader)?;
        let sorted_ids = IntVector::load(reader)?;
        Ok(Dictionary {
            strings: strings,
            sorted_ids: sorted_ids,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.strings.size_in_elements() + self.sorted_ids.size_in_elements()
    }
}

impl TryFrom<StringArray> for Dictionary {
    type Error = &'static str;

    fn try_from(source: StringArray) -> Result<Self, Self::Error> {
        // Sort the ids and check for duplicates.
        let mut sorted: Vec<usize> = Vec::with_capacity(source.len());
        for i in 0..source.len() {
            sorted.push(i);
        }
        sorted.sort_unstable_by(|a, b| source.bytes(*a).cmp(source.bytes(*b)));
        for i in 1..sorted.len() {
            if source.bytes(sorted[i - 1]) == source.bytes(sorted[i]) {
                return Err("Cannot build a dictionary from a source with duplicate strings");
            }
        }

        // Compact the sorted ids.
        let width = if sorted.is_empty() { 1 } else { bits::bit_len(sorted.len() as u64 - 1) };
        let mut sorted_ids = IntVector::with_capacity(sorted.len(), width).unwrap();
        sorted_ids.extend(sorted);

        Ok(Dictionary {
            strings: source,
            sorted_ids: sorted_ids,
        })
    }
}

impl<T: AsRef<str>> TryFrom<&[T]> for Dictionary {
    type Error = &'static str;

    fn try_from(source: &[T]) -> Result<Self, Self::Error> {
        Self::try_from(StringArray::from(source))
    }
}

impl<T: AsRef<str>> TryFrom<Vec<T>> for Dictionary {
    type Error = &'static str;

    fn try_from(source: Vec<T>) -> Result<Self, Self::Error> {
        Self::try_from(StringArray::from(source))
    }
}

impl AsRef<StringArray> for Dictionary {
    #[inline]
    fn as_ref(&self) -> &StringArray {
        &(self.strings)
    }
}

//-----------------------------------------------------------------------------

/// A key-value structure with strings as both keys and values.
///
/// The keys are case insensitive.
/// This structure is a simple wrapper over [`BTreeMap`]`<`[`String`]`, `[`String`]`>` that converts all keys to lower case.
///
/// # Examples
///
/// ```
/// use gbwt::support::Tags;
///
/// let mut tags = Tags::new();
/// tags.insert("first-key", "first-value");
/// tags.insert("second-key", "second-value");
/// assert!(tags.contains_key("First-Key"));
/// assert_eq!(*tags.get("second-key").unwrap(), "second-value");
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Tags {
    tags: BTreeMap<String, String>,
}

impl Tags {
    /// Creates an empty `Tags` structure.
    pub fn new() -> Tags {
        Tags::default()
    }

    /// Returns the number of tags.
    pub fn len(&self) -> usize {
        self.tags.len()
    }

    /// Returns `true` if the structure is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the value corresponding to the key, or [`None`] no such tag exists.
    pub fn get(&self, key: &str) -> Option<&String> {
        let key = key.to_lowercase();
        self.tags.get(&key)
    }

    /// Returns `true` if there is a tag with the given key.
    pub fn contains_key(&self, key: &str) -> bool {
        let key = key.to_lowercase();
        self.tags.contains_key(&key)
    }

    /// Inserts a new tag, overwriting the possible old value associated with the same key.
    ///
    /// # Arguments
    ///
    /// * `key`: Key of the tag. The key is converted to lower case before it is inserted into the hash table.
    /// * `value`: Value of the tag.
    pub fn insert(&mut self, key: &str, value: &str) {
        let key = key.to_lowercase();
        let _ = self.tags.insert(key, value.to_string());
    }

    /// Returns an iterator that visits all tags in sorted order by keys.
    ///
    /// The type of `Item` is `(&`[`String`]`, &`[`String`]`)`.
    pub fn iter(&self) -> TagIter<'_, String, String> {
        self.tags.iter()
    }

    // Returns the array of keys and values in serialized order.
    fn linearize(&self) -> StringArray {
        let mut linearized: Vec<&str> = Vec::with_capacity(2 * self.len());
        for (key, value) in self.iter() {
            linearized.push(key); linearized.push(value);
        }
        StringArray::from(linearized)
    }
}

impl Serialize for Tags {
    fn serialize_header<T: io::Write>(&self, _: &mut T) -> io::Result<()> {
        Ok(())
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        let linearized = self.linearize();
        linearized.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let linearized = StringArray::load(reader)?;
        if linearized.len() % 2 != 0 {
            return Err(Error::new(ErrorKind::InvalidData, "Tags: Key without a value"));
        }
        let mut result = Tags::new();
        for i in 0..linearized.len() / 2 {
            let key = linearized.str(2 * i).map_err(|_| Error::new(ErrorKind::InvalidData, "Tags: Invalid UTF-8 in a key"))?;
            let value = linearized.str(2 * i + 1).map_err(|_| Error::new(ErrorKind::InvalidData, "Tags: Invalid UTF-8 in a value"))?;
            result.insert(key, value);
        }
        if result.len() != linearized.len() / 2 {
            return Err(Error::new(ErrorKind::InvalidData, "Tags: Duplicate keys"));
        }
        Ok(result)
    }

    fn size_in_elements(&self) -> usize {
        let linearized = self.linearize();
        linearized.size_in_elements()
    }
}

impl AsRef<BTreeMap<String, String>> for Tags {
    #[inline]
    fn as_ref(&self) -> &BTreeMap<String, String> {
        &(self.tags)
    }
}

impl From<ByteCode> for Vec<u8> {
    fn from(source: ByteCode) -> Self {
        source.bytes
    }
}

//-----------------------------------------------------------------------------

/// A variable-length encoder for unsigned integers.
///
/// `ByteCode` encodes an integer as a sequence of bytes in little-endian order and stores it in the internal [`Vec`].
/// Each byte contains 7 bits of data, and the high bit indicates whether the encoding continues.
/// The bytes can be accessed with [`AsRef`] or extracted with [`From`], and [`ByteCodeIter`] can be used for decoding the integers.
/// Raw bytes can be appended to the encoding using [`ByteCode::write_byte`].
///
/// # Examples
///
/// ```
/// use gbwt::support::ByteCode;
///
/// let mut encoder = ByteCode::new();
/// encoder.write(123); encoder.write(456); encoder.write(789);
/// let bytes = encoder.as_ref();
/// assert_eq!(*bytes, [123, 72 + 128, 3, 21 + 128, 6]);
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ByteCode {
    bytes: Vec<u8>
}

impl ByteCode {
    const MASK: u8 = 0x7F;
    const FLAG: u8 = 0x80;
    const SHIFT: usize = 7;

    /// Creates a new encoder.
    pub fn new() -> Self {
        ByteCode::default()
    }

    /// Encodes `value` and stores the encoding.
    pub fn write(&mut self, value: usize) {
        let mut value = value;
        while value > (Self::MASK as usize) {
            self.bytes.push(((value as u8) & Self::MASK) | Self::FLAG);
            value >>= Self::SHIFT;
        }
        self.bytes.push(value as u8);
    }

    /// Appends a byte to the encoding.
    pub fn write_byte(&mut self, byte: u8) {
        self.bytes.push(byte);
    }

    /// Returns the total number of bytes in the encoding.
    #[inline]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Returns `true` if the encoding is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl AsRef<[u8]> for ByteCode {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.bytes
    }
}

/// An iterator that decodes integers from a byte slice encoded by [`ByteCode`].
///
/// The type of `Item` is [`usize`].
/// Raw bytes can be read from the encoding using [`ByteCodeIter::byte`].
///
/// # Examples
///
/// ```
/// use gbwt::support::{ByteCode, ByteCodeIter};
///
/// let mut source = ByteCode::new();
/// source.write(123); source.write(456); source.write(789);
///
/// let mut iter = ByteCodeIter::new(source.as_ref());
/// assert_eq!(iter.next(), Some(123));
/// assert_eq!(iter.next(), Some(456));
/// assert_eq!(iter.next(), Some(789));
/// assert_eq!(iter.next(), None);
/// ```
#[derive(Clone, Debug)]
pub struct ByteCodeIter<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> ByteCodeIter<'a> {
    /// Returns an iterator over the byte slice.
    pub fn new(bytes: &'a [u8]) -> Self {
        ByteCodeIter {
            bytes: bytes,
            offset: 0,
        }
    }

    /// Returns the next byte from the slice, or `None` if there are no more bytes left.
    pub fn byte(&mut self) -> Option<u8> {
        if self.offset >= self.bytes.len() {
            return None;
        }
        let result = Some(self.bytes[self.offset]);
        self.offset += 1;
        result
    }

    /// Returns the first unvisited offset in the byte slice.
    #[inline]
    pub fn offset(&self) -> usize {
        self.offset
    }
}

impl<'a> Iterator for ByteCodeIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let mut offset = 0;
        let mut result = 0;
        while self.offset < self.bytes.len() {
            let value = unsafe { *self.bytes.get_unchecked(self.offset) };
            self.offset += 1;
            result += ((value & ByteCode::MASK) as usize) << offset;
            offset += ByteCode::SHIFT;
            if value & ByteCode::FLAG == 0 {
                return Some(result);
            }
        }
        None
    }
}

impl<'a> FusedIterator for ByteCodeIter<'a> {}

//-----------------------------------------------------------------------------

/// A run-length encoder for non-empty runs of unsigned integers.
///
/// The exact encoding depends on alphabet size `sigma`.
/// If `sigma` is small, the encoder tries to encode short runs as a single byte.
/// For long runs, the remaining run length is encoded using [`ByteCode`].
/// For a large `sigma`, both the value and the run length are encoded using [`ByteCode`].
/// Alphabet size `sigma == 0` indicates a large alphabet of unknown size.
///
/// The bytes can be accessed with [`AsRef`] or extracted with [`From`], and [`RLEIter`] can be used for decoding the integers.
/// Raw bytes and [`ByteCode`]-encoded integers can be appended to the encoding using [`RLE::write_byte`] and [`RLE::write_int`].
/// The following functions can for creating a byte stream with various encodings:
///
/// # Examples
///
/// ```
/// use gbwt::support::RLE;
///
/// let mut encoder = RLE::with_sigma(4);
/// encoder.write(3, 12); encoder.write(2, 721); encoder.write(0, 34);
/// assert_eq!(*encoder.as_ref(), [3 + 4 * 11, 2 + 4 * 63, 17 + 128, 5, 0 + 4 * 33]);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RLE {
    bytes: ByteCode,
    sigma: usize,
    threshold: usize,
}

impl RLE {
    const THRESHOLD: usize = 255;
    const UNIVERSE: usize = 256;

    /// Creates a new encoder with alphabet size `0`.
    pub fn new() -> Self {
        RLE::default()
    }

    /// Creates a new encoder with the given alphabet size.
    pub fn with_sigma(sigma: usize) -> Self {
        let (sigma, threshold) = Self::sanitize(sigma);
        RLE {
            bytes: ByteCode::new(),
            sigma: sigma,
            threshold: threshold,
        }
    }

    /// Encodes a run of `len` copies of value `value` and stores the encoding.
    ///
    /// Does nothing if `len == 0`.
    ///
    /// # Panics
    ///
    /// Panics if `value >= self.sigma()`.
    pub fn write(&mut self, value: usize, len: usize) {
        if len == 0 {
            return;
        }
        if value >= self.sigma {
            panic!("RLE: Cannot encode value {} with alphabet size {}", value, self.sigma);
        }
        unsafe { self.write_unchecked(value, len); }
    }

    /// Encodes a run of `len` copies of value `value` and stores the encoding.
    ///
    /// Behavior is undefined if `len == 0` or `value >= self.sigma()`.
    pub unsafe fn write_unchecked(&mut self, value: usize, len: usize) {
        if self.sigma >= Self::THRESHOLD {
            self.bytes.write(value);
            self.bytes.write(len - 1);
        } else if len < self.threshold {
            self.write_basic(value, len);
        } else {
            self.write_basic(value, self.threshold);
            self.bytes.write(len - self.threshold);
        }
    }

    /// Appends a byte to the encoding.
    pub fn write_byte(&mut self, byte: u8) {
        self.bytes.write_byte(byte);
    }

    /// Encodes `value` using [`ByteCode`] and stores the encoding.
    pub fn write_int(&mut self, value: usize) {
        self.bytes.write(value);
    }

    /// Returns the total number of bytes in the encoding.
    #[inline]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Returns `true` if the encoding is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the alphabet size.
    #[inline]
    pub fn sigma(&self) -> usize {
        self.sigma
    }

    /// Changes the alphabet size to `sigma`.
    pub fn set_sigma(&mut self, sigma: usize) {
        let (sigma, threshold) = Self::sanitize(sigma);
        self.sigma = sigma;
        self.threshold = threshold;
    }

    // Writes a single-byte run.
    fn write_basic(&mut self, value: usize, len: usize) {
        let code = value + self.sigma * (len - 1);
        self.bytes.write_byte(code as u8);
    }

    // Returns (effective sigma, threshold for short runs).
    fn sanitize(sigma: usize) -> (usize, usize) {
        let sigma = if sigma == 0 { usize::MAX } else { sigma };
        let threshold = if sigma < Self::THRESHOLD { Self::UNIVERSE / sigma } else { 0 };
        (sigma, threshold)
    }
}

impl Default for RLE {
    fn default() -> Self {
        let (sigma, threshold) = Self::sanitize(0);
        RLE {
            bytes: ByteCode::new(),
            sigma: sigma,
            threshold: threshold,
        }
    }
}

impl AsRef<[u8]> for RLE {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        self.bytes.as_ref()
    }
}

impl From<RLE> for Vec<u8> {
    fn from(source: RLE) -> Self {
        Self::from(source.bytes)
    }
}

//-----------------------------------------------------------------------------

/// An iterator that decodes runs from a byte slice encoded by [`RLE`].
///
/// The type of `Item` is `(`[`usize`]`, `[`usize`]`)`.
/// Alphabet size `sigma == 0` indicates a large alphabet of unknown size.
/// Raw bytes and [`ByteCode`]-encoded integers can be read from the encoding using [`RLEIter::byte`] and [`RLEIter::int`].
///
/// # Examples
///
/// ```
/// use gbwt::support::{RLE, RLEIter};
///
/// let mut source = RLE::with_sigma(4);
/// source.write(3, 12); source.write(2, 721); source.write(0, 34);
///
/// let mut iter = RLEIter::with_sigma(source.as_ref(), 4);
/// assert_eq!(iter.next(), Some((3, 12)));
/// assert_eq!(iter.next(), Some((2, 721)));
/// assert_eq!(iter.next(), Some((0, 34)));
/// assert_eq!(iter.next(), None);
/// ```
#[derive(Clone, Debug)]
pub struct RLEIter<'a> {
    source: ByteCodeIter<'a>,
    sigma: usize,
    threshold: usize,
}

impl<'a> RLEIter<'a> {
    /// Creates a new iterator over the byte slice with alphabet size `0`.
    pub fn new(bytes: &'a [u8]) -> Self {
        let (sigma, threshold) = RLE::sanitize(0);
        RLEIter {
            source: ByteCodeIter::new(bytes),
            sigma: sigma,
            threshold: threshold,
        }
    }

    /// Creates a new iterator.
    ///
    /// # Arguments
    ///
    /// * `bytes`: Byte slice.
    /// * `sigma`: Alphabet size.
    pub fn with_sigma(bytes: &'a [u8], sigma: usize) -> Self {
        let (sigma, threshold) = RLE::sanitize(sigma);
        RLEIter {
            source: ByteCodeIter::new(bytes),
            sigma: sigma,
            threshold: threshold,
        }
    }

    /// Returns the next byte from the slice, or `None` if there are no more bytes left.
    pub fn byte(&mut self) -> Option<u8> {
        self.source.byte()
    }

    /// Returns the next [`ByteCode`]-encoded integer from the slice, or `None` if no more integers can be decoded.
    pub fn int(&mut self) -> Option<usize> {
        self.source.next()
    }

    /// Returns the first unvisited offset in the byte slice.
    #[inline]
    pub fn offset(&self) -> usize {
        self.source.offset()
    }

    /// Returns the alphabet size.
    #[inline]
    pub fn sigma(&self) -> usize {
        self.sigma
    }

    /// Changes the alphabet size to `sigma`.
    pub fn set_sigma(&mut self, sigma: usize) {
        let (sigma, threshold) = RLE::sanitize(sigma);
        self.sigma = sigma;
        self.threshold = threshold;
    }
}

impl<'a> Iterator for RLEIter<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let mut run = (0, 0);
        if self.sigma >= RLE::THRESHOLD {
            if let Some(value) = self.source.next() { run.0 = value; } else { return None; }
            if let Some(len) = self.source.next() { run.1 = len + 1; } else { return None; }
        } else {
            if let Some(byte) = self.source.byte() {
                run.0 = (byte as usize) % self.sigma;
                run.1 = (byte as usize) / self.sigma + 1;
            } else {
                return None;
            }
            if run.1 == self.threshold {
                if let Some(len) = self.source.next() { run.1 += len; } else { return None; }
            }
        }
        return Some(run);
    }
}

impl<'a> FusedIterator for RLEIter<'a> {}

//-----------------------------------------------------------------------------
