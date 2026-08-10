use super::*;

use crate::GBZ;

use simple_sds::{bits, serialize};

use rand::Rng;
use rand::seq::SliceRandom;
use rand::rngs::ThreadRng;

use std::collections::HashSet;
use std::fs::{self, File, OpenOptions};

//-----------------------------------------------------------------------------

#[test]
fn reverse_sequences() {
    {
        let sequence = b"";
        let truth = b"";
        let rc = reverse_complement(sequence);
        assert_eq!(&rc, truth, "Invalid reverse complement for an empty sequence");
    }

    {
        let sequence = b"C";
        let truth = b"G";
        let rc = reverse_complement(sequence);
        assert_eq!(&rc, truth, "Invalid reverse complement for a sequence of length 1");
    }

    {
        let sequence = b"GATTACA";
        let truth = b"TGTAATC";
        let rc = reverse_complement(sequence);
        assert_eq!(&rc, truth, "Invalid reverse complement for a sequence with odd length");
    }

    {
        let sequence = b"GATTACAT";
        let truth = b"ATGTAATC";
        let rc = reverse_complement(sequence);
        assert_eq!(&rc, truth, "Invalid reverse complement for a sequence with even length");
    }
}

#[test]
fn reverse_paths() {
    let mut original = vec![1, 2, 4, 6];
    let reversed = reverse_path(&original);
    assert_eq!(reversed, vec![7, 5, 3, 0], "Failed to reverse the path correctly");

    reverse_path_in_place(&mut original);
    assert_eq!(original, reversed, "Failed to reverse the path in place correctly");

    let mut buffer = Vec::new();
    reverse_path_into(&original, &mut buffer);
    reverse_path_into(&original, &mut buffer);
    assert_eq!(buffer, vec![1, 2, 4, 6, 1, 2, 4, 6], "Failed to reverse the path into a buffer correctly");
}

#[test]
fn intersections() {
    assert!(intersect(&(3..5), &(7..8)).is_empty(), "Failed when a is before b");
    assert_eq!(intersect(&(4..8), &(5..11)), 5..8, "Failed when a overlaps the start of b");
    assert_eq!(intersect(&(5..8), &(4..9)), 5..8, "Failed when a is contained in b");
    assert_eq!(intersect(&(3..12), &(7..11)), 7..11, "Failed when a contains b");
    assert_eq!(intersect(&(4..10), &(1..6)), 4..6, "Failed when a overlaps the end of b");
    assert!(intersect(&(8..13), &(1..7)).is_empty(), "Failed when a is after b");

    assert!(intersect(&(3..3), &(2..4)).is_empty(), "Failed when a is empty");
    assert!(intersect(&(2..4), &(2..2)).is_empty(), "Failed when b is empty");
    assert!(intersect(&(3..3), &(4..4)).is_empty(), "Failed when both are empty");
}

//-----------------------------------------------------------------------------

fn check_array(array: &StringArray, truth: &[&str]) {
    // Statistics.
    assert_eq!(array.len(), truth.len(), "Incorrect array length");
    assert_eq!(array.is_empty(), truth.is_empty(), "Incorrect array emptiness");
    assert_eq!(array.iter().len(), truth.len(), "Invalid iterator length");

    // Access.
    for i in 0..array.len() {
        assert_eq!(array.str_len(i), truth[i].len(), "Incorrect length for string {}", i);
        assert_eq!(array.bytes(i), truth[i].as_bytes(), "Incorrect bytes for string {}", i);
        assert_eq!(array.str(i).unwrap(), truth[i], "Incorrect string slice {}", i);
        assert_eq!(array.string(i).unwrap(), truth[i], "Incorrect string {}", i);
    }

    // Access with ranges.
    for start in 0..=array.len() {
        let mut expected_len = 0;
        for end in start..=array.len() {
            let len = array.range_len(start..end);
            let bytes = array.range(start..end);
            if start < end {
                expected_len += truth[end - 1].len();
                assert_eq!(len, expected_len, "Invalid slice length for range {}..{}", start, end);
                assert_eq!(bytes, truth[start..end].concat().as_bytes(), "Invalid slice for range {}..{}", start, end);
            } else {
                assert_eq!(len, 0, "Invalid slice length for empty range {}..{}", start, end);
                assert!(bytes.is_empty(), "Non-empty slice for empty range {}..{}", start, end);
            }
        }
    }

    // Iterate forward.
    for (index, bytes) in array.iter().enumerate() {
        assert_eq!(bytes, truth[index].as_bytes(), "Invalid bytes for string {} from iterator (forward)", index);
    }

    // Iterate backward.
    let mut next = array.len();
    let mut iter = array.iter();
    while let Some(bytes) = iter.next_back() {
        next -= 1;
        assert_eq!(bytes, truth[next].as_bytes(), "Invalid bytes for string {} from iterator (backward)", next);
    }

    // Meet in the middle.
    let mut next = 0;
    let mut limit = array.len();
    let mut iter = array.iter();
    while iter.len() > 0 {
        assert_eq!(iter.next().unwrap(), truth[next].as_bytes(), "Invalid bytes for string {} from iterator (forward, bidirectional)", next);
        next += 1;
        if iter.len() == 0 {
            break;
        }
        limit -= 1;
        assert_eq!(iter.next_back().unwrap(), truth[limit].as_bytes(), "Invalid bytes for string {} from iterator (backward, bidirectional)", limit);
    }
}

fn check_compression(array: &StringArray) {
    let filename = serialize::temp_file_name("string-array-compression");

    let mut options = OpenOptions::new();
    let file = options.write(true).create(true).open(&filename);
    assert!(file.is_ok(), "Failed to create a temporary file: {}", file.err().unwrap());
    let mut file = file.unwrap();

    let result =array.compress(&mut file, None);
    assert!(result.is_ok(), "Failed to compress the string array: {}", result.err().unwrap());

    let metadata = fs::metadata(&filename);
    assert!(metadata.is_ok(), "Failed to get metadata for the temporary file: {}", metadata.err().unwrap());
    let metadata = metadata.unwrap();
    let file_size = metadata.len() as usize;
    let reported_size = bits::words_to_bytes(array.compressed_size_in_elements(None));
    assert_eq!(reported_size, file_size, "Reported compressed size does not match the actual file size");

    let mut options = OpenOptions::new();
    let file = options.read(true).open(&filename);
    assert!(file.is_ok(), "Failed to open the temporary file: {}", file.err().unwrap());
    let mut file = file.unwrap();

    let decompressed = StringArray::decompress(&mut file);
    assert!(decompressed.is_ok(), "Failed to decompress the string array: {}", decompressed.err().unwrap());
    let decompressed = decompressed.unwrap();
    assert_eq!(&decompressed, array, "Decompressed string array does not match the original");

    fs::remove_file(&filename).unwrap();
}

#[test]
fn empty_string_array() {
    let truth: Vec<&str> = Vec::new();
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    let _ = serialize::test(&array, "empty-string-array", None, true);
    check_compression(&array);
}

#[test]
fn non_empty_string_array() {
    let truth = vec!["first", "second", "third", "fourth"];
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    let _ = serialize::test(&array, "non-empty-string-array", None, true);
    check_compression(&array);
}

#[test]
fn string_array_with_empty_strings() {
    // Serialization with an empty string at the end used to fail in the original GBWT implementation.
    let truth = vec!["first", "second", "", "fourth", ""];
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    let _ = serialize::test(&array, "string-array-with-empty", None, true);
    check_compression(&array);
}

//-----------------------------------------------------------------------------

fn check_dict(dict: &Dictionary, truth: &[&str], missing: &[&str]) {
    // Statistics.
    assert_eq!(dict.len(), truth.len(), "Incorrect dictionary length");
    assert_eq!(dict.is_empty(), truth.is_empty(), "Incorrect dictionary emptiness");

    // Id -> string.
    for i in 0..dict.len() {
        assert_eq!(dict.bytes(i), truth[i].as_bytes(), "Incorrect bytes for string {}", i);
        assert_eq!(dict.str(i).unwrap(), truth[i], "Incorrect string slice {}", i);
        assert_eq!(dict.string(i).unwrap(), truth[i], "Incorrect string {}", i);
    }

    // String -> id.
    for i in 0..truth.len() {
        assert_eq!(dict.id(truth[i]), Some(i), "Invalid id for original string {}: {}", i, truth[i]);
    }
    for string in missing.iter() {
        assert_eq!(dict.id(string), None, "String {} should not be present", string);
    }
}

#[test]
fn empty_dict() {
    let truth: Vec<&str> = Vec::new();
    let missing = vec!["this", "should", "not", "exist"];
    let dict = Dictionary::try_from(truth.as_slice()).unwrap();
    check_dict(&dict, &truth, &missing);
    let _ = serialize::test(&dict, "empty-dict", None, true);
}

#[test]
fn non_empty_dict() {
    let truth = vec!["first", "second", "third", "fourth"];
    let missing = vec!["this", "should", "not", "exist"];
    let dict = Dictionary::try_from(truth.as_slice()).unwrap();
    check_dict(&dict, &truth, &missing);
    let _ = serialize::test(&dict, "non-empty-dict", None, true);
}

#[test]
fn dict_from_duplicates() {
    let source = vec!["this", "contains", "contains", "many", "duplicates", "many", "this"];
    let result = Dictionary::try_from(source);
    assert!(result.is_err(), "Did not get an error from a source with duplicate strings");
}

//-----------------------------------------------------------------------------

fn check_tags(tags: &Tags, truth: &BTreeMap<&str, &str>, missing: &[&str]) {
    // Statistics.
    assert_eq!(tags.len(), truth.len(), "Incorrect tags length");
    assert_eq!(tags.is_empty(), truth.is_empty(), "Incorrect tags emptiness");

    // Truth is present.
    for (key, value) in truth.iter() {
        assert!(tags.contains_key(key), "Key {} is missing", key);
        assert_eq!(tags.get(key).unwrap(), value, "Invalid value for key {}", key);
    }

    // Keys and values are correct.
    for (key, value) in tags.iter() {
        assert!(truth.contains_key(key.as_str()), "Key {} is incorrect", key);
        assert_eq!(truth.get(key.as_str()).unwrap(), value, "Incorrect value for key {}", key);
    }

    // Missing keys.
    for key in missing.iter() {
        assert!(!tags.contains_key(key), "Key {} should not be present", key);
    }
}

#[test]
fn empty_tags() {
    let truth: BTreeMap<&str, &str> = BTreeMap::new();
    let missing = vec!["this", "should", "not", "exist"];
    let tags = Tags::new();
    check_tags(&tags, &truth, &missing);
    let _ = serialize::test(&tags, "empty-tags", None, true);
}

#[test]
fn non_empty_tags() {
    let mut truth: BTreeMap<&str, &str> = BTreeMap::new();
    truth.insert("first-key", "first-value");
    truth.insert("second-key", "second-value");
    truth.insert("third-key", "third-value");
    truth.insert("fourth-key", "fourth-value");
    let missing = vec!["this", "should", "not", "exist"];
    let mut tags = Tags::new();
    for (key, value) in truth.iter() {
        tags.insert(key, value);
    }
    check_tags(&tags, &truth, &missing);
    let _ = serialize::test(&tags, "non-empty-tags", None, true);
}

#[test]
fn case_insensitive_tags() {
    let mut truth: BTreeMap<&str, &str> = BTreeMap::new();
    truth.insert("first-key", "first-value");
    truth.insert("second-key", "second-value");
    truth.insert("third-key", "third-value");
    truth.insert("fourth-key", "fourth-value");
    let missing = vec!["this", "should", "not", "exist"];
    let mut tags = Tags::new();
    tags.insert("First-Key", "first-value");
    tags.insert("second-Key", "second-value");
    tags.insert("Third-key", "third-value");
    tags.insert("fourth-key", "fourth-value");
    check_tags(&tags, &truth, &missing);
    let _ = serialize::test(&tags, "case-insensitive-tags", None, true);
}

#[test]
fn duplicate_tags() {
    let mut truth: BTreeMap<&str, &str> = BTreeMap::new();
    truth.insert("first-key", "first-value");
    truth.insert("second-key", "second-value");
    truth.insert("third-key", "third-value");
    truth.insert("fourth-key", "fourth-value");
    let missing = vec!["this", "should", "not", "exist"];
    let mut tags = Tags::new();
    tags.insert("second-key", "incorrect-value");
    tags.insert("Fourth-Key", "incorrect-value");
    for (key, value) in truth.iter() {
        tags.insert(key, value);
    }
    check_tags(&tags, &truth, &missing);
    let _ = serialize::test(&tags, "duplicate-tags", None, true);
}

#[test]
fn remove_tags() {
    let mut truth: BTreeMap<&str, &str> = BTreeMap::new();
    truth.insert("first-key", "first-value");
    truth.insert("second-key", "second-value");
    truth.insert("third-key", "third-value");
    truth.insert("fourth-key", "fourth-value");
    let missing = vec!["this", "should", "not", "exist"];
    let mut tags = Tags::new();
    for (key, value) in truth.iter() {
        tags.insert(key, value);
    }
    truth.remove("second-key");
    assert_eq!(tags.remove("second-key"), Some(String::from("second-value")), "Failed to remove a tag");
    for key in missing.iter() {
        assert_eq!(tags.remove(key), None, "Failed to remove a non-existing tag");
    }
    check_tags(&tags, &truth, &missing);
    let _ = serialize::test(&tags, "remove-tags", None, true);
}

//-----------------------------------------------------------------------------

// Generate a random value, with the width (almost) geometrically distributed (p = 0.5) in blocks of `w` bits.
fn generate_value(rng: &mut ThreadRng, w: usize) -> usize {
    let len = (rng.random::<u64>() | 1).leading_zeros() as usize; // 0 to 63
    let width = cmp::min((len + 1) * w, bits::WORD_BITS);
    let mask = bits::low_set(width) as usize;
    rng.random::<u64>() as usize & mask
}

// Generate `n` random values, with the widths (almost) geometrically distributed (p = 0.5) in blocks of `w` bits.
fn generate_values(n: usize, w: usize) -> Vec<usize> {
    let mut result = Vec::with_capacity(n);
    let mut rng = rand::rng();
    for _ in 0..n {
        result.push(generate_value(&mut rng, w));
    }
    result
}

#[test]
fn random_byte_code() {
    let values = generate_values(647, 4);
    let mut encoder = ByteCode::new();
    assert_eq!(encoder.len(), 0, "Newly created encoder contains bytes");
    assert!(encoder.is_empty(), "Newly created encoder is not empty");
    for value in values.iter() {
        encoder.write(*value);
    }
    assert!(encoder.len() >= values.len(), "The encoding is shorter than the number of values");
    assert!(!encoder.is_empty(), "The encoding is empty");

    let mut iter = ByteCodeIter::new(encoder.as_ref());
    assert_eq!(iter.offset(), 0, "Newly creater iterator is not at offset 0");
    let mut i = 0;
    while let Some(value) = iter.next() {
        assert!(i < values.len(), "Too many values from the iterator");
        assert_eq!(value, values[i], "Invalid value {}", i);
        i += 1;
    }
    assert_eq!(i, values.len(), "Too few values from the iterator");
    assert_eq!(iter.offset(), encoder.len(), "Iterator did not consume all bytes");
}

//-----------------------------------------------------------------------------

// Generate `n` random runs from an alphabet of size `sigma`.
// The widths of run lengths are (almost) geometrically distributed (p = 0.5) in blocks of `w` bits.
fn generate_runs(n: usize, sigma: usize, w: usize) -> Vec<Run> {
    let sigma = if sigma == 0 { usize::MAX } else { sigma };
    let mut result = Vec::with_capacity(n);
    let mut rng = rand::rng();
    for _ in 0..n {
        let c: usize = rng.random_range(0..sigma);
        let len = generate_value(&mut rng, w) + 1;
        result.push(Run::new(c, len));
    }
    result
}

fn encode_runs(encoder: &mut RLE, runs: &[Run], name: &str) {
    assert_eq!(encoder.len(), 0, "[{}]: Newly created encoder contains runs", name);
    assert!(encoder.is_empty(), "[{}]: Newly created encoder is not empty", name);
    for run in runs {
        encoder.write(*run);
    }
    assert!(encoder.len() >= runs.len(), "[{}]: The encoding is shorter than the number of runs", name);
    assert!(!encoder.is_empty(), "[{}]: The encoding is empty", name);
}

fn check_runs(encoder: &RLE, truth: &[Run], name: &str) {
    let mut iter = RLEIter::with_sigma(encoder.as_ref(), encoder.sigma());
    assert_eq!(iter.offset(), 0, "[{}]: Newly creater iterator is not at offset 0", name);
    let mut i = 0;
    while let Some(run) = iter.next() {
        assert!(i < truth.len(), "[{}]: Too many runs from the iterator", name);
        assert_eq!(run, truth[i], "[{}]: Invalid run {}", name, i);
        i += 1;
    }
    assert_eq!(i, truth.len(), "[{}]: Too few runs from the iterator", name);
    assert_eq!(iter.offset(), encoder.len(), "[{}]: Iterator did not consume all bytes", name);
}

fn test_rle(n: usize, sigma: usize, name: &str) {
    let runs = generate_runs(n, sigma, 4);
    let mut encoder = RLE::with_sigma(sigma);
    encode_runs(&mut encoder, &runs, name);
    check_runs(&encoder, &runs, name);
}

fn add_run(encoder: &mut RLE, truth: &mut Vec<Run>, len: usize, bytes: usize, name: &str) {
    let old_len = encoder.len();
    encoder.write(Run::new(encoder.sigma() - 1, len));
    truth.push(Run::new(encoder.sigma() - 1, len));
    assert_eq!(encoder.len() - old_len, bytes, "[{}]: Run of length {} not encoded using {} byte(s)", name, len, bytes);
}

fn test_threshold(sigma: usize, name: &str) {
    let (sigma, threshold) = RLE::sanitize(sigma);
    let mut encoder = RLE::with_sigma(sigma);
    let mut truth: Vec<Run> = Vec::new();
    if threshold > 1 {
        add_run(&mut encoder, &mut truth, threshold - 1, 1, name);
    }
    if threshold > 0 {
        add_run(&mut encoder, &mut truth, threshold, 2, name);
    }
    check_runs(&encoder, &truth, name);
}

#[test]
fn runs_with_sigma() {
    test_rle(591, 4, "sigma == 4");
    test_rle(366, 254, "sigma == 254");
    test_rle(421, 255, "sigma == 255");
    test_rle(283, 14901, "sigma == 14901");
    test_rle(330, 0, "sigma == 0");
}

#[test]
fn run_length_thresholds() {
    test_threshold(1, "sigma == 1");
    test_threshold(4, "sigma == 4");
    test_threshold(5, "sigma == 5");
    test_threshold(128, "sigma == 128");
    test_threshold(129, "sigma == 129");
    test_threshold(254, "sigma == 254");
}

#[test]
fn gbwt_record() {
    // Original data for the record.
    let sigma = 4;
    let edges: Vec<(usize, usize)> = vec![(0, 0), (13, 7), (22, 1), (44, 0)];
    let runs = generate_runs(8, sigma, 4);

    // Encode the record.
    let mut encoder = RLE::new();
    encoder.write_int(sigma);
    let mut prev = 0;
    for (node, offset) in edges.iter() {
        encoder.write_int(*node - prev); encoder.write_int(*offset);
        prev = *node;
    }
    encoder.set_sigma(sigma);
    for run in runs.iter() {
        encoder.write(*run);
    }

    // Decompress the edges.
    let mut iter = RLEIter::new(encoder.as_ref());
    assert_eq!(iter.int(), Some(sigma), "Invalid alphabet size in the record");
    let mut prev = 0;
    for i in 0..sigma {
        let node = iter.int().unwrap() + prev;
        assert_eq!(node, edges[i].0, "Invalid successor node {}", i);
        prev = node;
        assert_eq!(iter.int(), Some(edges[i].1), "Invalid record offset for edge {}", i);
    }

    // Decompress the runs.
    iter.set_sigma(sigma);
    let mut decoded: Vec<Run> = Vec::new();
    while let Some(run) = iter.next() {
        decoded.push(run);
    }
    assert_eq!(decoded.len(), runs.len(), "Invalid number of runs");
    for i in 0..decoded.len() {
        assert_eq!(decoded[i], runs[i], "Invalid run {}", i);
    }
    assert_eq!(iter.offset(), encoder.len(), "Iterator did not consume all bytes");
}

//-----------------------------------------------------------------------------

#[test]
fn empty_disjoint_sets() {
    let mut sets = DisjointSets::new(0, 0);
    assert_eq!(sets.len(), 0, "Empty structure has non-zero length");
    assert!(sets.is_empty(), "Empty structure is not empty");
    assert_eq!(sets.offset(), 0, "Empty structure has non-zero offset");

    let sets = sets.extract(|_| true);
    assert!(sets.is_empty(), "Empty structure contains sets");
}

fn random_sets(len: usize, offset: usize, num_sets: usize, rng: &mut ThreadRng) -> Vec<Vec<usize>> {
    let mut values: Vec<usize> = (offset..len + offset).collect();
    values.shuffle(rng);

    let mut result: Vec<Vec<usize>> = Vec::new();
    for &value in values[0..num_sets].iter() {
        result.push(vec![value]);
    }
    for &value in values[num_sets..len].iter() {
        let set = rng.random::<u64>() as usize % num_sets;
        result[set].push(value);
    }

    result
}

fn join_sets(sets: &mut DisjointSets, source: &Vec<Vec<usize>>, rng: &mut ThreadRng) {
    for set in source.iter() {
        for right in 1..set.len() {
            let left = rng.random::<u64>() as usize % right;
            sets.union(set[left], set[right]);
        }
    }
}

fn filter_sets<F: Fn(usize) -> bool>(source: Vec<Vec<usize>>, include_value: F) -> Vec<Vec<usize>> {
    let mut result: Vec<Vec<usize>> = Vec::new();
    for set in source.iter() {
        let filtered: Vec<usize> = set.iter().filter(|value| include_value(**value)).copied().collect();
        if !filtered.is_empty() {
            result.push(filtered);
        }
    }
    result
}

fn sort_sets(sets: &mut Vec<Vec<usize>>) {
    for set in sets.iter_mut() {
        set.sort();
    }
    sets.sort();
}

#[test]
fn zero_offset_disjoint_sets() {
    let len = 38;
    let offset = 0;
    let num_sets = 5;

    let mut sets = DisjointSets::new(len, offset);
    assert_eq!(sets.len(), len, "Invalid length");
    assert_eq!(sets.is_empty(), len == 0, "Invalid emptiness");
    assert_eq!(sets.offset(), offset, "Invalid offset");

    let mut rng = rand::rng();
    let mut source = random_sets(len, offset, num_sets, &mut rng);
    join_sets(&mut sets, &source, &mut rng);
    sort_sets(&mut source);

    let extracted = sets.extract(|_| true);
    assert_eq!(extracted.len(), num_sets, "Invalid number of sets");
    for i in 0..num_sets {
        assert_eq!(extracted[i], source[i], "Invalid set {}", i);
    }
}

#[test]
fn non_zero_offset_disjoint_sets() {
    let len = 41;
    let offset = 22;
    let num_sets = 6;

    let mut sets = DisjointSets::new(len, offset);
    assert_eq!(sets.len(), len, "Invalid length");
    assert_eq!(sets.is_empty(), len == 0, "Invalid emptiness");
    assert_eq!(sets.offset(), offset, "Invalid offset");

    let mut rng = rand::rng();
    let mut source = random_sets(len, offset, num_sets, &mut rng);
    join_sets(&mut sets, &source, &mut rng);
    sort_sets(&mut source);

    let extracted = sets.extract(|_| true);
    assert_eq!(extracted.len(), num_sets, "Invalid number of sets");
    for i in 0..num_sets {
        assert_eq!(extracted[i], source[i], "Invalid set {}", i);
    }
}

#[test]
fn filtered_disjoint_sets() {
    let len = 52;
    let offset = 21;
    let num_sets = 7;

    let mut sets = DisjointSets::new(len, offset);
    assert_eq!(sets.len(), len, "Invalid length");
    assert_eq!(sets.is_empty(), len == 0, "Invalid emptiness");
    assert_eq!(sets.offset(), offset, "Invalid offset");

    let mut rng = rand::rng();
    let source = random_sets(len, offset, num_sets, &mut rng);
    join_sets(&mut sets, &source, &mut rng);
    let mut source = filter_sets(source, |value| value % 3 != 0);
    sort_sets(&mut source);

    let extracted = sets.extract(|value| value % 3 != 0);
    assert_eq!(extracted.len(), num_sets, "Invalid number of sets");
    for i in 0..num_sets {
        assert_eq!(extracted[i], source[i], "Invalid set {}", i);
    }
}

//-----------------------------------------------------------------------------

#[test]
fn edge_list_empty() {
    let edges = EdgeList::new();
    assert!(edges.is_empty(), "Edge list should be empty");
    assert_eq!(edges.len(), 0, "Edge list should have length 0");
    assert_eq!(edges.get(1), None, "Edge list should not contain any edges");
    assert_eq!(edges.iter().count(), 0, "Edge list iterator should be empty");
}

fn create_edge_list_and_sort(truth: &mut [Pos]) -> EdgeList {
    let mut edges = EdgeList::new();
    for edge in truth.iter() {
        edges.increment(edge.node, edge.offset);
    }
    truth.sort();
    edges
}

// Assumes that the truth is sorted and non-empty.
fn check_edge_list(edges: &EdgeList, truth: &[Pos]) {
    assert!(!edges.is_empty(), "Edge list should not be empty");
    assert_eq!(edges.len(), truth.len(), "Edge list should have length {}", truth.len());

    let mut max_node = 0;
    for edge in truth.iter() {
        let offset = edge.offset as u32;
        assert_eq!(edges.get(edge.node), Some(offset), "Edge list should contain edge ({}, {})", edge.node, edge.offset);
        if edge.node > max_node {
            max_node = edge.node;
        }
    }
    assert!(edges.get(max_node + 1).is_none(), "Edge list should not contain edge ({}, _)", max_node + 1);

    // Increment the offsets in a copy and check them again.
    let mut copy = edges.clone();
    for edge in truth.iter() {
        let offset = copy.get_mut(edge.node);
        assert!(offset.is_some(), "Edge list should contain edge ({}, _)", edge.node);
        let offset = offset.unwrap();
        *offset += 1;
    }
    for edge in truth.iter() {
        let offset = edge.offset as u32 + 1;
        assert_eq!(copy.get(edge.node), Some(offset), "Incremented edge list should contain edge ({}, {})", edge.node, offset);
    }
}

#[test]
fn edge_list_small() {
    let mut truth = vec![
        Pos::new(3, 10),
        Pos::new(1, 5),
        Pos::new(2, 7),
    ];
    let edges = create_edge_list_and_sort(&mut truth);
    check_edge_list(&edges, &truth);
}

#[test]
fn edge_list_large() {
    let mut truth = vec![
        Pos::new(3, 10),
        Pos::new(1, 5),
        Pos::new(2, 7),
        Pos::new(4, 12),
        Pos::new(10, 18),
    ];
    let edges = create_edge_list_and_sort(&mut truth);
    check_edge_list(&edges, &truth);
}

#[test]
fn edge_list_increment() {
    let mut truth = vec![
        Pos::new(3, 10),
        Pos::new(1, 6),
        Pos::new(2, 8),
        Pos::new(4, 12),
        Pos::new(10, 18),
    ];
    let mut edges = EdgeList::new();
    for edge in truth.iter() {
        edges.increment(edge.node, edge.offset / 2);
        edges.increment(edge.node, edge.offset / 2);
    }
    truth.sort();
    check_edge_list(&edges, &truth);
}

#[test]
fn edge_list_ranks() {
    let mut truth = vec![
        Pos::new(3, 2),
        Pos::new(1, 0),
        Pos::new(2, 1),
        Pos::new(4, 3),
        Pos::new(10, 4),
    ];
    let mut edges = EdgeList::new();
    for edge in truth.iter() {
        // Initialize the edge list with arbitrary offsets.
        edges.increment(edge.node, 10);
    }

    edges.set_ranks();
    truth.sort();
    check_edge_list(&edges, &truth);
}

//-----------------------------------------------------------------------------

fn load_chains(filename: &PathBuf) -> Vec<IntVector> {
    let file = File::open(&filename);
    assert!(file.is_ok(), "Failed to open chains file {}", filename.display());
    let mut file = file.unwrap();
    let data = Chains::read_data(&mut file);
    assert!(data.is_ok(), "Failed to read chains from {}", filename.display());
    data.unwrap()
}

fn check_chains(chains: &Chains, data: &[IntVector]) {
    assert_eq!(chains.len(), data.len(), "Wrong number of chains");
    assert!(chains.trivial_chains().is_none(), "Expected an unknown number of trivial chains");
    assert!(chains.components().is_none(), "Expected an unknown number of components");
    let mut expected_links = 0;
    for chain in data.iter() {
        if chain.len() > 1 {
            expected_links += chain.len() - 1;
        }
    }
    assert_eq!(chains.links(), expected_links, "Wrong number of links");

    // Handles and nodes should be present.
    for chain in data.iter() {
        for handle in chain.iter().map(|x| x as usize) {
            assert!(chains.has_handle(handle), "Missing handle {}", handle);
            let rev_handle = flip_node(handle);
            assert!(chains.has_handle(rev_handle), "Missing reverse handle {}", rev_handle);
            let (node_id, _) = decode_node(handle);
            assert!(chains.has_node(node_id), "Missing node {}", node_id);
        }
    }

    // Missing handles and nodes should not be present.
    let max_handle = data.iter().flat_map(|x| x.iter().map(|y| y as usize)).max().unwrap();
    assert!(!chains.has_handle(max_handle + 2), "Unexpected handle {}", max_handle + 1);
    let missing_node = decode_node(max_handle).0 + 1;
    assert!(!chains.has_node(missing_node), "Unexpected node {}", missing_node);
}

fn check_region(graph: &GBZ, chains: &Chains, from: usize, to: usize) {
    // Active handles. We proceed to their successors but not predecessors.
    let mut active = vec![from, flip_node(to)];
    // Visited node identifiers.
    let mut visited = HashSet::new();
    visited.insert(decode_node(from).0);
    visited.insert(decode_node(to).0);

    while !active.is_empty() {
        let curr = active.pop().unwrap();
        let (node_id, orientation) = decode_node(curr);
        for (next_id, next_o) in graph.successors(node_id, orientation).unwrap() {
            if visited.contains(&next_id) {
                continue;
            }
            let fw_handle = encode_node(next_id, next_o);
            let rev_handle = flip_node(fw_handle);
            assert!(!chains.has_node(next_id), "Reached boundary node {} ({}, {}) from region {}..{}", fw_handle, next_id, next_o, from, to);
            active.push(fw_handle); active.push(rev_handle);
            visited.insert(next_id);
        }
    }
}

#[test]
fn chains_empty() {
    let chains = Chains::new();
    assert_eq!(chains.len(), 0, "Expected empty chains");
    assert!(chains.trivial_chains().is_none(), "Expected an unknown number of trivial chains");
    assert!(chains.components().is_none(), "Expected an unknown number of components");
    assert_eq!(chains.links(), 0, "Expected no links");
    let _ = serialize::test(&chains, "empty-chains", None, true);
}

#[test]
fn chains_nonempty() {
    let chains_file = get_test_data("micb-kir3dl1.chains");
    let data = load_chains(&chains_file);
    let chains = serialize::load_from(&chains_file);
    if let Err(msg) = chains {
        panic!("Failed to read chains from {}: {}", chains_file.display(), msg);
    }
    let chains: Chains = chains.unwrap();
    check_chains(&chains, &data);

    let graph_file = get_test_data("micb-kir3dl1.gbz");
    let graph: GBZ = serialize::load_from(&graph_file).unwrap();
    for (from, to) in chains.iter() {
        check_region(&graph, &chains, from, to);
    }

    let _ = serialize::test(&chains, "nonempty-chains", None, true);
}

#[test]
fn chains_add_link() {
    let chains_file = get_test_data("micb-kir3dl1.chains");
    let data = load_chains(&chains_file);
    let truth = serialize::load_from(&chains_file);
    if let Err(msg) = truth {
        panic!("Failed to read chains from {}: {}", chains_file.display(), msg);
    }
    let truth: Chains = truth.unwrap();

    // Build the chains manually.
    let mut chains = Chains::new();
    for (i, chain) in data.iter().enumerate() {
        for j in 1..chain.len() {
            let from = chain.get(j - 1) as usize;
            let to = chain.get(j) as usize;
            let result = chains.add_link(from, to);
            assert!(result, "Failed to add link from {} to {} (chain {}, index {})", from, to, i, j);
        }
    }
    chains.count_chains();
    assert_eq!(chains.len(), truth.len(), "Mismatch in number of chains");
    assert!(chains.trivial_chains().is_none(), "Expected an unknown number of trivial chains");
    assert!(chains.components().is_none(), "Expected an unknown number of components");
    assert_eq!(chains.links(), truth.links(), "Mismatch in number of links");
    assert_eq!(chains, truth, "Mismatch in chains content");

    // Try adding duplicate links in either orientation.
    const NO_NODE: usize = 1_000_000;
    for (i, chain) in data.iter().enumerate() {
        for (j, handle) in chain.iter().enumerate() {
            let node = handle as usize;
            if j + 1 < chain.len() {
                let result = chains.add_link(node, NO_NODE);
                assert!(!result, "Added a duplicate link from {} (chain {}, index {})", node, i, j);
            }
            if j > 0 {
                let result = chains.add_link(flip_node(node), NO_NODE);
                assert!(!result, "Added a duplicate link from {} (chain {}, index {})", flip_node(node), i, j);
            }
        }
    }
}

//-----------------------------------------------------------------------------
