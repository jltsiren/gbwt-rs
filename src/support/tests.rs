use super::*;

use simple_sds::serialize;

use rand::Rng;
use rand::seq::SliceRandom;
use rand::rngs::ThreadRng;

use std::fs::{self, OpenOptions};

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
