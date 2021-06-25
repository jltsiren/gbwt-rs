use super::*;

use simple_sds::serialize;

use rand::Rng;
use rand::rngs::ThreadRng;

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

#[test]
fn empty_string_array() {
    let truth: Vec<&str> = Vec::new();
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    let _ = serialize::test(&array, "empty-string-array", None, true);
}

#[test]
fn non_empty_string_array() {
    let truth = vec!["first", "second", "third", "fourth"];
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    let _ = serialize::test(&array, "non-empty-string-array", None, true);
}

#[test]
fn array_with_empty_strings() {
    // Serialization with an empty string at the end used to fail in the original GBWT implementation.
    let truth = vec!["first", "second", "", "fourth", ""];
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    let _ = serialize::test(&array, "string-array-with-empty", None, true);
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

//-----------------------------------------------------------------------------

// Generate a random value, with the width (almost) geometrically distributed (p = 0.5) in blocks of `w` bits.
fn generate_value(rng: &mut ThreadRng, w: usize) -> usize {
    let len = (rng.gen::<usize>() | 1).leading_zeros() as usize; // 0 to 63
    let width = cmp::min((len + 1) * w, bits::WORD_BITS);
    let mask = bits::low_set(width) as usize;
    rng.gen::<usize>() & mask
}

// Generate `n` random values, with the widths (almost) geometrically distributed (p = 0.5) in blocks of `w` bits.
fn generate_values(n: usize, w: usize) -> Vec<usize> {
    let mut result = Vec::with_capacity(n);
    let mut rng = rand::thread_rng();
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
fn generate_runs(n: usize, sigma: usize, w: usize) -> Vec<(usize, usize)> {
    let sigma = if sigma == 0 { usize::MAX } else { sigma };
    let mut result = Vec::with_capacity(n);
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let c: usize = rng.gen_range(0, sigma);
        let len = generate_value(&mut rng, w) + 1;
        result.push((c, len));
    }
    result
}

fn encode_runs(encoder: &mut RLE, runs: &[(usize, usize)], name: &str) {
    assert_eq!(encoder.len(), 0, "[{}]: Newly created encoder contains runs", name);
    assert!(encoder.is_empty(), "[{}]: Newly created encoder is not empty", name);
    for (c, len) in runs.iter() {
        encoder.write(*c, *len);
    }
    assert!(encoder.len() >= runs.len(), "[{}]: The encoding is shorter than the number of runs", name);
    assert!(!encoder.is_empty(), "[{}]: The encoding is empty", name);
}

fn check_runs(encoder: &RLE, truth: &[(usize, usize)], name: &str) {
    let mut iter = RLEIter::new(encoder.as_ref(), encoder.sigma());
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
    let mut encoder = RLE::new(sigma);
    encode_runs(&mut encoder, &runs, name);
    check_runs(&encoder, &runs, name);
}

fn add_run(encoder: &mut RLE, truth: &mut Vec<(usize, usize)>, len: usize, bytes: usize, name: &str) {
    let old_len = encoder.len();
    encoder.write(encoder.sigma() - 1, len);
    truth.push((encoder.sigma() - 1, len));
    assert_eq!(encoder.len() - old_len, bytes, "[{}]: Run of length {} not encoded using {} byte(s)", name, len, bytes);
}

fn test_threshold(sigma: usize, name: &str) {
    let (sigma, threshold) = RLE::sanitize(sigma);
    let mut encoder = RLE::new(sigma);
    let mut truth: Vec<(usize, usize)> = Vec::new();
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
    let mut encoder = ByteCode::new();
    encoder.write(sigma);
    let mut prev = 0;
    for (node, offset) in edges.iter() {
        encoder.write(*node - prev); encoder.write(*offset);
        prev = *node;
    }
    let mut encoder = RLE::from_byte_code(encoder, sigma);
    for (c, len) in runs.iter() {
        encoder.write(*c, *len);
    }

    // Decompress the record.
    let mut iter = ByteCodeIter::new(encoder.as_ref());
    assert_eq!(iter.offset(), 0, "Newly created iterator is not at offset 0");
    assert_eq!(iter.next(), Some(sigma), "Invalid alphabet size in the record");
    let mut prev = 0;
    for i in 0..sigma {
        let node = iter.next().unwrap() + prev;
        assert_eq!(node, edges[i].0, "Invalid successor node {}", i);
        prev = node;
        assert_eq!(iter.next(), Some(edges[i].1), "Invalid record offset for edge {}", i);
    }
    let mut iter = RLEIter::from_byte_code(iter, sigma);
    let mut decoded: Vec<(usize, usize)> = Vec::new();
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
