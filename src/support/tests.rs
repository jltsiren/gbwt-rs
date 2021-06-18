use super::*;

use simple_sds::serialize;

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
