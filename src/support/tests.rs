use super::*;

use simple_sds::serialize;

use std::fs;

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

fn serialize_array(array: &StringArray) {
    let filename = serialize::temp_file_name("string-array");
    serialize::serialize_to(array, &filename).unwrap();

    let metadata = fs::metadata(&filename).unwrap();
    let len = metadata.len() as usize;
    assert_eq!(array.size_in_bytes(), len, "Invalid size estimate for the serialized StringArray");

    let copy: StringArray = serialize::load_from(&filename).unwrap();
    assert_eq!(copy, *array, "Serialization changed the StringArray");

    fs::remove_file(&filename).unwrap();
}

#[test]
fn empty_string_array() {
    let truth: Vec<&str> = Vec::new();
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    serialize_array(&array);
}

#[test]
fn non_empty_string_array() {
    let truth = vec!["first", "second", "third", "fourth"];
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    serialize_array(&array);
}

#[test]
fn array_with_empty_strings() {
    // Serialization with an empty string at the end used to fail in the original GBWT implementation.
    let truth = vec!["first", "second", "", "fourth", ""];
    let array = StringArray::from(truth.as_slice());
    check_array(&array, &truth);
    serialize_array(&array);
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
    for i in 0..missing.len() {
        assert_eq!(dict.id(missing[i]), None, "String {} should not be present", missing[i]);
    }
}

fn serialize_dict(dict: &Dictionary) {
    let filename = serialize::temp_file_name("dictionary");
    serialize::serialize_to(dict, &filename).unwrap();

    let metadata = fs::metadata(&filename).unwrap();
    let len = metadata.len() as usize;
    assert_eq!(dict.size_in_bytes(), len, "Invalid size estimate for the serialized Dictionary");

    let copy: Dictionary = serialize::load_from(&filename).unwrap();
    assert_eq!(copy, *dict, "Serialization changed the Dictionary");

    fs::remove_file(&filename).unwrap();
}

#[test]
fn empty_dict() {
    let truth: Vec<&str> = Vec::new();
    let missing = vec!["this", "should", "not", "exist"];
    let dict = Dictionary::try_from(truth.as_slice()).unwrap();
    check_dict(&dict, &truth, &missing);
    serialize_dict(&dict);
}

#[test]
fn non_empty_dict() {
    let truth = vec!["first", "second", "third", "fourth"];
    let missing = vec!["this", "should", "not", "exist"];
    let dict = Dictionary::try_from(truth.as_slice()).unwrap();
    check_dict(&dict, &truth, &missing);
    serialize_dict(&dict);
}

#[test]
fn dict_from_duplicates() {
    let source = vec!["this", "contains", "contains", "many", "duplicates", "many", "this"];
    let result = Dictionary::try_from(source);
    assert!(result.is_err(), "Did not get an error from a source with duplicate strings");
}

//-----------------------------------------------------------------------------
