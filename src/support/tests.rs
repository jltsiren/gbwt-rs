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
    assert_eq!(array.size_in_bytes(), len, "Invalid size estimate for the serialized array");

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
