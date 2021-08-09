use super::*;

use simple_sds::serialize;

//-----------------------------------------------------------------------------

#[test]
fn statistics() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();

    assert_eq!(index.len(), 68, "Invalid total length");
    assert!(!index.is_empty(), "Invalid emptiness");
    assert_eq!(index.sequences(), 12, "Invalid number of sequences");
    assert_eq!(index.alphabet_size(), 52, "Invalid alphabet size");
    assert_eq!(index.alphabet_offset(), 21, "Invalid alphabet offset");
    assert_eq!(index.effective_size(), 31, "Invalid effective alphabet size");
    assert_eq!(index.first_node(), 22, "Invalid first node id");
    assert!(index.is_bidirectional(), "Index is not bidirectional");

    for i in 0..index.first_node() {
        assert!(!index.has_node(i), "Index should not contain node {}", i);
    }
    for i in index.first_node()..index.alphabet_size() {
        assert!(index.has_node(i), "Index should contain node {}", i);
    }
    assert!(!index.has_node(index.alphabet_size()), "Index contains a node past the end");
}

#[test]
fn serialize() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    serialize::test(&index, "gbwt", None, true);
}

//-----------------------------------------------------------------------------

fn extract_sequence(index: &GBWT, id: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut pos = index.start(id);
    while pos != None {
        result.push(pos.unwrap().0);
        pos = index.forward(pos.unwrap());
    }
    result
}

fn extract_backward(index: &GBWT, id: usize) -> Vec<usize> {
    let mut last = None;
    let mut pos = index.start(id);
    while pos != None {
        last = pos;
        pos = index.forward(pos.unwrap());
    }

    let mut result = Vec::new();
    pos = last;
    while pos != None {
        result.push(pos.unwrap().0);
        pos = index.backward(pos.unwrap());
    }

    result
}

fn true_paths() -> Vec<Vec<usize>> {
    vec![
        vec![support::encode_node(11, false), support::encode_node(12, false), support::encode_node(14, false), support::encode_node(15, false), support::encode_node(17, false)],
        vec![support::encode_node(21, false), support::encode_node(22, false), support::encode_node(24, false), support::encode_node(25, false)],
        vec![support::encode_node(11, false), support::encode_node(12, false), support::encode_node(14, false), support::encode_node(15, false), support::encode_node(17, false)],
        vec![support::encode_node(11, false), support::encode_node(13, false), support::encode_node(14, false), support::encode_node(16, false), support::encode_node(17, false)],
        vec![support::encode_node(21, false), support::encode_node(22, false), support::encode_node(24, false), support::encode_node(23, true), support::encode_node(21, true)],
        vec![support::encode_node(21, false), support::encode_node(22, false), support::encode_node(24, false), support::encode_node(25, false)],
    ]
}

#[test]
fn extract() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let truth = true_paths();

    for i in 0..index.sequences() / 2 {
        let forward = extract_sequence(&index, support::encode_node(i, false));
        assert_eq!(forward, truth[i], "Invalid forward path {}", i);
        let reverse = extract_sequence(&index, support::encode_node(i, true));
        assert_eq!(reverse.len(), forward.len(), "Invalid reverse path {} length", i);
        for j in 0..reverse.len() {
            let expected = support::flip_node(forward[forward.len() - j - 1]);
            assert_eq!(reverse[j], expected, "Invalid node {} on reverse path {}", j, i);
        }
    }
}

#[test]
fn backward() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();

    for i in 0..index.sequences() {
        let forward = extract_sequence(&index, i);
        let reverse = extract_backward(&index, i);
        assert_eq!(reverse.len(), forward.len(), "Invalid reverse sequence {} length", i);
        for j in 0..reverse.len() {
            let expected = forward[forward.len() - j - 1];
            assert_eq!(reverse[j], expected, "Invalid node {} on reverse sequence {}", j, i);
        }
    }
}


// Iter

//-----------------------------------------------------------------------------
