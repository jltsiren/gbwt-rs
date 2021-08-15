use super::*;

use simple_sds::serialize;

use std::collections::HashSet;

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

#[test]
fn sequence() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();

    for i in 0..index.sequences() {
        let extracted = extract_sequence(&index, i);
        let iterated: Vec<usize> = index.sequence(i).collect();
        assert_eq!(iterated, extracted, "Invalid sequence {} from an iterator", i);
    }
}

//-----------------------------------------------------------------------------

fn true_nodes() -> HashSet<usize> {
    let nodes: Vec<usize> = vec![11, 12, 13, 14, 15, 16, 17, 21, 22, 23, 24, 25];
    let mut result: HashSet<usize> = HashSet::new();
    for node in nodes.iter() {
        result.insert(support::encode_node(*node, false));
        result.insert(support::encode_node(*node, true));
    }
    result
}

fn count_occurrences(paths: &[Vec<usize>], subpath: &[usize]) -> usize {
    let mut result = 0;
    let reverse = support::reverse_path(subpath);
    for path in paths {
        for i in 0..path.len() {
            if path[i..].starts_with(subpath) {
                result += 1;
            }
            if path[..i + 1].ends_with(&reverse) {
                result += 1;
            }
        }
    }
    result
}

#[test]
fn find() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();

    for i in 0..index.alphabet_size() + 1 {
        if let Some(state) = index.find(i) {
            assert!(nodes.contains(&i), "Found a search state for a nonexistent node {}", i);
            assert_eq!(state.node, i, "Found an invalid search state for node {}", i);
            assert!(!state.is_empty(), "Found an empty search state for node {}", i);
        } else {
            assert!(!nodes.contains(&i), "Did not find a search state for node {}", i);
        }
    }
}

#[test]
fn extend() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();
    let paths = true_paths();

    // Check all possible and impossible extensions of the initial node.
    for &first in nodes.iter() {
        let start = index.find(first).unwrap();
        for i in 0..index.alphabet_size() + 1 {
            let count = count_occurrences(&paths, &[first, i]);
            if let Some(state) = index.extend(&start, i) {
                assert_eq!(state.len(), count, "Invalid number of occurrences for substring {} to {}", first, i);
            } else {
                assert_eq!(count, 0, "Could not find the occurrences of substring {} to {}", first, i);
            }
        }
    }

    // Search for all existing subpaths.
    for i in 0..paths.len() {
        let path = &paths[i];
        for j in 0..path.len() {
            let mut forward = index.find(path[j]).unwrap();
            for k in j + 1..path.len() {
                if let Some(state) = index.extend(&forward, path[k]) {
                    let count = count_occurrences(&paths, &path[j..k + 1]);
                    assert_eq!(state.len(), count, "Invalid number of occurrences for path {} at {}..{}", i, j, k + 1);
                    forward = state;
                } else {
                    panic!("Could not find occurrences of path {} at {}..{}", i, j, k + 1);
                }
            }

            let mut backward = index.find(support::flip_node(path[j])).unwrap();
            for k in (0..j).rev() {
                if let Some(state) = index.extend(&backward, support::flip_node(path[k])) {
                    let count = count_occurrences(&paths, &path[k..j + 1]); // No need to reverse the pattern here.
                    assert_eq!(state.len(), count, "Invalid number of occurrences for path {} at {}..{} (reversed)", i, k, j + 1);
                    backward = state;
                } else {
                    panic!("Could not find occurrences of path {} at {}..{} (reversed)", i, k, j + 1);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

fn bd_search(index: &GBWT, path: &[usize], first: usize, range: &Range<usize>) -> Option<BidirectionalState> {
    let mut state = index.bd_find(path[first])?;
    for i in first + 1..range.end {
        state = index.extend_forward(&state, path[i])?;
    }
    for i in (range.start..first).rev() {
        state = index.extend_backward(&state, path[i])?;
    }
    Some(state)
}

#[test]
fn bd_find() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();

    for i in 0..index.alphabet_size() + 1 {
        if let Some(state) = index.bd_find(i) {
            assert!(nodes.contains(&i), "Found a search state for a nonexistent node {}", i);
            assert_eq!(state.forward.node, i, "Found an invalid search state for node {}", i);
            assert!(!state.is_empty(), "Found an empty search state for node {}", i);
            assert_eq!(state.reverse.node, support::flip_node(i), "Found an invalid reverse node for node {}", i);
            assert_eq!(state.reverse.len(), state.forward.len(), "Invalid reverse range length for node {}", i);
        } else {
            assert!(!nodes.contains(&i), "Did not find a search state for node {}", i);
        }
    }
}

#[test]
fn bd_extend() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();
    let paths = true_paths();

    // Check all possible and impossible extensions of the initial node.
    for &first in nodes.iter() {
        let start = index.bd_find(first).unwrap();
        for i in 0..index.alphabet_size() + 1 {
            // Forward.
            let count = count_occurrences(&paths, &[first, i]);
            if let Some(state) = index.extend_forward(&start, i) {
                assert_eq!(state.len(), count, "Invalid number of occurrences for substring {} to {} (forward)", first, i);
            } else {
                assert_eq!(count, 0, "Could not find the occurrences of substring {} to {} (forward)", first, i);
            }
            // Backward.
            let count = count_occurrences(&paths, &[i, first]);
            if let Some(state) = index.extend_backward(&start, i) {
                assert_eq!(state.len(), count, "Invalid number of occurrences for substring {} to {} (backward)", i, first);
            } else {
                assert_eq!(count, 0, "Could not find the occurrences of substring {} to {} (backward)", i, first);
            }
        }
    }

    // Search for all existing subpaths.
    for i in 0..paths.len() {
        let path = &paths[i];
        for first in 0..path.len() {
            for start in 0..first + 1 {
                for end in first + 1..path.len() + 1 {
                    // Forward path.
                    let count = count_occurrences(&paths, &path[start..end]);
                    if let Some(state) = bd_search(&index, &path, first, &(start..end)) {
                        assert_eq!(state.len(), count, "Invalid number of occurrences for path {} at {}..{} from {}", i, start, end, first);
                        assert_eq!(state.reverse.len(), state.len(), "Invalid reverse state length for path {} at {}..{} from {}", i, start, end, first);
                        assert_eq!(state.forward.node, path[end - 1], "Invalid final node for path {} at {}..{} from {}", i, start, end, first);
                        assert_eq!(state.reverse.node, support::flip_node(path[start]), "Invalid initial node for path {} at {}..{} from {}", i, start, end, first);
                    } else {
                        panic!("Could not find occurrences of path {} at {}..{} from {}", i, start, end, first);
                    }

                    // Reverse path.
                    let reversed = support::reverse_path(&path);
                    let count = count_occurrences(&paths, &reversed[start..end]);
                    if let Some(state) = bd_search(&index, &reversed, first, &(start..end)) {
                        assert_eq!(state.len(), count, "Invalid number of occurrences for path {} at {}..{} from {} (reversed)", i, start, end, first);
                        assert_eq!(state.reverse.len(), state.len(), "Invalid reverse state length for path {} at {}..{} from {} (reversed)", i, start, end, first);
                        assert_eq!(state.forward.node, reversed[end - 1], "Invalid final node for path {} at {}..{} from {} (reversed)", i, start, end, first);
                        assert_eq!(state.reverse.node, support::flip_node(reversed[start]), "Invalid initial node for path {} at {}..{} from {} (reversed)", i, start, end, first);
                    } else {
                        panic!("Could not find occurrences of path {} at {}..{} from {} (reversed)", i, start, end, first);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
