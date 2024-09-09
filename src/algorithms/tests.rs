use super::*;

use crate::Orientation;
use rand::Rng;
use simple_sds::serialize;

//use simple_sds::serialize;

//-----------------------------------------------------------------------------

#[test]
fn empty_lcs() {
    let a: Vec<usize> = Vec::new();
    let b: Vec<usize> = Vec::new();
    assert!(lcs(&a, &b).is_empty(), "Non-empty LCS for empty inputs");

    let non_empty: Vec<usize> = vec![1, 2, 3];
    assert!(lcs(&a, &non_empty).is_empty(), "Non-empty LCS for empty first input");
    assert!(lcs(&non_empty, &b).is_empty(), "Non-empty LCS for empty second input");
}

#[test]
fn first_last_lcs() {
    let a: Vec<usize> = vec![1, 2, 3, 4, 5];

    let no_no: Vec<usize> = vec![6, 2, 4, 7];
    let no_no_left: Vec<(usize, usize)> = vec![(1, 1), (2, 3)];
    assert_eq!(lcs(&no_no, &a), no_no_left, "Incorrect LCS for no first / no last / left");
    let no_no_right: Vec<(usize, usize)> = vec![(1, 1), (3, 2)];
    assert_eq!(lcs(&a, &no_no), no_no_right, "Incorrect LCS for no first / no last / right");

    let no_yes: Vec<usize> = vec![6, 2, 4, 5];
    let no_yes_left: Vec<(usize, usize)> = vec![(1, 1), (2, 3), (3, 4)];
    assert_eq!(lcs(&no_yes, &a), no_yes_left, "Incorrect LCS for no first / yes last / left");
    let no_yes_right: Vec<(usize, usize)> = vec![(1, 1), (3, 2), (4, 3)];
    assert_eq!(lcs(&a, &no_yes), no_yes_right, "Incorrect LCS for no first / yes last / right");

    let yes_no: Vec<usize> = vec![1, 2, 4, 7];
    let yes_no_left: Vec<(usize, usize)> = vec![(0, 0), (1, 1), (2, 3)];
    assert_eq!(lcs(&yes_no, &a), yes_no_left, "Incorrect LCS for yes first / no last / left");
    let yes_no_right: Vec<(usize, usize)> = vec![(0, 0), (1, 1), (3, 2)];
    assert_eq!(lcs(&a, &yes_no), yes_no_right, "Incorrect LCS for yes first / no last / right");

    let yes_yes: Vec<usize> = vec![1, 2, 4, 5];
    let yes_yes_left: Vec<(usize, usize)> = vec![(0, 0), (1, 1), (2, 3), (3, 4)];
    assert_eq!(lcs(&yes_yes, &a), yes_yes_left, "Incorrect LCS for yes first / yes last / left");
    let yes_yes_right: Vec<(usize, usize)> = vec![(0, 0), (1, 1), (3, 2), (4, 3)];
    assert_eq!(lcs(&a, &yes_yes), yes_yes_right, "Incorrect LCS for yes first / yes last / right");
}

#[test]
fn minimal_maximal_lcs() {
    let four_odd: Vec<usize> = vec![1, 3, 5, 7];
    let five_odd: Vec<usize> = vec![1, 3, 5, 7, 9];
    let four_even: Vec<usize> = vec![4, 6, 8, 10];
    let five_even: Vec<usize> = vec![2, 4, 6, 8, 10];

    assert!(lcs(&four_odd, &four_even).is_empty(), "Non-empty LCS for disjoint inputs");

    let four_four: Vec<(usize, usize)> = vec![(0, 0), (1, 1), (2, 2), (3, 3)];
    assert_eq!(lcs(&four_odd, &four_odd), four_four, "Incorrect LCS for equal inputs");

    let prefix: Vec<(usize, usize)> = vec![(0, 0), (1, 1), (2, 2), (3, 3)];
    assert_eq!(lcs(&four_odd, &five_odd), prefix, "Incorrect LCS for prefix / full");
    assert_eq!(lcs(&five_odd, &four_odd), prefix, "Incorrect LCS for full / prefix");

    let left_suffix: Vec<(usize, usize)> = vec![(0, 1), (1, 2), (2, 3), (3, 4)];
    assert_eq!(lcs(&four_even, &five_even), left_suffix, "Incorrect LCS for suffix / full");
    let right_suffix: Vec<(usize, usize)> = vec![(1, 0), (2, 1), (3, 2), (4, 3)];
    assert_eq!(lcs(&five_even, &four_even), right_suffix, "Incorrect LCS for full / suffix");
}

fn random_sequence(len: usize, sigma: usize) -> Vec<usize> {
    let mut rng = rand::thread_rng();
    (0..len).map(|_| rng.gen_range(0..sigma)).collect()
}

fn random_lcs_instance(len: usize, sigma: usize) {
    let left = random_sequence(len, sigma);
    let right = random_sequence(len, sigma);
    let result = lcs(&left, &right);

    for i in 0..result.len() {
        assert!(result[i].0 < left.len(), "Invalid left LCS position {} (len {}, sigma {})", i, len, sigma);
        assert!(result[i].1 < right.len(), "Invalid right LCS position {} (len {}, sigma {})", i, len, sigma);
        if i > 0 {
            assert!(result[i].0 > result[i - 1].0, "Non-increasing left LCS position {} (len {}, sigma {})", i, len, sigma);
            assert!(result[i].1 > result[i - 1].1, "Non-increasing right LCS position {} (len {}, sigma {})", i, len, sigma);
        }
        assert_eq!(left[result[i].0], right[result[i].1], "Incorrect LCS element {} (len {}, sigma {})", i, len, sigma);
    }
}

#[test]
fn random_lcs() {
    random_lcs_instance(10, 5);
    random_lcs_instance(100, 10);
    random_lcs_instance(300, 20);
}

// TODO: Large almost equal sequences

//-----------------------------------------------------------------------------

fn extract_path(graph: &GBZ, path_id: usize, orientation: Orientation) -> Vec<usize> {
    graph.path(path_id, orientation)
        .unwrap()
        .map(|(id, o)| support::encode_node(id, o))
        .collect()
}

fn node_len(graph: &GBZ, handle: usize) -> usize {
    graph.sequence_len(support::node_id(handle)).unwrap_or(0)
}

fn path_len(graph: &GBZ, path: &[usize]) -> usize {
    path.iter().map(|handle| node_len(graph, *handle)).sum()
}

#[test]
fn empty_path_lcs() {
    let filename = support::get_test_data("translation.gbz");
    let graph: GBZ = serialize::load_from(&filename).unwrap();
    let empty: Vec<usize> = Vec::new();
    let path = extract_path(&graph, 0, Orientation::Forward);

    let truth = (Vec::new(), 0);
    assert_eq!(path_lcs(&empty, &empty, &graph), truth, "Non-empty LCS for empty inputs");
    assert_eq!(path_lcs(&empty, &path, &graph), truth, "Non-empty LCS for empty first input");
    assert_eq!(path_lcs(&path, &empty, &graph), truth, "Non-empty LCS for empty second input");

    let reverse = extract_path(&graph, 0, Orientation::Reverse);
    assert_eq!(path_lcs(&path, &reverse, &graph), truth, "Non-empty LCS non-overlapping paths");
}

#[test]
fn identical_path_lcs() {
    let filename = support::get_test_data("translation.gbz");
    let graph: GBZ = serialize::load_from(&filename).unwrap();

    for path_id in 0..graph.paths() {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let path = extract_path(&graph, path_id, orientation);
            let truth = path.iter().enumerate().map(|(i, _)| (i, i)).collect();
            let len = path_len(&graph, &path);
            assert_eq!(path_lcs(&path, &path, &graph), (truth, len), "Incorrect LCS for identical paths {} {}", path_id, orientation);
        }
    }
}

fn path_lcs_instance(graph: &GBZ, left_id: usize, right_id: usize, orientation: Orientation) {
    let left = extract_path(graph, left_id, orientation);
    let right = extract_path(graph, right_id, orientation);
    let (result, len) = path_lcs(&left, &right, graph);

    let mut total_len = 0;
    let test = format!("(paths {} and {} {})", left_id, right_id, orientation);
    for i in 0..result.len() {
        assert!(result[i].0 < left.len(), "Invalid left LCS position {} {}", i, test);
        assert!(result[i].1 < right.len(), "Invalid right LCS position {} {}", i, test);
        if i > 0 {
            assert!(result[i].0 > result[i - 1].0, "Non-increasing left LCS position {} {}", i, test);
            assert!(result[i].1 > result[i - 1].1, "Non-increasing right LCS position {} {}", i, test);
        }
        assert_eq!(left[result[i].0], right[result[i].1], "Incorrect LCS element {} {}", i, test);
        total_len += node_len(graph, left[result[i].0]);
    }
    assert_eq!(total_len, len, "Incorrect LCS length {}", test);
}

#[test]
fn real_path_lcs() {
    let filename = support::get_test_data("translation.gbz");
    let graph: GBZ = serialize::load_from(&filename).unwrap();

    for left_id in 0..graph.paths() {
        for right_id in 0..graph.paths() {
            for orientation in [Orientation::Forward, Orientation::Reverse] {
                path_lcs_instance(&graph, left_id, right_id, orientation);
            }
        }
    }
}

#[test]
fn path_lcs_weights() {
    let filename = support::get_test_data("translation.gbz");
    let graph: GBZ = serialize::load_from(&filename).unwrap();

    // This is an artificial scenario where the weights matter.
    // 1, 5, and 11 are weight 2 nodes and form the path LCS.
    // Normal LCS would return 1, 2, 3, 4.
    let left: Vec<usize> = vec![
        support::encode_node(1, Orientation::Forward),
        support::encode_node(2, Orientation::Forward),
        support::encode_node(3, Orientation::Forward),
        support::encode_node(4, Orientation::Forward),
        support::encode_node(5, Orientation::Forward),
        support::encode_node(11, Orientation::Forward),
    ];
    let right: Vec<usize> = vec![
        support::encode_node(1, Orientation::Forward),
        support::encode_node(5, Orientation::Forward),
        support::encode_node(11, Orientation::Forward),
        support::encode_node(2, Orientation::Forward),
        support::encode_node(3, Orientation::Forward),
        support::encode_node(4, Orientation::Forward),
    ];

    let unweighted_truth: Vec<(usize, usize)> = vec![(0, 0), (1, 3), (2, 4), (3, 5)];
    let path_truth: Vec<(usize, usize)> = vec![(0, 0), (4, 1), (5, 2)];
    let path_len = 6;

    assert_eq!(lcs(&left, &right), unweighted_truth, "Incorrect unweighted LCS");
    assert_eq!(path_lcs(&left, &right, &graph), (path_truth, path_len), "Incorrect path LCS");
}

//-----------------------------------------------------------------------------
