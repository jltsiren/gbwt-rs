use super::*;

use crate::support;

use simple_sds::serialize;

//-----------------------------------------------------------------------------

#[test]
fn statistics() {
    let filename = support::get_test_data("example.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();

    assert_eq!(graph.nodes(), 12, "Invalid number of nodes");
    assert!(!graph.is_empty(), "The graph should not be empty");
    assert_eq!(graph.sequences(), 15, "Invalid number of sequences");
    assert_eq!(graph.segments(), 0, "The graph should not contain segments");
}

#[test]
fn sequences() {
    let filename = support::get_test_data("example.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();
    let truth: Vec<&str> = vec!["G", "A", "T", "T", "A", "C", "A", "", "", "", "G", "A", "T", "T", "A"];

    for i in 0..graph.sequences() {
        assert_eq!(graph.sequence(i), truth[i].as_bytes(), "Invalid sequence {}", i);
    }

    assert!(graph.iter().eq(truth.iter().map(|x| x.as_bytes())), "Invalid sequences from an iterator");
}

#[test]
fn serialize() {
    let filename = support::get_test_data("example.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();
    serialize::test(&graph, "graph", None, true);
}

#[test]
fn no_translation() {
    let filename = support::get_test_data("example.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();

    assert!(!graph.has_translation(), "The graph should not contain a node-to-segment translation");
    let mut iter = graph.segment_iter();
    assert_eq!(iter.next(), None, "Got a segment from an iterator");
}

//-----------------------------------------------------------------------------

#[test]
fn statistics_trans() {
    let filename = support::get_test_data("translation.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();

    assert_eq!(graph.nodes(), 9, "Invalid number of nodes");
    assert!(!graph.is_empty(), "The graph should not be empty");
    assert_eq!(graph.sequences(), 9, "Invalid number of sequences");
    assert_eq!(graph.segments(), 7, "Invalid number of segments");
}

#[test]
fn sequences_trans() {
    let filename = support::get_test_data("translation.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();
    let truth: Vec<&str> = vec!["GA", "T", "T", "A", "CA", "G", "A", "T", "TA"];

    for i in 0..graph.sequences() {
        assert_eq!(graph.sequence(i), truth[i].as_bytes(), "Invalid sequence {}", i);
    }

    assert_eq!(graph.iter().len(), graph.sequences(), "Invalid number of sequences from an iterator");
    assert!(graph.iter().eq(truth.iter().map(|x| x.as_bytes())), "Invalid sequences from an iterator");
}

#[test]
fn serialize_trans() {
    let filename = support::get_test_data("translation.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();
    serialize::test(&graph, "translation", None, true);
}

#[test]
fn translation() {
    let filename = support::get_test_data("translation.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();
    let truth: Vec<(&str, Range<usize>, &str)> = vec![
        ("s11", 1..3, "GAT"),
        ("s12", 3..4, "T"),
        ("s13", 4..5, "A"),
        ("s14", 5..7, "CAG"),
        ("s15", 7..8, "A"),
        ("s16", 8..9, "T"),
        ("s17", 9..10, "TA"),
    ];

    assert!(graph.has_translation(), "The graph does not contain a node-to-segment translation");

    for i in 0..graph.segments() {
        assert_eq!(graph.segment_name(i), truth[i].0.as_bytes(), "Invalid name for segment {}", i);
        assert_eq!(graph.segment_nodes(i), truth[i].1, "Invalid node range for segment {}", i);
        assert_eq!(graph.segment_sequence(i), truth[i].2.as_bytes(), "Invalid sequence for segment {}", i);
    }

    assert_eq!(graph.segment_iter().len(), graph.segments(), "Invalid number of segments from an iterator");
    assert!(graph.segment_iter().eq(truth.iter().map(|(n, r, s)| (n.as_bytes(), r.clone(), s.as_bytes()))), "Invalid segments from an iterator");
}

//-----------------------------------------------------------------------------
