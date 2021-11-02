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
        assert_eq!(graph.sequence_len(i), truth[i].len(), "Invalid sequence length {}", i);
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
    assert_eq!(graph.sequences(), 11, "Invalid number of sequences");
    assert_eq!(graph.segments(), 8, "Invalid number of segments");
}

#[test]
fn sequences_trans() {
    let filename = support::get_test_data("translation.gg");
    let graph: Graph = serialize::load_from(&filename).unwrap();
    let truth: Vec<&str> = vec!["GA", "T", "T", "A", "CA", "G", "", "", "A", "T", "TA"];

    for i in 0..graph.sequences() {
        assert_eq!(graph.sequence(i), truth[i].as_bytes(), "Invalid sequence {}", i);
        assert_eq!(graph.sequence_len(i), truth[i].len(), "Invalid sequence length {}", i);
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
    let truth: Vec<Segment> = vec![
        Segment::from_fields(0, "s11".as_bytes(), 1..3, "GAT".as_bytes()),
        Segment::from_fields(1, "s12".as_bytes(), 3..4, "T".as_bytes()),
        Segment::from_fields(2, "s13".as_bytes(), 4..5, "A".as_bytes()),
        Segment::from_fields(3, "s14".as_bytes(), 5..7, "CAG".as_bytes()),
        Segment::from_fields(4, "".as_bytes(), 7..9, "".as_bytes()),
        Segment::from_fields(5, "s15".as_bytes(), 9..10, "A".as_bytes()),
        Segment::from_fields(6, "s16".as_bytes(), 10..11, "T".as_bytes()),
        Segment::from_fields(7, "s17".as_bytes(), 11..12, "TA".as_bytes()),
    ];

    assert!(graph.has_translation(), "The graph does not contain a node-to-segment translation");

    for i in 0..graph.segments() {
        assert_eq!(graph.segment(i), truth[i], "Invalid segment {}", i);
        for j in truth[i].nodes.clone() {
            assert_eq!(graph.node_to_segment(j), truth[i], "Invalid segment for node {}", j);
        }
        assert_eq!(graph.segment_name(i), truth[i].name, "Invalid name for segment {}", i);
        assert_eq!(graph.segment_nodes(i), truth[i].nodes, "Invalid node range for segment {}", i);
        assert_eq!(graph.segment_sequence(i), truth[i].sequence, "Invalid sequence for segment {}", i);
        assert_eq!(graph.segment_len(i), truth[i].sequence.len(), "Invalid sequence length for segment {}", i);
    }

    assert_eq!(graph.segment_iter().len(), graph.segments(), "Invalid number of segments from an iterator");
    assert!(graph.segment_iter().eq(truth.iter().cloned()), "Invalid segments from an iterator");
}

//-----------------------------------------------------------------------------
