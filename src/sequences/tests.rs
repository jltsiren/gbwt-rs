use super::*;

use crate::{GBWT, GBWTBuilder, GraphBuilder, Orientation};
use crate::support;

use simple_sds::serialize;

//-----------------------------------------------------------------------------

#[test]
fn statistics() {
    let filename = support::get_test_data("example.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();

    assert_eq!(sequences.nodes(), 12, "Invalid number of nodes");
    assert!(!sequences.is_empty(), "The graph should not be empty");
    assert_eq!(sequences.sequences(), 15, "Invalid number of sequences");
    assert_eq!(sequences.segments(), 0, "The graph should not contain segments");
}

#[test]
fn sequences() {
    let filename = support::get_test_data("example.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();
    let truth: Vec<&str> = vec!["G", "A", "T", "T", "A", "C", "A", "", "", "", "G", "A", "T", "T", "A"];

    for i in 0..sequences.sequences() {
        assert_eq!(sequences.sequence(i), truth[i].as_bytes(), "Invalid sequence {}", i);
        assert_eq!(sequences.sequence_len(i), truth[i].len(), "Invalid sequence length {}", i);
    }

    assert_eq!(sequences.iter().len(), sequences.sequences(), "Invalid number of sequences from an iterator");
    assert!(sequences.iter().eq(truth.iter().map(|x| x.as_bytes())), "Invalid sequences from an iterator");
}

#[test]
fn serialize() {
    let filename = support::get_test_data("example.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();
    serialize::test(&sequences, "sequences", None, true);
}

#[test]
fn serialize_versions() {
    let filename = support::get_test_data("example.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();
    let lossless_versions: Vec<usize> = (Sequences::MIN_VERSION..=Sequences::MAX_VERSION).collect();
    serialize::test_versions(&sequences, "sequences", &lossless_versions);
}

#[test]
fn no_translation() {
    let filename = support::get_test_data("example.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();

    assert!(!sequences.has_translation(), "The graph should not contain a node-to-segment translation");
    let mut iter = sequences.segment_iter();
    assert_eq!(iter.next(), None, "Got a segment from an iterator");
}

//-----------------------------------------------------------------------------

#[test]
fn statistics_trans() {
    let filename = support::get_test_data("translation.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();

    assert_eq!(sequences.nodes(), 9, "Invalid number of nodes");
    assert!(!sequences.is_empty(), "The graph should not be empty");
    assert_eq!(sequences.sequences(), 11, "Invalid number of sequences");
    assert_eq!(sequences.segments(), 8, "Invalid number of segments");
}

#[test]
fn sequences_trans() {
    let filename = support::get_test_data("translation.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();
    let truth: Vec<&str> = vec!["GA", "T", "T", "A", "CA", "G", "", "", "A", "T", "TA"];

    for i in 0..sequences.sequences() {
        assert_eq!(sequences.sequence(i), truth[i].as_bytes(), "Invalid sequence {}", i);
        assert_eq!(sequences.sequence_len(i), truth[i].len(), "Invalid sequence length {}", i);
    }

    assert_eq!(sequences.iter().len(), sequences.sequences(), "Invalid number of sequences from an iterator");
    assert!(sequences.iter().eq(truth.iter().map(|x| x.as_bytes())), "Invalid sequences from an iterator");
}

#[test]
fn serialize_trans() {
    let filename = support::get_test_data("translation.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();
    serialize::test(&sequences, "sequences-translation", None, true);
}

#[test]
fn serialize_trans_versions() {
    let filename = support::get_test_data("translation.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();
    let lossless_versions: Vec<usize> = (Sequences::MIN_VERSION..=Sequences::MAX_VERSION).collect();
    serialize::test_versions(&sequences, "sequences-translation", &lossless_versions);
}

#[test]
fn translation() {
    let filename = support::get_test_data("translation.gg");
    let sequences: Sequences = serialize::load_from(&filename).unwrap();
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

    assert!(sequences.has_translation(), "The graph does not contain a node-to-segment translation");

    for i in 0..sequences.segments() {
        assert_eq!(sequences.segment(i), truth[i], "Invalid segment {}", i);
        for j in truth[i].nodes.clone() {
            assert_eq!(sequences.node_to_segment(j), truth[i], "Invalid segment for node {}", j);
        }
        assert_eq!(sequences.segment_name(i), truth[i].name, "Invalid name for segment {}", i);
        assert_eq!(sequences.segment_nodes(i), truth[i].nodes, "Invalid node range for segment {}", i);
        assert_eq!(sequences.segment_sequence(i), truth[i].sequence, "Invalid sequence for segment {}", i);
        assert_eq!(sequences.segment_len(i), truth[i].sequence.len(), "Invalid sequence length for segment {}", i);
    }

    // Iterate forward.
    assert_eq!(sequences.segment_iter().len(), sequences.segments(), "Invalid number of segments from an iterator");
    assert!(sequences.segment_iter().eq(truth.iter().cloned()), "Invalid segments from an iterator (forward)");

    // Iterate backward.
    assert!(sequences.segment_iter().rev().eq(truth.iter().rev().cloned()), "Invalid segments from an iterator (backward)");
}

//-----------------------------------------------------------------------------

#[test]
fn from_graph_no_translation() {
    let filename = support::get_test_data("example.gbwt");
    let gbwt: GBWT = serialize::load_from(&filename).unwrap();

    let mut graph = GraphBuilder::new();
    graph.add_node(11, b"G").unwrap();
    graph.add_node(12, b"A").unwrap();
    graph.add_node(13, b"T").unwrap();
    graph.add_node(14, b"T").unwrap();
    graph.add_node(15, b"A").unwrap();
    graph.add_node(16, b"C").unwrap();
    graph.add_node(17, b"A").unwrap();
    graph.add_node(21, b"G").unwrap();
    graph.add_node(22, b"A").unwrap();
    graph.add_node(23, b"T").unwrap();
    graph.add_node(24, b"T").unwrap();
    graph.add_node(25, b"A").unwrap();

    let sequences = Sequences::from_graph(&gbwt, &graph).unwrap();

    let truth_filename = support::get_test_data("example.gg");
    let truth: Sequences = serialize::load_from(&truth_filename).unwrap();

    // The empty mapping may use a different internal integer width than the
    // loaded reference, so we compare the meaningful parts instead of the
    // whole structure.
    assert_eq!(sequences.nodes(), truth.nodes(), "Invalid number of nodes");
    assert_eq!(sequences.sequences(), truth.sequences(), "Invalid number of sequences");
    assert!(!sequences.has_translation(), "The graph should not contain a translation");
    assert_eq!(sequences.segments(), 0, "The graph should not contain segments");
    for i in 0..sequences.sequences() {
        assert_eq!(sequences.sequence(i), truth.sequence(i), "Invalid sequence {}", i);
    }
}

#[test]
fn from_graph_translation() {
    let filename = support::get_test_data("translation.gbwt");
    let gbwt: GBWT = serialize::load_from(&filename).unwrap();

    let mut graph = GraphBuilder::new();
    graph.translate_segment("s11", b"GAT", 2).unwrap();
    graph.translate_segment("s12", b"T", 2).unwrap();
    graph.translate_segment("s13", b"A", 2).unwrap();
    graph.translate_segment("s14", b"CAG", 2).unwrap();
    graph.translate_segment("unused", b"GAA", 2).unwrap();
    graph.translate_segment("s15", b"A", 2).unwrap();
    graph.translate_segment("s16", b"T", 2).unwrap();
    graph.translate_segment("s17", b"TA", 2).unwrap();

    let sequences = Sequences::from_graph(&gbwt, &graph).unwrap();

    let truth_filename = support::get_test_data("translation.gg");
    let truth: Sequences = serialize::load_from(&truth_filename).unwrap();
    assert_eq!(sequences, truth, "Sequences built from the graph do not match the reference");
}

//-----------------------------------------------------------------------------

// Builds a bidirectional GBWT with a single path over the given nodes.
fn single_path_gbwt(nodes: &[usize]) -> GBWT {
    let path: Vec<usize> = nodes.iter().map(|&id| support::encode_node(id, Orientation::Forward)).collect();
    let mut builder = GBWTBuilder::new(true, false, 100);
    builder.insert(&path, None).unwrap();
    builder.build().unwrap()
}

#[test]
fn from_graph_node_offset() {
    // Without a translation, the node ids do not have to start from 1.
    let gbwt = single_path_gbwt(&[5, 6]);
    let mut graph = GraphBuilder::new();
    graph.add_node(5, b"A").unwrap();
    graph.add_node(6, b"CG").unwrap();
    let sequences = Sequences::from_graph(&gbwt, &graph).unwrap();
    assert_eq!(sequences.nodes(), 2, "Invalid number of nodes");
    assert_eq!(sequences.sequences(), 2, "Invalid number of sequences");
    assert_eq!(sequences.sequence(0), b"A", "Invalid sequence 0");
    assert_eq!(sequences.sequence(1), b"CG", "Invalid sequence 1");
}

#[test]
fn from_graph_translation_offset() {
    // The path only visits the last segment, so the first node of the GBWT is not 1.
    let gbwt = single_path_gbwt(&[3]);
    let mut graph = GraphBuilder::new();
    graph.translate_segment("a", b"A", 0).unwrap();
    graph.translate_segment("b", b"C", 0).unwrap();
    graph.translate_segment("c", b"G", 0).unwrap();
    assert!(Sequences::from_graph(&gbwt, &graph).is_err(), "Built a translation that does not start from node 1");
}

#[test]
fn from_graph_missing_node() {
    let gbwt = single_path_gbwt(&[1, 2]);
    let mut graph = GraphBuilder::new();
    graph.add_node(1, b"A").unwrap();
    assert!(Sequences::from_graph(&gbwt, &graph).is_err(), "Built sequences with a missing node");
}

#[test]
fn from_graph_empty_gbwt() {
    let gbwt = GBWTBuilder::new(true, false, 100).build().unwrap();
    let graph = GraphBuilder::new();
    let sequences = Sequences::from_graph(&gbwt, &graph).unwrap();
    assert!(sequences.is_empty(), "Sequences built from an empty GBWT are not empty");
    assert_eq!(sequences.nodes(), 0, "Invalid number of nodes");
    assert!(!sequences.has_translation(), "Unexpected translation");
}

//-----------------------------------------------------------------------------
