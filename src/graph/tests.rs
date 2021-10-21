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

    assert!(!graph.has_translation(), "The graph should not contain a node-to-segment translation");
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

//-----------------------------------------------------------------------------

// FIXME translation

//-----------------------------------------------------------------------------
