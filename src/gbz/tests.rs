use super::*;

use crate::support;

use simple_sds::serialize;

use std::collections::{BTreeSet, BTreeMap};

//-----------------------------------------------------------------------------

type Neighbors = BTreeMap<(usize, Orientation), BTreeSet<(usize, Orientation)>>;

fn get_pred_succ(gbz: &GBZ, edges: &Vec<(usize, Orientation, usize, Orientation)>) -> (Neighbors, Neighbors) {
    let mut predecessors: BTreeMap<(usize, Orientation), BTreeSet<(usize, Orientation)>> = BTreeMap::new();
    let mut successors: BTreeMap<(usize, Orientation), BTreeSet<(usize, Orientation)>> = BTreeMap::new();
    for node_id in gbz.node_iter() {
        predecessors.insert((node_id, Orientation::Forward), BTreeSet::new());
        predecessors.insert((node_id, Orientation::Reverse), BTreeSet::new());
        successors.insert((node_id, Orientation::Forward), BTreeSet::new());
        successors.insert((node_id, Orientation::Reverse), BTreeSet::new());
    }

    for (from, from_o, to, to_o) in edges.iter() {
        if let Some(pred) = predecessors.get_mut(&(*to, *to_o)) {
            pred.insert((*from, *from_o));
        }
        if let Some(pred) = predecessors.get_mut(&(*from, from_o.flip())) {
            pred.insert((*to, to_o.flip()));
        }
        if let Some(succ) = successors.get_mut(&(*from, *from_o)) {
            succ.insert((*to, *to_o));
        }
        if let Some(succ) = successors.get_mut(&(*to, to_o.flip())) {
            succ.insert((*from, from_o.flip()));
        }
    }

    (predecessors, successors)
}

fn check_pred_succ(gbz: &GBZ, predecessors: &Neighbors, successors: &Neighbors, node_id: usize, orientation: Orientation) {
    let name = if orientation == Orientation::Forward { "(forward)" } else { "(reverse)" };
    if gbz.has_node(node_id) {
        let truth = predecessors.get(&(node_id, orientation)).unwrap();
        let iter = gbz.predecessors(node_id, orientation);
        assert!(iter.is_some(), "Could not find predecessors for node {} {}", node_id, name);
        assert!(iter.unwrap().eq(truth.iter().cloned()), "Invalid predecessors for node {} {}", node_id, name);
        let truth = successors.get(&(node_id, orientation)).unwrap();
        let iter = gbz.successors(node_id, orientation);
        assert!(iter.is_some(), "Could not find successors for node {} {}", node_id, name);
        assert!(iter.unwrap().eq(truth.iter().cloned()), "Invalid successors for node {} {}", node_id, name);
    } else {
        assert!(gbz.predecessors(node_id, orientation).is_none(), "Found predecessors for non-existent node {} {}", node_id, name);
        assert!(gbz.successors(node_id, orientation).is_none(), "Found successors for non-existent node {} {}", node_id, name);
    }
}

//-----------------------------------------------------------------------------

#[test]
fn statistics() {
    let filename = support::get_test_data("example.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    assert_eq!(gbz.nodes(), 12, "Invalid number of nodes");
    assert_eq!(gbz.min_node(), 11, "Invalid minimum node identifier");
    assert_eq!(gbz.max_node(), 25, "Invalid maximum node identifier");
    assert!(!gbz.has_translation(), "The graph should not contain a translation");
}

#[test]
fn serialize() {
    let filename = support::get_test_data("example.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();
    serialize::test(&gbz, "gbz", None, true);
}

#[test]
fn nodes() {
    let filename = support::get_test_data("example.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    let mut truth: BTreeMap<usize, String> = BTreeMap::new();
    let true_nodes = [
        (11, "G"), (12, "A"), (13, "T"), (14, "T"), (15, "A"), (16, "C"), (17, "A"),
        (21, "G"), (22, "A"), (23, "T"), (24, "T"), (25, "A"),
    ];
    for (key, value) in true_nodes.iter() {
        truth.insert(*key, value.to_string());
    }

    // Random access.
    for node_id in 0..gbz.max_node() + 2 {
        let should_exist = truth.contains_key(&node_id);
        assert_eq!(gbz.has_node(node_id), should_exist, "Invalid has_node({}) result", node_id);
        if should_exist {
            assert_eq!(gbz.sequence(node_id), truth.get(&node_id).map(|s| s.as_bytes()), "Invalid sequence for node {}", node_id);
        } else {
            assert!(gbz.sequence(node_id).is_none(), "Got a sequence for non-existent node {}", node_id);
        }
    }

    // Iterate forward.
    assert_eq!(gbz.node_iter().len(), truth.len(), "Invalid iterator length");
    assert!(gbz.node_iter().eq(truth.iter().map(|(k, _)| *k)), "Invalid iterator (forward)");

    // Iterate backward.
    assert!(gbz.node_iter().rev().eq(truth.iter().rev().map(|(k, _)| *k)), "Invalid iterator (backward)");
}

#[test]
fn edges() {
    let filename = support::get_test_data("example.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    let edges: Vec<(usize, Orientation, usize, Orientation)> = vec![
        (11, Orientation::Forward, 12, Orientation::Forward),
        (11, Orientation::Forward, 13, Orientation::Forward),
        (12, Orientation::Forward, 14, Orientation::Forward),
        (13, Orientation::Forward, 14, Orientation::Forward),
        (14, Orientation::Forward, 15, Orientation::Forward),
        (14, Orientation::Forward, 16, Orientation::Forward),
        (15, Orientation::Forward, 17, Orientation::Forward),
        (16, Orientation::Forward, 17, Orientation::Forward),
        (21, Orientation::Forward, 22, Orientation::Forward),
        (21, Orientation::Forward, 23, Orientation::Forward),
        (22, Orientation::Forward, 24, Orientation::Forward),
        (23, Orientation::Forward, 24, Orientation::Reverse),
        (24, Orientation::Forward, 25, Orientation::Forward),
    ];
    let (predecessors, successors) = get_pred_succ(&gbz, &edges);

    for node_id in 0..gbz.max_node() + 2 {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            check_pred_succ(&gbz, &predecessors, &successors, node_id, orientation);
        }
    }
}

//-----------------------------------------------------------------------------

// FIXME translation.gbz: statistics, serialize, nodes, edges, segments, links (include reverse iterators)

//-----------------------------------------------------------------------------
