use super::*;

use crate::support;

use simple_sds::serialize;

use std::collections::{BTreeSet, BTreeMap};

//-----------------------------------------------------------------------------

fn check_nodes(gbz: &GBZ, true_nodes: &[(usize, &str)]) {
    let mut truth: BTreeMap<usize, String> = BTreeMap::new();
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

type Neighbors = BTreeMap<(usize, Orientation), BTreeSet<(usize, Orientation)>>;

fn get_pred_succ<T: Iterator<Item = usize>>(edges: &Vec<(usize, Orientation, usize, Orientation)>, iter: T) -> (Neighbors, Neighbors) {
    let mut predecessors: BTreeMap<(usize, Orientation), BTreeSet<(usize, Orientation)>> = BTreeMap::new();
    let mut successors: BTreeMap<(usize, Orientation), BTreeSet<(usize, Orientation)>> = BTreeMap::new();

    // Initialize the predecessor/successor lists for all valid ids.
    for id in iter {
        predecessors.insert((id, Orientation::Forward), BTreeSet::new());
        predecessors.insert((id, Orientation::Reverse), BTreeSet::new());
        successors.insert((id, Orientation::Forward), BTreeSet::new());
        successors.insert((id, Orientation::Reverse), BTreeSet::new());
    }

    // Add edges in both orientations.
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

    let true_nodes = [
        (11, "G"), (12, "A"), (13, "T"), (14, "T"), (15, "A"), (16, "C"), (17, "A"),
        (21, "G"), (22, "A"), (23, "T"), (24, "T"), (25, "A"),
    ];
    check_nodes(&gbz, &true_nodes);
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
    let (predecessors, successors) = get_pred_succ(&edges, gbz.node_iter());

    for node_id in 0..gbz.max_node() + 2 {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            check_pred_succ(&gbz, &predecessors, &successors, node_id, orientation);
        }
    }
}

#[test]
fn no_translation() {
    let filename = support::get_test_data("example.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    assert!(!gbz.has_translation(), "The graph should not contain a translation");
    assert!(gbz.segment_iter().is_none(), "The graph should not have a segment iterator");
    for node_id in gbz.min_node()..=gbz.max_node() {
        assert!(gbz.node_to_segment(node_id).is_none(), "The graph should not have a segment for node {}", node_id);
    }
}

//-----------------------------------------------------------------------------

#[test]
fn statistics_trans() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    assert_eq!(gbz.nodes(), 9, "Invalid number of nodes");
    assert_eq!(gbz.min_node(), 1, "Invalid minimum node identifier");
    assert_eq!(gbz.max_node(), 11, "Invalid maximum node identifier");
}

#[test]
fn serialize_trans() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();
    serialize::test(&gbz, "gbz-trans", None, true);
}

#[test]
fn nodes_trans() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    let true_nodes = [
        (1, "GA"), (2, "T"), (3, "T"), (4, "A"), (5, "CA"), (6, "G"),
        (9, "A"), (10, "T"), (11, "TA"),
    ];
    check_nodes(&gbz, &true_nodes);
}

#[test]
fn edges_trans() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    let edges: Vec<(usize, Orientation, usize, Orientation)> = vec![
        (1, Orientation::Forward, 2, Orientation::Forward),
        (2, Orientation::Forward, 3, Orientation::Forward),
        (2, Orientation::Forward, 4, Orientation::Forward),
        (3, Orientation::Forward, 5, Orientation::Forward),
        (4, Orientation::Forward, 5, Orientation::Forward),
        (5, Orientation::Forward, 6, Orientation::Forward),
        (6, Orientation::Forward, 9, Orientation::Forward),
        (6, Orientation::Forward, 10, Orientation::Forward),
        (9, Orientation::Forward, 11, Orientation::Forward),
        (10, Orientation::Forward, 11, Orientation::Forward),
    ];
    let (predecessors, successors) = get_pred_succ(&edges, gbz.node_iter());

    for node_id in 0..gbz.max_node() + 2 {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            check_pred_succ(&gbz, &predecessors, &successors, node_id, orientation);
        }
    }
}

#[test]
fn segments() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    assert!(gbz.has_translation(), "The graph does not contain node-to-segment translation");

    let truth = vec![
        Segment::from_fields(0, "s11".as_bytes(), 1..3, "GAT".as_bytes()),
        Segment::from_fields(1, "s12".as_bytes(), 3..4, "T".as_bytes()),
        Segment::from_fields(2, "s13".as_bytes(), 4..5, "A".as_bytes()),
        Segment::from_fields(3, "s14".as_bytes(), 5..7, "CAG".as_bytes()),
        Segment::from_fields(5, "s15".as_bytes(), 9..10, "A".as_bytes()),
        Segment::from_fields(6, "s16".as_bytes(), 10..11, "T".as_bytes()),
        Segment::from_fields(7, "s17".as_bytes(), 11..12, "TA".as_bytes()),
    ];

    // Iterate forward.
    if let Some(iter) = gbz.segment_iter() {
        assert!(iter.eq(truth.iter().cloned()), "Invalid iterator (forward)");
    } else {
        panic!("Could not create a segment iterator");
    }

    // Iterate backward.
    assert!(gbz.segment_iter().unwrap().eq(truth.iter().cloned()), "Invalid iterator (backward)");
}

#[test]
fn links() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    assert!(gbz.has_translation(), "The graph does not contain node-to-segment translation");

    // Use segment identifiers as a shorthand for segments.
    let links: Vec<(usize, Orientation, usize, Orientation)> = vec![
        (0, Orientation::Forward, 1, Orientation::Forward),
        (0, Orientation::Forward, 2, Orientation::Forward),
        (1, Orientation::Forward, 3, Orientation::Forward),
        (2, Orientation::Forward, 3, Orientation::Forward),
        (3, Orientation::Forward, 5, Orientation::Forward),
        (3, Orientation::Forward, 6, Orientation::Forward),
        (5, Orientation::Forward, 7, Orientation::Forward),
        (6, Orientation::Forward, 7, Orientation::Forward),
    ];
    let (predecessors, successors) = get_pred_succ(&links, gbz.segment_iter().unwrap().map(|s| s.id));

    // Validate the links.
    for segment in gbz.segment_iter().unwrap() {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let name = if orientation == Orientation::Forward { "(forward)" } else { "(reverse)" };
            let truth = predecessors.get(&(segment.id, orientation)).unwrap();

            if let Some(iter) = gbz.segment_predecessors(&segment, orientation) {
                assert!(iter.map(|(s, o)| (s.id, o)).eq(truth.iter().cloned()), "Invalid predecessors for segment {} {}", segment.id, name);
            } else {
                panic!("Could not get predecessors for segment {} {}", segment.id, name);
            }
            for (s, _) in gbz.segment_predecessors(&segment, orientation).unwrap() {
                assert_eq!(s, gbz.graph.segment(s.id), "Invalid predecessor segment {} for segment {} {}", s.id, segment.id, name);
            }

            let truth = successors.get(&(segment.id, orientation)).unwrap();
            if let Some(iter) = gbz.segment_successors(&segment, orientation) {
                assert!(iter.map(|(s, o)| (s.id, o)).eq(truth.iter().cloned()), "Invalid successors for segment {} {}", segment.id, name);
            } else {
                panic!("Could not get successors for segment {} {}", segment.id, name);
            }
            for (s, _) in gbz.segment_successors(&segment, orientation).unwrap() {
                assert_eq!(s, gbz.graph.segment(s.id), "Invalid successor segment {} for segment {} {}", s.id, segment.id, name);
            }
        }
    }

}

//-----------------------------------------------------------------------------
