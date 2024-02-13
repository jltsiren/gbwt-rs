use super::*;

use crate::support;

use simple_sds::serialize;

use std::collections::{BTreeSet, BTreeMap, HashSet};

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

fn get_pred_succ<T: Iterator<Item = usize>>(edges: &[(usize, Orientation, usize, Orientation)], iter: T) -> (Neighbors, Neighbors) {
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
    if gbz.has_node(node_id) {
        let truth = predecessors.get(&(node_id, orientation)).unwrap();
        let iter = gbz.predecessors(node_id, orientation);
        assert!(iter.is_some(), "Could not find predecessors for node {} ({})", node_id, orientation);
        assert!(iter.unwrap().eq(truth.iter().cloned()), "Invalid predecessors for node {} ({})", node_id, orientation);
        let truth = successors.get(&(node_id, orientation)).unwrap();
        let iter = gbz.successors(node_id, orientation);
        assert!(iter.is_some(), "Could not find successors for node {} ({})", node_id, orientation);
        assert!(iter.unwrap().eq(truth.iter().cloned()), "Invalid successors for node {} ({})", node_id, orientation);
    } else {
        assert!(gbz.predecessors(node_id, orientation).is_none(), "Found predecessors for non-existent node {} ({})", node_id, orientation);
        assert!(gbz.successors(node_id, orientation).is_none(), "Found successors for non-existent node {} ({})", node_id, orientation);
    }
}

fn check_paths(gbz: &GBZ, truth: &[Vec<(usize, Orientation)>]) {
    assert_eq!(gbz.paths(), truth.len(), "Invalid number of paths");
    for id in 0..gbz.paths() {
        let iter = gbz.path(id, Orientation::Forward);
        assert!(iter.is_some(), "Could not get path {} (forward)", id);
        assert!(iter.unwrap().eq(truth[id].iter().cloned()), "Invalid path {} (forward)", id);
        let iter = gbz.path(id, Orientation::Reverse);
        assert!(iter.is_some(), "Could not get path {} (reverse)", id);
        assert!(iter.unwrap().eq(truth[id].iter().rev().map(|(v, o)| (*v, o.flip()))), "Invalid path {} (reverse)", id);
    }

    assert!(gbz.path(truth.len(), Orientation::Forward).is_none(), "Got a past-the-end path {} (forward)", truth.len());
    assert!(gbz.path(truth.len(), Orientation::Reverse).is_none(), "Got a past-the-end path {} (reverse)", truth.len());
}

fn check_states(gbz: &GBZ) {
    let mut stack: Vec<BidirectionalState> = Vec::new();
    let mut visited: HashSet<BidirectionalState> = HashSet::new();

    // Start from every node in every orientation.
    for node_id in 0..gbz.max_node() + 2 {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let state = gbz.search_state(node_id, orientation);
            if gbz.has_node(node_id) {
                let truth = gbz.index.bd_find(support::encode_node(node_id, orientation));
                assert_eq!(state, truth, "Invalid search state for node {} ({})", node_id, orientation);
                if let Some(state) = state {
                    stack.push(state.clone());
                    visited.insert(state);
                }
            } else {
                assert!(state.is_none(), "Got a search state for nonexistent node {} ({})", node_id, orientation);
            }
        }
    }

    // Find all extensions of each state and compare them to those obtained from the index.
    // Put unvisited states into stack.
    while !stack.is_empty() {
        let state = stack.pop().unwrap();

        // Forward extensions.
        let mut truth: HashSet<BidirectionalState> = HashSet::new();
        let (last_id, last_o) = state.to();
        for (node_id, orientation) in gbz.successors(last_id, last_o).unwrap() {
            if let Some(successor) = gbz.index.extend_forward(&state, support::encode_node(node_id, orientation)) {
                truth.insert(successor);
            }
        }
        if let Some(iter) = gbz.follow_forward(&state) {
            let found: HashSet<BidirectionalState> = iter.collect();
            assert_eq!(found, truth, "Invalid successors from {:?}", state);
            for succ in found.iter() {
                if !visited.contains(succ) {
                    stack.push(succ.clone());
                    visited.insert(succ.clone());
                }
            }
        } else {
            panic!("Could not follow paths forward from {:?}", state);
        }

        // Backward extensions.
        let mut truth: HashSet<BidirectionalState> = HashSet::new();
        let (first_id, first_o) = state.from();
        for (node_id, orientation) in gbz.predecessors(first_id, first_o).unwrap() {
            if let Some(predecessor) = gbz.index.extend_backward(&state, support::encode_node(node_id, orientation)) {
                truth.insert(predecessor);
            }
        }
        if let Some(iter) = gbz.follow_backward(&state) {
            let found: HashSet<BidirectionalState> = iter.collect();
            assert_eq!(found, truth, "Invalid predecessors from {:?}", state);
            for succ in found.iter() {
                if !visited.contains(succ) {
                    stack.push(succ.clone());
                    visited.insert(succ.clone());
                }
            }
        } else {
            panic!("Could not follow paths backward from {:?}", state);
        }
    }
}

//-----------------------------------------------------------------------------

#[test]
fn reference_samples() {
    let filename = support::get_test_data("example.gbz");
    let mut gbz: GBZ = serialize::load_from(&filename).unwrap();

    // No reference samples in the original graph.
    assert!(gbz.reference_sample_ids(false).is_empty(), "The graph should not have reference sample ids");
    assert!(gbz.reference_sample_names(false).is_empty(), "The graph should not have reference sample names");
    assert_eq!(gbz.reference_sample_ids(true), vec![0], "Invalid reference + generic sample ids (default)");
    assert_eq!(gbz.reference_sample_names(true), vec![String::from(REF_SAMPLE)], "Invalid reference + generic sample names (default)");

    // Add a reference sample.
    let ref_samples = vec![String::from("sample")];
    let ref_and_generic = vec![String::from("sample"), String::from(REF_SAMPLE)];
    assert_eq!(gbz.set_reference_samples(&ref_samples), 1, "Could not set reference samples (correct)");
    assert_eq!(gbz.reference_sample_ids(false), vec![1], "Invalid reference sample ids (correct)");
    assert_eq!(gbz.reference_sample_names(false), ref_samples, "Invalid reference sample names (correct)");
    assert_eq!(gbz.reference_sample_ids(true), vec![1, 0], "Invalid reference + generic sample ids (correct)");
    assert_eq!(gbz.reference_sample_names(true), ref_and_generic, "Invalid reference + generic sample names (correct)");

    // Try setting duplicate names.
    let duplicates = vec![String::from("sample"), String::from("sample")];
    assert_eq!(gbz.set_reference_samples(&duplicates), 1, "Could not set reference samples (duplicates)");
    assert_eq!(gbz.reference_sample_ids(false), vec![1], "Invalid reference sample ids (duplicates)");
    assert_eq!(gbz.reference_sample_names(false), ref_samples, "Invalid reference sample names (duplicates)");
    assert_eq!(gbz.reference_sample_ids(true), vec![1, 0], "Invalid reference + generic sample ids (duplicates)");
    assert_eq!(gbz.reference_sample_names(true), ref_and_generic, "Invalid reference + generic sample names (duplicates)");

    // Try nonexistent names.
    let nonexistent = vec![String::from("nonexistent")];
    let generic_only = vec![String::from(REF_SAMPLE)];
    assert_eq!(gbz.set_reference_samples(&nonexistent), 0, "Could not set reference samples (nonexistent)");
    assert_eq!(gbz.reference_sample_ids(false), Vec::new(), "Invalid reference sample ids (nonexistent)");
    assert_eq!(gbz.reference_sample_names(false), Vec::<String>::new(), "Invalid reference sample names (nonexistent)");
    assert_eq!(gbz.reference_sample_ids(true), vec![0], "Invalid reference + generic sample ids (nonexistent)");
    assert_eq!(gbz.reference_sample_names(true), generic_only, "Invalid reference + generic sample names (nonexistent)");

    // Try setting the generic sample as a reference sample.
    assert_eq!(gbz.set_reference_samples(&generic_only), 0, "Could not set reference samples (generic)");
    assert_eq!(gbz.reference_sample_ids(false), Vec::new(), "Invalid reference sample ids (generic)");
    assert_eq!(gbz.reference_sample_names(false), Vec::<String>::new(), "Invalid reference sample names (generic)");
    assert_eq!(gbz.reference_sample_ids(true), vec![0], "Invalid reference + generic sample ids (generic)");
    assert_eq!(gbz.reference_sample_names(true), generic_only, "Invalid reference + generic sample names (generic)");
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
    assert!(GBZ::is_gbz(&filename), "File {} is not a GBZ file", filename.to_string_lossy());
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
fn paths() {
    let filename = support::get_test_data("example.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    let paths: Vec<Vec<(usize, Orientation)>> = vec![
        vec![(11, Orientation::Forward), (12, Orientation::Forward), (14, Orientation::Forward), (15, Orientation::Forward), (17, Orientation::Forward)],
        vec![(21, Orientation::Forward), (22, Orientation::Forward), (24, Orientation::Forward), (25, Orientation::Forward)],
        vec![(11, Orientation::Forward), (12, Orientation::Forward), (14, Orientation::Forward), (15, Orientation::Forward), (17, Orientation::Forward)],
        vec![(11, Orientation::Forward), (13, Orientation::Forward), (14, Orientation::Forward), (16, Orientation::Forward), (17, Orientation::Forward)],
        vec![(21, Orientation::Forward), (22, Orientation::Forward), (24, Orientation::Forward), (23, Orientation::Reverse), (21, Orientation::Reverse)],
        vec![(21, Orientation::Forward), (22, Orientation::Forward), (24, Orientation::Forward), (25, Orientation::Forward)],
    ];
    check_paths(&gbz, &paths);
}

#[test]
fn states() {
    let filename = support::get_test_data("example.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();
    check_states(&gbz);
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
fn paths_trans() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    let paths: Vec<Vec<(usize, Orientation)>> = vec![
        vec![(1, Orientation::Forward), (2, Orientation::Forward), (3, Orientation::Forward), (5, Orientation::Forward), (6, Orientation::Forward), (9, Orientation::Forward), (11, Orientation::Forward)],
        vec![(1, Orientation::Forward), (2, Orientation::Forward), (3, Orientation::Forward), (5, Orientation::Forward), (6, Orientation::Forward), (9, Orientation::Forward), (11, Orientation::Forward)],
        vec![(1, Orientation::Forward), (2, Orientation::Forward), (4, Orientation::Forward), (5, Orientation::Forward), (6, Orientation::Forward), (10, Orientation::Forward), (11, Orientation::Forward)],
    ];
    check_paths(&gbz, &paths);
}

#[test]
fn states_trans() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();
    check_states(&gbz);
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
            let truth = predecessors.get(&(segment.id, orientation)).unwrap();

            if let Some(iter) = gbz.segment_predecessors(&segment, orientation) {
                assert!(iter.map(|(s, o)| (s.id, o)).eq(truth.iter().cloned()), "Invalid predecessors for segment {} ({})", segment.id, orientation);
            } else {
                panic!("Could not get predecessors for segment {} {}", segment.id, orientation);
            }
            for (s, _) in gbz.segment_predecessors(&segment, orientation).unwrap() {
                assert_eq!(s, gbz.graph.segment(s.id), "Invalid predecessor segment {} for segment {} ({})", s.id, segment.id, orientation);
            }

            let truth = successors.get(&(segment.id, orientation)).unwrap();
            if let Some(iter) = gbz.segment_successors(&segment, orientation) {
                assert!(iter.map(|(s, o)| (s.id, o)).eq(truth.iter().cloned()), "Invalid successors for segment {} ({})", segment.id, orientation);
            } else {
                panic!("Could not get successors for segment {} {}", segment.id, orientation);
            }
            for (s, _) in gbz.segment_successors(&segment, orientation).unwrap() {
                assert_eq!(s, gbz.graph.segment(s.id), "Invalid successor segment {} for segment {} ({})", s.id, segment.id, orientation);
            }
        }
    }

}

#[test]
fn segment_paths() {
    let filename = support::get_test_data("translation.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    // Paths as sequences of segment identifiers.
    let paths: Vec<Vec<(usize, Orientation)>> = vec![
        vec![(0, Orientation::Forward), (1, Orientation::Forward), (3, Orientation::Forward), (5, Orientation::Forward), (7, Orientation::Forward)],
        vec![(0, Orientation::Forward), (1, Orientation::Forward), (3, Orientation::Forward), (5, Orientation::Forward), (7, Orientation::Forward)],
        vec![(0, Orientation::Forward), (2, Orientation::Forward), (3, Orientation::Forward), (6, Orientation::Forward), (7, Orientation::Forward)],
    ];

    // Check the paths.
    assert_eq!(gbz.paths(), paths.len(), "Invalid number of paths");
    for id in 0..gbz.paths() {
        let iter = gbz.segment_path(id, Orientation::Forward);
        assert!(iter.is_some(), "Could not get path {} (forward)", id);
        assert!(iter.unwrap().map(|(s, o)| (s.id, o)).eq(paths[id].iter().cloned()), "Invalid path {} (forward)", id);
        let iter = gbz.segment_path(id, Orientation::Reverse);
        assert!(iter.is_some(), "Could not get path {} (reverse)", id);
        assert!(iter.unwrap().map(|(s, o)| (s.id, o)).eq(paths[id].iter().rev().map(|(v, o)| (*v, o.flip()))), "Invalid path {} (reverse)", id);
    }
    assert!(gbz.segment_path(paths.len(), Orientation::Forward).is_none(), "Got a past-the-end path {} (forward)", paths.len());
    assert!(gbz.segment_path(paths.len(), Orientation::Reverse).is_none(), "Got a past-the-end path {} (reverse)", paths.len());

    // Check the segments on the paths.
    for id in 0..gbz.paths() {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            for (s, _) in gbz.segment_path(id, orientation).unwrap() {
                assert_eq!(s, gbz.graph.segment(s.id), "Invalid segment {} on path {} ({})", s.id, id, orientation);
            }
        }
    }
}

//-----------------------------------------------------------------------------

#[test]
fn weakly_connected_components() {
    let filename = support::get_test_data("example.gbz");
    let gbz: GBZ = serialize::load_from(&filename).unwrap();

    let truth: Vec<Vec<usize>> = vec![
        vec![11, 12, 13, 14, 15, 16, 17],
        vec![21, 22, 23, 24, 25],
    ];

    let components = gbz.weakly_connected_components();
    assert_eq!(components.len(), truth.len(), "Invalid number of components");
    for i in 0..components.len() {
        assert_eq!(components[i], truth[i], "Invalid component {}", i);
    }
}

#[test]
fn reference_positions() {
    let filename = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&filename).unwrap();
    let index: &GBWT = graph.as_ref();
    let metadata = graph.metadata().unwrap();
    let ref_samples: BTreeSet<usize> = graph.reference_sample_ids(true).into_iter().collect();

    // For each reference path, collect the starting positions of all nodes.
    let mut node_starts: Vec<(usize, Vec<(usize, Pos)>)> = Vec::new();
    for (path_id, path_name) in metadata.path_iter().enumerate() {
        if !ref_samples.contains(&path_name.sample()) {
            continue;
        }
        let mut path_offset = 0;
        let mut curr = index.start(support::encode_path(path_id, Orientation::Forward));
        let mut positions: Vec<(usize, Pos)> = Vec::new();
        while let Some(p) = curr {
            positions.push((path_offset, p));
            path_offset += graph.sequence_len(support::node_id(p.node)).unwrap();
            curr = index.forward(p);
        }
        node_starts.push((path_id, positions));
    }

    // Check the indexed positions with various intervals.
    for interval in 0..10 {
        let paths = graph.reference_positions(interval, false);
        assert_eq!(paths.len(), node_starts.len(), "Wrong number of reference positions for interval {}", interval);
        for i in 0..paths.len() {
            assert_eq!(paths[i].0, node_starts[i].0, "Wrong path id at offset {} with interval {}", i, interval);
            let mut next = 0;
            let mut iter = paths[i].1.iter();
            for (offset, pos) in node_starts[i].1.iter() {
                if *offset >= next {
                    let indexed = iter.next();
                    assert!(indexed.is_some(), "Missing indexed position for path {} with interval {}", paths[i].0, interval);
                    let indexed = indexed.unwrap();
                    assert_eq!(indexed.0, *offset, "Wrong indexed offset for path {} with interval {}", paths[i].0, interval);
                    assert_eq!(indexed.1, *pos, "Wrong indexed GBWT position for path {} with interval {}", paths[i].0, interval);
                    next = offset + interval;
                }
            }
            assert!(iter.next().is_none(), "Too many indexed positions for path {} with interval {}", paths[i].0, interval);
        }
    }
}

//-----------------------------------------------------------------------------
