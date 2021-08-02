use super::*;

use simple_sds::serialize;

//-----------------------------------------------------------------------------

// GBWT example from the paper: (edges, runs, invalid_node)
fn get_edges_runs() -> (Vec<Vec<(usize, usize)>>, Vec<Vec<(usize, usize)>>, usize) {
    let edges = vec![
        vec![(1, 0)],
        vec![(2, 0), (3, 0)],
        vec![(4, 0), (5, 0)],
        vec![(4, 1)],
        vec![(5, 1), (6, 0)],
        vec![(7, 0)],
        vec![(7, 2)],
        vec![(0, 0)],
    ];
    let runs = vec![
        vec![(0, 3)],
        vec![(0, 2), (1, 1)],
        vec![(0, 1), (1, 1)],
        vec![(0, 1)],
        vec![(1, 1), (0, 1)],
        vec![(0, 2)],
        vec![(0, 1)],
        vec![(0, 3)],
    ];
    (edges, runs, 8)
}

// Bidirectional version of the example: (edges, runs, invalid_node)
fn get_bidirectional() -> (Vec<Vec<(usize, usize)>>, Vec<Vec<(usize, usize)>>, usize) {
    let edges = vec![
        // ENDMARKER
        vec![(2, 0), (15, 0)],
        // 1
        vec![(4, 0), (6, 0)],
        vec![(0, 0)],
        // 2
        vec![(8, 0), (10, 0)],
        vec![(3, 0)],
        // 3
        vec![(8, 1)],
        vec![(3, 2)],
        // 4
        vec![(10, 1), (12, 0)],
        vec![(5, 0), (7, 0)],
        // 5
        vec![(14, 0)],
        vec![(5, 1), (9, 0)],
        // 6
        vec![(14, 2)],
        vec![(9, 1)],
        // 7
        vec![(0, 0)],
        vec![(11, 0), (13, 0)],
    ];
    let runs = vec![
        // ENDMARKER
        vec![(0, 3), (1, 3)],
        // 1
        vec![(0, 2), (1, 1)],
        vec![(0, 3)],
        // 2
        vec![(0, 1), (1, 1)],
        vec![(0, 2)],
        // 3
        vec![(0, 1)],
        vec![(0, 1)],
        // 4
        vec![(1, 1), (0, 1)],
        vec![(1, 1), (0, 1)],
        // 5
        vec![(0, 2)],
        vec![(0, 1), (1, 1)],
        // 6
        vec![(0, 1)],
        vec![(0, 1)],
        // 7
        vec![(0, 3)],
        vec![(1, 1), (0, 2)],
    ];
    (edges, runs, 16)
}

fn create_bwt(edges: &[Vec<(usize, usize)>], runs: &[Vec<(usize, usize)>]) -> BWT {
    let mut builder = BWTBuilder::new();
    assert_eq!(builder.len(), 0, "Newly created builder has non-zero length");
    assert!(builder.is_empty(), "Newly created builder is not empty");

    for i in 0..edges.len() {
        builder.append(&edges[i], &runs[i]);
    }
    assert_eq!(builder.len(), edges.len(), "Invalid number of records in the builder");
    assert_eq!(builder.is_empty(), edges.is_empty(), "Invalid builder emptiness");

    BWT::from(builder)
}

// Check records in the BWT, using the provided edges as the source of truth.
// Also checks that `id()` works correctly.
fn check_records(bwt: &BWT, edges: &[Vec<(usize, usize)>]) {
    assert_eq!(bwt.len(), edges.len(), "Invalid number of records in the BWT");
    assert_eq!(bwt.is_empty(), edges.is_empty(), "Invalid BWT emptiness");

    // Edges.
    for i in 0..bwt.len() {
        let record = bwt.record(i);
        let curr_edges = &edges[i];
        assert_eq!(record.is_none(), curr_edges.is_empty(), "Invalid record {} existence", i);
        if let Some(record) = record {
            assert_eq!(record.id(), i, "Invalid id for record {}", i);
            assert_eq!(record.outdegree(), curr_edges.len(), "Invalid outdegree in record {}", i);
            for j in 0..record.outdegree() {
                assert_eq!(record.successor(j), curr_edges[j].0, "Invalid successor {} in record {}", j, i);
                assert_eq!(record.offset(j), curr_edges[j].1, "Invalid offset {} in record {}", j, i);
            }
        }
    }
}

// Check that the iterator finds the correct records.
fn check_iter(bwt: &BWT) {
    let mut iter = bwt.iter();
    for i in 0..bwt.len() {
        if let Some(truth) = bwt.record(i) {
            if let Some(record) = iter.next() {
                assert_eq!(record.id(), truth.id(), "Invalid record id from the iterator");
            } else {
                panic!("Iterator did not find record {}", i);
            }
        }
    }
    assert!(iter.next().is_none(), "Iterator found a record past the end");
}

// Check all `lf()` results in the BWT, using the provided edges and runs as the source of truth.
// Also checks that `offset_to()` works in positive cases and that `len()` is correct.
fn check_lf(bwt: &BWT, edges: &[Vec<(usize, usize)>], runs: &[Vec<(usize, usize)>]) {
    // `lf()` at each offset of each record.
    for i in 0..bwt.len() {
        if let Some(record) = bwt.record(i) {
            let mut offset = 0;
            let mut curr_edges = edges[i].clone();
            let curr_runs = &runs[i];
            for (rank, len) in curr_runs {
                for _ in 0..*len {
                    let edge = curr_edges[*rank];
                    let expected = if edge.0 == ENDMARKER { None } else { Some(edge) };
                    assert_eq!(record.lf(offset), expected, "Invalid lf({}) in record {}", offset, i);
                    let expected = if edge.0 == ENDMARKER { None } else { Some(offset) };
                    assert_eq!(record.offset_to(edge), expected, "Invalid offset_to(({}, {})) in record {}", edge.0, edge.1, i);
                    offset += 1;
                    curr_edges[*rank].1 += 1;
                }
            }
            assert_eq!(record.len(), offset, "Invalid record {} length", i);
            assert_eq!(record.lf(offset), None, "Got an lf() result past the end in record {}", i);
        }
    }
}

// Check all `follow()` results in the BWT, using `lf()` as the source of truth.
fn check_follow(bwt: &BWT, invalid_node: usize) {
    for record in bwt.iter() {
        let i = record.id();
        // Check all ranges, including empty and past-the-end ones.
        let len = record.len();
        for start in 0..len + 1 {
            for limit in start..len + 1 {
                // With an endmarker.
                assert_eq!(record.follow(start..limit, ENDMARKER), None, "Got a follow({}..{}, endmarker) result in record {}", start, limit, i);

                // With each successor node.
                for rank in 0..record.outdegree() {
                    let successor = record.successor(rank);
                    if successor == ENDMARKER {
                        continue;
                    }
                    if let Some(result) = record.follow(start..limit, successor) {
                        let mut found = result.start..result.start;
                        for j in start..limit {
                            if let Some((node, offset)) = record.lf(j) {
                                if node == successor && offset == found.end {
                                    found.end += 1;
                                }
                            }
                        }
                        assert_eq!(result, found, "follow({}..{}, {}) did not find the correct range in record {}", start, limit, successor, i);
                    } else {
                        for j in start..limit {
                            if let Some((node, _)) = record.lf(j) {
                                assert_ne!(node, successor, "follow({}..{}, {}) did not follow offset {} in record {}", start, limit, successor, j, i);
                            }
                        }
                    }
                }

                // With an invalid node.
                assert_eq!(record.follow(start..limit, invalid_node), None, "Got a follow({}..{}, invalid) result in record {}", start, limit, i);
            }
        }
    }
}

// Check negative cases for `offset_to()`.
fn negative_offset_to(bwt: &BWT, invalid_node: usize) {
    for record in bwt.iter() {
        assert_eq!(record.offset_to((ENDMARKER, 0)), None, "Got an offset to the endmarker from record {}", record.id());
        assert_eq!(record.offset_to((invalid_node, 0)), None, "Got an offset to an invalid node from record {}", record.id());
        for rank in 0..record.outdegree() {
            let successor = record.successor(rank);
            if successor == ENDMARKER {
                continue;
            }
            let offset = record.offset(rank);
            if offset > 0 {
                assert_eq!(record.offset_to((successor, offset - 1)), None, "Got an offset from record {} to a too small position in {}", record.id(), successor);
            }
            let count = record.follow(0..record.len(), successor).unwrap().len();
            assert_eq!(record.offset_to((successor, offset + count)), None, "Got an offset from record {} to a too large position in {}", record.id(), successor);
        }
    }
}

//-----------------------------------------------------------------------------

#[test]
fn empty_bwt() {
    let edges = Vec::new();
    let runs = Vec::new();
    let invalid_node = 0;
    let bwt = create_bwt(&edges, &runs);
    check_records(&bwt, &edges);
    check_iter(&bwt);
    check_lf(&bwt, &edges, &runs);
    check_follow(&bwt, invalid_node);
    negative_offset_to(&bwt, invalid_node);
    serialize::test(&bwt, "empty-bwt", None, true);
}

#[test]
fn non_empty_bwt() {
    let (edges, runs, invalid_node) = get_edges_runs();
    let bwt = create_bwt(&edges, &runs);
    check_records(&bwt, &edges);
    check_iter(&bwt);
    check_lf(&bwt, &edges, &runs);
    check_follow(&bwt, invalid_node);
    negative_offset_to(&bwt, invalid_node);
    serialize::test(&bwt, "non-empty-bwt", None, true);
}

#[test]
fn empty_records() {
    let (mut edges, mut runs, invalid_node) = get_edges_runs();
    edges[2] = Vec::new();
    edges[6] = Vec::new();
    runs[2] = Vec::new();
    runs[6] = Vec::new();
 
    let bwt = create_bwt(&edges, &runs);
    check_records(&bwt, &edges);
    check_iter(&bwt);
    check_lf(&bwt, &edges, &runs);
    check_follow(&bwt, invalid_node);
    negative_offset_to(&bwt, invalid_node);
    serialize::test(&bwt, "bwt-with-empty", None, true);
}

#[test]
fn bidirectional_bwt() {
    let (edges, runs, invalid_node) = get_bidirectional();
    let bwt = create_bwt(&edges, &runs);
    check_records(&bwt, &edges);
    check_iter(&bwt);
    check_lf(&bwt, &edges, &runs);
    check_follow(&bwt, invalid_node);
    negative_offset_to(&bwt, invalid_node);
    serialize::test(&bwt, "bidirectional-bwt", None, true);
}

//-----------------------------------------------------------------------------
