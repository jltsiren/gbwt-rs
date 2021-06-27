use super::*;

use simple_sds::serialize;

//-----------------------------------------------------------------------------

fn get_edges() -> Vec<Vec<(usize, usize)>> {
    vec![
        vec![(1, 0)],
        vec![(2, 0), (3, 0)],
        vec![(4, 0), (5, 0)],
        vec![(4, 1)],
        vec![(5, 1), (6, 0)],
        vec![(7, 0)],
        vec![(7, 2)],
        vec![(0, 0)],
    ]
}

fn get_runs() -> Vec<Vec<(usize, usize)>> {
    vec![
        vec![(0, 3)],
        vec![(0, 2), (1, 1)],
        vec![(0, 1), (1, 1)],
        vec![(0, 1)],
        vec![(1, 1), (0, 1)],
        vec![(0, 2)],
        vec![(0, 1)],
        vec![(0, 3)],
    ]
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

// Check all edges in the BWT, using the provided edges as the source of truth.
fn check_edges(bwt: &BWT, edges: &[Vec<(usize, usize)>]) {
    assert_eq!(bwt.len(), edges.len(), "Invalid number of records in the BWT");
    assert_eq!(bwt.is_empty(), edges.is_empty(), "Invalid BWT emptiness");

    // Edges.
    for i in 0..bwt.len() {
        let record = bwt.record(i);
        let curr_edges = &edges[i];
        assert_eq!(record.is_none(), curr_edges.is_empty(), "Invalid record {} existence", i);
        if let Some(record) = record {
            assert_eq!(record.outdegree(), curr_edges.len(), "Invalid outdegree in record {}", i);
            for j in 0..record.outdegree() {
                assert_eq!(record.successor(j), curr_edges[j].0, "Invalid successor {} in record {}", j, i);
                assert_eq!(record.offset(j), curr_edges[j].1, "Invalid offset {} in record {}", j, i);
            }
        }
    }
}

// Check all `lf()` results in the BWT, using the provided edges and runs as the source of truth.
fn check_lf(bwt: &BWT, edges: &[Vec<(usize, usize)>], runs: &[Vec<(usize, usize)>]) {
    // `lf()` at each offset of each record.
    for i in 0..bwt.len() {
        if let Some(record) = bwt.record(i) {
            let mut offset = 0;
            let mut curr_edges = edges[i].clone();
            let curr_runs = &runs[i];
            for (rank, len) in curr_runs {
                for _ in 0..*len {
                    let expected = if curr_edges[*rank].0 == ENDMARKER { None } else { Some(curr_edges[*rank]) };
                    assert_eq!(record.lf(offset), expected, "Invalid lf({}) in record {}", offset, i);
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
fn check_follow(bwt: &BWT) {
    for i in 0..bwt.len() {
        if let Some(record) = bwt.record(i) {
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
                    assert_eq!(record.follow(start..limit, 8), None, "Got a follow({}..{}, invalid) result in record {}", start, limit, i);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

#[test]
fn empty_bwt() {
    let edges = Vec::new();
    let runs = Vec::new();
    let bwt = create_bwt(&edges, &runs);
    check_edges(&bwt, &edges);
    check_lf(&bwt, &edges, &runs);
    check_follow(&bwt);
    serialize::test(&bwt, "empty-bwt", None, true);
}

#[test]
fn non_empty_bwt() {
    let edges = get_edges();
    let runs = get_runs();
    let bwt = create_bwt(&edges, &runs);
    check_edges(&bwt, &edges);
    check_lf(&bwt, &edges, &runs);
    check_follow(&bwt);
    serialize::test(&bwt, "non-empty-bwt", None, true);
}

#[test]
fn empty_records() {
    let mut edges = get_edges();
    edges[2] = Vec::new();
    edges[6] = Vec::new();

    let mut runs = get_runs();
    runs[2] = Vec::new();
    runs[6] = Vec::new();
 
    let bwt = create_bwt(&edges, &runs);
    check_edges(&bwt, &edges);
    check_lf(&bwt, &edges, &runs);
    check_follow(&bwt);
    serialize::test(&bwt, "bwt-with-empty", None, true);
}

//-----------------------------------------------------------------------------
