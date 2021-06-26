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

fn check_bwt(bwt: &BWT, edges: &[Vec<(usize, usize)>], _: &[Vec<(usize, usize)>]) {
    assert_eq!(bwt.len(), edges.len(), "Invalid number of records in the BWT");
    assert_eq!(bwt.is_empty(), edges.is_empty(), "Invalid BWT emptiness");

    for i in 0..bwt.len() {
        let record = bwt.record(i);
        let curr_edges = &edges[i];
        assert_eq!(record.is_none(), curr_edges.is_empty(), "Invalid record {} existence", i);
        if let Some(record) = record {
            assert_eq!(record.outdegree(), curr_edges.len(), "Invalid outdegree for record {}", i);
            for j in 0..record.outdegree() {
                assert_eq!(record.successor(j), curr_edges[j].0, "Invalid successor {} for record {}", j, i);
                assert_eq!(record.offset(j), curr_edges[j].1, "Invalid offset {} for record {}", j, i);
            }
        }
    }

    // TODO: Check runs once we have the functionality.
}

//-----------------------------------------------------------------------------

#[test]
fn empty_bwt() {
    let edges = Vec::new();
    let runs = Vec::new();
    let bwt = create_bwt(&edges, &runs);
    check_bwt(&bwt, &edges, &runs);
    serialize::test(&bwt, "empty-bwt", None, true);
}

#[test]
fn non_empty_bwt() {
    let edges = get_edges();
    let runs = get_runs();
    let bwt = create_bwt(&edges, &runs);
    check_bwt(&bwt, &edges, &runs);
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
    check_bwt(&bwt, &edges, &runs);
    serialize::test(&bwt, "bwt-with-empty", None, true);
}

//-----------------------------------------------------------------------------
