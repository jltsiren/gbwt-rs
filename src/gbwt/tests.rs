use super::*;

use simple_sds::serialize;

use std::collections::HashSet;
use std::convert::TryFrom;

//-----------------------------------------------------------------------------

#[test]
fn statistics() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();

    assert_eq!(index.len(), 68, "Invalid total length");
    assert!(!index.is_empty(), "Invalid emptiness");
    assert_eq!(index.sequences(), 12, "Invalid number of sequences");
    assert_eq!(index.alphabet_size(), 52, "Invalid alphabet size");
    assert_eq!(index.alphabet_offset(), 21, "Invalid alphabet offset");
    assert_eq!(index.effective_size(), 31, "Invalid effective alphabet size");
    assert_eq!(index.first_node(), 22, "Invalid first node id");
    assert!(index.is_bidirectional(), "Index is not bidirectional");

    for i in 0..index.first_node() {
        assert!(!index.has_node(i), "Index should not contain node {}", i);
    }
    for i in index.first_node()..index.alphabet_size() {
        assert!(index.has_node(i), "Index should contain node {}", i);
    }
    assert!(!index.has_node(index.alphabet_size()), "Index contains a node past the end");
}

#[test]
fn index_metadata() {
    let gbwt_filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&gbwt_filename).unwrap();
    assert!(index.has_metadata(), "Index does not contain metadata");

    let metadata_filename = support::get_test_data("example.meta");
    let metadata: Metadata = serialize::load_from(&metadata_filename).unwrap();
    assert_eq!(index.metadata().unwrap(), &metadata, "Invalid metadata in the index");
}

#[test]
fn serialize() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    serialize::test(&index, "gbwt", None, true);
}

//-----------------------------------------------------------------------------

fn extract_sequence(index: &GBWT, id: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut pos = index.start(id);
    while pos != None {
        result.push(pos.unwrap().0);
        pos = index.forward(pos.unwrap());
    }
    result
}

fn extract_backward(index: &GBWT, id: usize) -> Vec<usize> {
    let mut last = None;
    let mut pos = index.start(id);
    while pos != None {
        last = pos;
        pos = index.forward(pos.unwrap());
    }

    let mut result = Vec::new();
    pos = last;
    while pos != None {
        result.push(pos.unwrap().0);
        pos = index.backward(pos.unwrap());
    }

    result
}

fn true_paths() -> Vec<Vec<usize>> {
    vec![
        vec![support::encode_node(11, false), support::encode_node(12, false), support::encode_node(14, false), support::encode_node(15, false), support::encode_node(17, false)],
        vec![support::encode_node(21, false), support::encode_node(22, false), support::encode_node(24, false), support::encode_node(25, false)],
        vec![support::encode_node(11, false), support::encode_node(12, false), support::encode_node(14, false), support::encode_node(15, false), support::encode_node(17, false)],
        vec![support::encode_node(11, false), support::encode_node(13, false), support::encode_node(14, false), support::encode_node(16, false), support::encode_node(17, false)],
        vec![support::encode_node(21, false), support::encode_node(22, false), support::encode_node(24, false), support::encode_node(23, true), support::encode_node(21, true)],
        vec![support::encode_node(21, false), support::encode_node(22, false), support::encode_node(24, false), support::encode_node(25, false)],
    ]
}

#[test]
fn extract() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let truth = true_paths();

    for i in 0..index.sequences() / 2 {
        let forward = extract_sequence(&index, support::encode_node(i, false));
        assert_eq!(forward, truth[i], "Invalid forward path {}", i);
        let reverse = extract_sequence(&index, support::encode_node(i, true));
        assert_eq!(reverse.len(), forward.len(), "Invalid reverse path {} length", i);
        for j in 0..reverse.len() {
            let expected = support::flip_node(forward[forward.len() - j - 1]);
            assert_eq!(reverse[j], expected, "Invalid node {} on reverse path {}", j, i);
        }
    }
}

#[test]
fn backward() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();

    for i in 0..index.sequences() {
        let forward = extract_sequence(&index, i);
        let reverse = extract_backward(&index, i);
        assert_eq!(reverse.len(), forward.len(), "Invalid reverse sequence {} length", i);
        for j in 0..reverse.len() {
            let expected = forward[forward.len() - j - 1];
            assert_eq!(reverse[j], expected, "Invalid node {} on reverse sequence {}", j, i);
        }
    }
}

#[test]
fn sequence() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();

    for i in 0..index.sequences() {
        let extracted = extract_sequence(&index, i);
        let iterated: Vec<usize> = index.sequence(i).collect();
        assert_eq!(iterated, extracted, "Invalid sequence {} from an iterator", i);
    }
}

//-----------------------------------------------------------------------------

fn true_nodes() -> HashSet<usize> {
    let nodes: Vec<usize> = vec![11, 12, 13, 14, 15, 16, 17, 21, 22, 23, 24, 25];
    let mut result: HashSet<usize> = HashSet::new();
    for node in nodes.iter() {
        result.insert(support::encode_node(*node, false));
        result.insert(support::encode_node(*node, true));
    }
    result
}

fn count_occurrences(paths: &[Vec<usize>], subpath: &[usize]) -> usize {
    let mut result = 0;
    let reverse = support::reverse_path(subpath);
    for path in paths {
        for i in 0..path.len() {
            if path[i..].starts_with(subpath) {
                result += 1;
            }
            if path[..i + 1].ends_with(&reverse) {
                result += 1;
            }
        }
    }
    result
}

#[test]
fn find() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();

    for i in 0..index.alphabet_size() + 1 {
        if let Some(state) = index.find(i) {
            assert!(nodes.contains(&i), "Found a search state for a nonexistent node {}", i);
            assert_eq!(state.node, i, "Found an invalid search state for node {}", i);
            assert!(!state.is_empty(), "Found an empty search state for node {}", i);
        } else {
            assert!(!nodes.contains(&i), "Did not find a search state for node {}", i);
        }
    }
}

#[test]
fn extend() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();
    let paths = true_paths();

    // Check all possible and impossible extensions of the initial node.
    for &first in nodes.iter() {
        let start = index.find(first).unwrap();
        for i in 0..index.alphabet_size() + 1 {
            let count = count_occurrences(&paths, &[first, i]);
            if let Some(state) = index.extend(&start, i) {
                assert_eq!(state.len(), count, "Invalid number of occurrences for substring {} to {}", first, i);
            } else {
                assert_eq!(count, 0, "Could not find the occurrences of substring {} to {}", first, i);
            }
        }
    }

    // Search for all existing subpaths.
    for i in 0..paths.len() {
        let path = &paths[i];
        for j in 0..path.len() {
            let mut forward = index.find(path[j]).unwrap();
            for k in j + 1..path.len() {
                if let Some(state) = index.extend(&forward, path[k]) {
                    let count = count_occurrences(&paths, &path[j..k + 1]);
                    assert_eq!(state.len(), count, "Invalid number of occurrences for path {} at {}..{}", i, j, k + 1);
                    forward = state;
                } else {
                    panic!("Could not find occurrences of path {} at {}..{}", i, j, k + 1);
                }
            }

            let mut backward = index.find(support::flip_node(path[j])).unwrap();
            for k in (0..j).rev() {
                if let Some(state) = index.extend(&backward, support::flip_node(path[k])) {
                    let count = count_occurrences(&paths, &path[k..j + 1]); // No need to reverse the pattern here.
                    assert_eq!(state.len(), count, "Invalid number of occurrences for path {} at {}..{} (reversed)", i, k, j + 1);
                    backward = state;
                } else {
                    panic!("Could not find occurrences of path {} at {}..{} (reversed)", i, k, j + 1);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

fn bd_search(index: &GBWT, path: &[usize], first: usize, range: Range<usize>) -> Option<BidirectionalState> {
    let mut state = index.bd_find(path[first])?;
    for i in first + 1..range.end {
        state = index.extend_forward(&state, path[i])?;
    }
    for i in (range.start..first).rev() {
        state = index.extend_backward(&state, path[i])?;
    }
    Some(state)
}

#[test]
fn bd_find() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();

    for i in 0..index.alphabet_size() + 1 {
        if let Some(state) = index.bd_find(i) {
            assert!(nodes.contains(&i), "Found a search state for a nonexistent node {}", i);
            assert_eq!(state.forward.node, i, "Found an invalid search state for node {}", i);
            assert!(!state.is_empty(), "Found an empty search state for node {}", i);
            assert_eq!(state.reverse.node, support::flip_node(i), "Found an invalid reverse node for node {}", i);
            assert_eq!(state.reverse.len(), state.forward.len(), "Invalid reverse range length for node {}", i);
        } else {
            assert!(!nodes.contains(&i), "Did not find a search state for node {}", i);
        }
    }
}

#[test]
fn bd_extend() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();
    let paths = true_paths();

    // Check all possible and impossible extensions of the initial node.
    for &first in nodes.iter() {
        let start = index.bd_find(first).unwrap();
        for i in 0..index.alphabet_size() + 1 {
            // Forward.
            let count = count_occurrences(&paths, &[first, i]);
            if let Some(state) = index.extend_forward(&start, i) {
                assert_eq!(state.len(), count, "Invalid number of occurrences for substring {} to {} (forward)", first, i);
            } else {
                assert_eq!(count, 0, "Could not find the occurrences of substring {} to {} (forward)", first, i);
            }
            // Backward.
            let count = count_occurrences(&paths, &[i, first]);
            if let Some(state) = index.extend_backward(&start, i) {
                assert_eq!(state.len(), count, "Invalid number of occurrences for substring {} to {} (backward)", i, first);
            } else {
                assert_eq!(count, 0, "Could not find the occurrences of substring {} to {} (backward)", i, first);
            }
        }
    }

    // Search for all existing subpaths.
    for i in 0..paths.len() {
        let path = &paths[i];
        for first in 0..path.len() {
            for start in 0..first + 1 {
                for end in first + 1..path.len() + 1 {
                    // Forward path.
                    let count = count_occurrences(&paths, &path[start..end]);
                    if let Some(state) = bd_search(&index, &path, first, start..end) {
                        assert_eq!(state.len(), count, "Invalid number of occurrences for path {} at {}..{} from {}", i, start, end, first);
                        assert_eq!(state.reverse.len(), state.len(), "Invalid reverse state length for path {} at {}..{} from {}", i, start, end, first);
                        assert_eq!(state.forward.node, path[end - 1], "Invalid final node for path {} at {}..{} from {}", i, start, end, first);
                        assert_eq!(state.reverse.node, support::flip_node(path[start]), "Invalid initial node for path {} at {}..{} from {}", i, start, end, first);
                    } else {
                        panic!("Could not find occurrences of path {} at {}..{} from {}", i, start, end, first);
                    }

                    // Reverse path.
                    let reversed = support::reverse_path(&path);
                    let count = count_occurrences(&paths, &reversed[start..end]);
                    if let Some(state) = bd_search(&index, &reversed, first, start..end) {
                        assert_eq!(state.len(), count, "Invalid number of occurrences for path {} at {}..{} from {} (reversed)", i, start, end, first);
                        assert_eq!(state.reverse.len(), state.len(), "Invalid reverse state length for path {} at {}..{} from {} (reversed)", i, start, end, first);
                        assert_eq!(state.forward.node, reversed[end - 1], "Invalid final node for path {} at {}..{} from {} (reversed)", i, start, end, first);
                        assert_eq!(state.reverse.node, support::flip_node(reversed[start]), "Invalid initial node for path {} at {}..{} from {} (reversed)", i, start, end, first);
                    } else {
                        panic!("Could not find occurrences of path {} at {}..{} from {} (reversed)", i, start, end, first);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

const SAMPLES: usize = 5;
const CONTIGS: usize = 4;
const PHASES: usize = 2;

fn create_metadata(paths: bool, samples: bool, contigs: bool) -> Metadata {
    let mut header = Header::<MetadataPayload>::new();
    header.payload_mut().sample_count = SAMPLES;
    header.payload_mut().haplotype_count = SAMPLES * PHASES;
    header.payload_mut().contig_count = CONTIGS;

    let mut path_names = Vec::<PathName>::new();
    if paths {
        header.set(MetadataPayload::FLAG_PATH_NAMES);
        for sample in 0..SAMPLES {
            for contig in 0..CONTIGS {
                for phase in 0..PHASES {
                    path_names.push(PathName::from_fields(sample, contig, phase, 0));
                }
            }
        }
    }

    let mut sample_names = Vec::<String>::new();
    if samples {
        header.set(MetadataPayload::FLAG_SAMPLE_NAMES);
        for sample in 0..SAMPLES {
            sample_names.push(format!("sample_{}", sample));
        }
    }

    let mut contig_names = Vec::<String>::new();
    if contigs {
        header.set(MetadataPayload::FLAG_CONTIG_NAMES);
        for contig in 0..CONTIGS {
            contig_names.push(format!("contig_{}", contig));
        }
    }

    Metadata {
        header: header,
        path_names: path_names,
        sample_names: Dictionary::try_from(sample_names).unwrap(),
        contig_names: Dictionary::try_from(contig_names).unwrap(),
    }
}

fn test_metadata(paths: bool, samples: bool, contigs: bool, name: &str) {
    let metadata = create_metadata(paths, samples, contigs);

    // Contents.
    assert_eq!(metadata.has_path_names(), paths, "{}: Invalid path name flag", name);
    assert_eq!(metadata.has_sample_names(), samples, "{}: Invalid sample name flag", name);
    assert_eq!(metadata.has_contig_names(), contigs, "{}: Invalid contig name flag", name);

    // Statistics.
    if paths {
        assert_eq!(metadata.paths(), SAMPLES * CONTIGS * PHASES, "{}: Invalid path count", name);
    } else {
        assert_eq!(metadata.paths(), 0, "{}: Invalid path count", name);
    }
    assert_eq!(metadata.samples(), SAMPLES, "{}: Invalid sample count", name);
    assert_eq!(metadata.haplotypes(), SAMPLES * PHASES, "{}: Invalid haplotype count", name);
    assert_eq!(metadata.contigs(), CONTIGS, "{}: Invalid contig count", name);

    // Path names.
    if paths {
        let mut index = 0;
        let mut iter = metadata.path_iter();
        for sample in 0..SAMPLES {
            for contig in 0..CONTIGS {
                for phase in 0..PHASES {
                    let path = PathName::from_fields(sample, contig, phase, 0);
                    assert_eq!(metadata.path(index), Some(path), "{}: Invalid path name {}", name, index);
                    assert_eq!(iter.next(), Some(&path), "{}: Invalid path name {} from iterator", name, index);
                    index += 1;
                }
            }
        }
        assert_eq!(metadata.path(index), None, "{}: Got a path name past the end", name);
        assert_eq!(iter.next(), None, "{}: Got a path name past the end from iterator", name);
    }

    // Sample names.
    if samples {
        let mut iter = metadata.sample_iter();
        for sample in 0..SAMPLES {
            let sample_name = format!("sample_{}", sample);
            assert_eq!(metadata.sample(sample), Some(sample_name.as_str()), "{}: Invalid sample name {}", name, sample);
            assert_eq!(iter.next(), Some(sample_name.as_bytes()), "{}: Invalid sample name {} from iterator", name, sample);
            assert_eq!(metadata.sample_id(&sample_name), Some(sample), "{}: Invalid id for sample {}", name, sample_name);
        }
        assert_eq!(metadata.sample(SAMPLES), None, "{}: Got a sample name past the end", name);
        assert_eq!(iter.next(), None, "{}: Got a sample name past the end from iterator", name);
    }
    assert!(metadata.sample_id("invalid").is_none(), "{}: Got an id for an invalid sample", name);

    // Contig names.
    if contigs {
        let mut iter = metadata.contig_iter();
        for contig in 0..CONTIGS {
            let contig_name = format!("contig_{}", contig);
            assert_eq!(metadata.contig(contig), Some(contig_name.as_str()), "{}: Invalid contig name {}", name, contig);
            assert_eq!(iter.next(), Some(contig_name.as_bytes()), "{}: Invalid contig name {} from iterator", name, contig);
            assert_eq!(metadata.contig_id(&contig_name), Some(contig), "{}: Invalid id for contig {}", name, contig_name);
        }
        assert_eq!(metadata.contig(CONTIGS), None, "{}: Got a contig name past the end", name);
        assert_eq!(iter.next(), None, "{}: Got a contig name past the end from iterator", name);
    }
    assert!(metadata.contig_id("invalid").is_none(), "{}: Got an id for an invalid contig", name);

    serialize::test(&metadata, name, None, true);
}

#[test]
fn metadata() {
    test_metadata(true, false, false, "Paths");
    test_metadata(false, true, false, "Samples");
    test_metadata(false, false, true, "Contigs");
}

#[test]
fn path_names() {
    let name = PathName::new();
    assert_eq!(name.size_in_elements(), 2, "Invalid serialized size for a path name");
}

//-----------------------------------------------------------------------------
