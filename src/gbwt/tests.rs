use super::*;

use crate::Orientation;

use simple_sds::serialize;

use std::collections::HashSet;
use std::convert::TryFrom;

//-----------------------------------------------------------------------------

// Assumes a non-empty bidirectional GBWT.
fn check_statistics(index: &GBWT, len: usize, sequences: usize, alphabet_size: usize, alphabet_offset: usize) {
    assert_eq!(index.len(), len, "Invalid total length");
    assert!(!index.is_empty(), "Invalid emptiness");
    assert_eq!(index.sequences(), sequences, "Invalid number of sequences");
    assert_eq!(index.alphabet_size(), alphabet_size, "Invalid alphabet size");
    assert_eq!(index.alphabet_offset(), alphabet_offset, "Invalid alphabet offset");
    assert_eq!(index.effective_size(), alphabet_size - alphabet_offset, "Invalid effective alphabet size");
    assert_eq!(index.first_node(), alphabet_offset + 1, "Invalid first node id");
    assert!(index.is_bidirectional(), "Index is not bidirectional");

    for i in 0..index.first_node() {
        assert!(!index.has_node(i), "Index should not contain node {}", i);
    }
    for i in index.first_node()..index.alphabet_size() {
        assert!(index.has_node(i), "Index should contain node {}", i);
    }
    assert!(!index.has_node(index.alphabet_size()), "Index contains a node past the end");

    // Node / record id conversions.
    for node_id in index.first_node()..index.alphabet_size() {
        let record_id = index.node_to_record(node_id);
        let converted = index.record_to_node(record_id);
        assert_eq!(converted, node_id, "Node -> record -> node conversion failed for {}", node_id);
    }
    for record_id in 1..index.effective_size() {
        let node_id = index.record_to_node(record_id);
        let converted = index.node_to_record(node_id);
        assert_eq!(converted, record_id, "Record -> node -> record conversion failed for {}", record_id);
    }
}

#[test]
fn statistics() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    check_statistics(&index, 68, 12, 52, 21);
}

#[test]
fn statistics_with_empty() {
    let filename = support::get_test_data("with-empty.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    check_statistics(&index, 70, 14, 52, 21);
}

#[test]
fn index_metadata() {
    let gbwt_filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&gbwt_filename).unwrap();
    assert!(index.has_metadata(), "Index does not contain metadata");

    let metadata_filename = support::get_test_data("example.meta");
    let metadata: Metadata = serialize::load_from(&metadata_filename).unwrap();
    assert_eq!(index.metadata().unwrap(), &metadata, "Invalid metadata in the index");

    let tags = index.tags();
    assert!(tags.contains_key(SOURCE_KEY), "Mandatory source key {} is missing", SOURCE_KEY);
}

#[test]
fn serialize() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    serialize::test(&index, "gbwt", None, true);
}

#[test]
fn serialize_with_empty() {
    let filename = support::get_test_data("with-empty.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    serialize::test(&index, "gbwt-with-empty", None, true);
}

//-----------------------------------------------------------------------------

fn extract_sequence(index: &GBWT, id: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut pos = index.start(id);
    while let Some(p) = pos {
        result.push(p.node);
        pos = index.forward(p);
    }
    result
}

fn extract_backward(index: &GBWT, id: usize) -> Vec<usize> {
    let mut last = None;
    let mut pos = index.start(id);
    while let Some(p) = pos {
        last = pos;
        pos = index.forward(p);
    }

    let mut result = Vec::new();
    pos = last;
    while let Some(p) = pos {
        result.push(p.node);
        pos = index.backward(p);
    }

    result
}

fn true_paths(with_empty: bool) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    result.push(vec![
        support::encode_node(11, Orientation::Forward),
        support::encode_node(12, Orientation::Forward),
        support::encode_node(14, Orientation::Forward),
        support::encode_node(15, Orientation::Forward),
        support::encode_node(17, Orientation::Forward)
    ]);
    result.push(vec![
        support::encode_node(21, Orientation::Forward),
        support::encode_node(22, Orientation::Forward),
        support::encode_node(24, Orientation::Forward),
        support::encode_node(25, Orientation::Forward)
    ]);
    result.push(vec![
        support::encode_node(11, Orientation::Forward),
        support::encode_node(12, Orientation::Forward),
        support::encode_node(14, Orientation::Forward),
        support::encode_node(15, Orientation::Forward),
        support::encode_node(17, Orientation::Forward)
    ]);
    result.push(vec![
        support::encode_node(11, Orientation::Forward),
        support::encode_node(13, Orientation::Forward),
        support::encode_node(14, Orientation::Forward),
        support::encode_node(16, Orientation::Forward),
        support::encode_node(17, Orientation::Forward)
    ]);
    if with_empty {
        result.push(Vec::new());
    }
    result.push(vec![
        support::encode_node(21, Orientation::Forward),
        support::encode_node(22, Orientation::Forward),
        support::encode_node(24, Orientation::Forward),
        support::encode_node(23, Orientation::Reverse),
        support::encode_node(21, Orientation::Reverse)
    ]);
    result.push(vec![
        support::encode_node(21, Orientation::Forward),
        support::encode_node(22, Orientation::Forward),
        support::encode_node(24, Orientation::Forward),
        support::encode_node(25, Orientation::Forward)
    ]);
    result
}

fn test_extract(test_file: &'static str, with_empty: bool) {
    let filename = support::get_test_data(test_file);
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let truth = true_paths(with_empty);

    for i in 0..index.sequences() / 2 {
        let forward = extract_sequence(&index, support::encode_node(i, Orientation::Forward));
        assert_eq!(forward, truth[i], "Invalid forward path {}", i);
        let reverse = extract_sequence(&index, support::encode_node(i, Orientation::Reverse));
        assert_eq!(reverse.len(), forward.len(), "Invalid reverse path {} length", i);
        for j in 0..reverse.len() {
            let expected = support::flip_node(forward[forward.len() - j - 1]);
            assert_eq!(reverse[j], expected, "Invalid node {} on reverse path {}", j, i);
        }
    }
}

#[test]
fn extract() {
    test_extract("example.gbwt", false);
}

#[test]
fn extract_with_empty() {
    test_extract("with-empty.gbwt", true);
}

fn test_backward(test_file: &'static str) {
    let filename = support::get_test_data(test_file);
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
fn backward() {
    test_backward("example.gbwt");
}

#[test]
fn backward_with_empty() {
    test_backward("with-empty.gbwt");
}

fn test_sequence(test_file: &'static str) {
    let filename = support::get_test_data(test_file);
    let index: GBWT = serialize::load_from(&filename).unwrap();

    for i in 0..index.sequences() {
        let extracted = extract_sequence(&index, i);
        let iter = index.sequence(i);
        assert!(iter.is_some(), "Could not get an iterator for sequence {}", i);
        let iterated: Vec<usize> = iter.unwrap().collect();
        assert_eq!(iterated, extracted, "Invalid sequence {} from an iterator", i);
    }
    assert!(index.sequence(index.sequences()).is_none(), "Got an iterator for a past-the-end sequence id");
}

#[test]
fn sequence() {
    test_sequence("example.gbwt");
}

#[test]
fn sequence_with_empty() {
    test_sequence("with-empty.gbwt");
}

//-----------------------------------------------------------------------------

fn true_nodes() -> HashSet<usize> {
    let nodes: Vec<usize> = vec![11, 12, 13, 14, 15, 16, 17, 21, 22, 23, 24, 25];
    let mut result: HashSet<usize> = HashSet::new();
    for node in nodes.iter() {
        result.insert(support::encode_node(*node, Orientation::Forward));
        result.insert(support::encode_node(*node, Orientation::Reverse));
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

fn test_find(test_file: &'static str) {
    let filename = support::get_test_data(test_file);
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
fn find() {
    test_find("example.gbwt");
}

#[test]
fn find_with_empty() {
    test_find("with-empty.gbwt");
}

fn test_extend(test_file: &'static str, with_empty: bool) {
    let filename = support::get_test_data(test_file);
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();
    let paths = true_paths(with_empty);

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

#[test]
fn extend() {
    test_extend("example.gbwt", false);
}

#[test]
fn extend_with_empty() {
    test_extend("with-empty.gbwt", true);
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

fn test_bd_find(test_file: &'static str) {
    let filename = support::get_test_data(test_file);
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
fn bd_find() {
    test_bd_find("example.gbwt");
}

#[test]
fn bd_find_with_empty() {
    test_bd_find("with-empty.gbwt");
}

fn test_bd_extend(test_file: &'static str, with_empty: bool) {
    let filename = support::get_test_data(test_file);
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let nodes = true_nodes();
    let paths = true_paths(with_empty);

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

#[test]
fn bd_extend() {
    test_bd_extend("example.gbwt", false);
}

#[test]
fn bd_extend_with_empty() {
    test_bd_extend("with-empty.gbwt", true);
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
                    let pan_sn = format!("{}#{}#{}", metadata.sample_name(sample), phase, metadata.contig_name(contig));
                    assert_eq!(metadata.path(index), Some(path), "{}: Invalid path name {}", name, index);
                    assert_eq!(metadata.pan_sn_path(index), Some(pan_sn), "{}: Invalid PanSN path name {}", name, index);
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
            assert_eq!(metadata.sample_name(sample), sample_name, "{}: Invalid forced sample name {}", name, sample);
            assert_eq!(iter.next(), Some(sample_name.as_bytes()), "{}: Invalid sample name {} from iterator", name, sample);
            assert_eq!(metadata.sample_id(&sample_name), Some(sample), "{}: Invalid id for sample {}", name, sample_name);
        }
        assert_eq!(metadata.sample(SAMPLES), None, "{}: Got a sample name past the end", name);
        assert_eq!(iter.next(), None, "{}: Got a sample name past the end from iterator", name);
    } else {
        for sample in 0..SAMPLES {
            let sample_name = sample.to_string();
            assert_eq!(metadata.sample_name(sample), sample_name, "{}: Invalid forced sample name {}", name, sample);
        }
    }
    assert!(metadata.sample_id("invalid").is_none(), "{}: Got an id for an invalid sample", name);

    // Contig names.
    if contigs {
        let mut iter = metadata.contig_iter();
        for contig in 0..CONTIGS {
            let contig_name = format!("contig_{}", contig);
            assert_eq!(metadata.contig(contig), Some(contig_name.as_str()), "{}: Invalid contig name {}", name, contig);
            assert_eq!(metadata.contig_name(contig), contig_name, "{}: Invalid forced contig name {}", name, contig);
            assert_eq!(iter.next(), Some(contig_name.as_bytes()), "{}: Invalid contig name {} from iterator", name, contig);
            assert_eq!(metadata.contig_id(&contig_name), Some(contig), "{}: Invalid id for contig {}", name, contig_name);
        }
        assert_eq!(metadata.contig(CONTIGS), None, "{}: Got a contig name past the end", name);
        assert_eq!(iter.next(), None, "{}: Got a contig name past the end from iterator", name);
    } else {
        for contig in 0..CONTIGS {
            let contig_name = contig.to_string();
            assert_eq!(metadata.contig_name(contig), contig_name, "{}: Invalid forced contig name {}", name, contig);
        }
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

fn check_path_name(path_name: &FullPathName, sample: &str, contig: &str, haplotype: usize, fragment: usize, len: usize) {
    assert_eq!(path_name.sample, sample, "Wrong sample name");
    assert_eq!(path_name.contig, contig, "Wrong contig name");
    assert_eq!(path_name.haplotype, haplotype, "Wrong haplotype number");
    assert_eq!(path_name.fragment, fragment, "Wrong fragment number");

    let pan_sn_name = format!("{}#{}#{}", sample, haplotype, contig);
    assert_eq!(path_name.pan_sn_name(), pan_sn_name, "Wrong PanSN name");

    let path_fragment_name = format!("{}#{}#{}[{}-{}]", sample, haplotype, contig, fragment, len);
    assert_eq!(path_name.path_fragment_name(len), path_fragment_name, "Wrong path fragment name");
}

#[test]
fn full_path_name_from_metadata() {
    let filename = support::get_test_data("example.meta");
    let metadata: Metadata = serialize::load_from(&filename).unwrap();

    for (path_id, path_name) in metadata.path_iter().enumerate() {
        let from_metadata = FullPathName::from_metadata(&metadata, path_id);
        assert!(from_metadata.is_some(), "Failed to create FullPathName from metadata for path {}", path_id);
        let from_metadata = from_metadata.unwrap();
        let truth = FullPathName {
            sample: metadata.sample_name(path_name.sample()),
            contig: metadata.contig_name(path_name.contig()),
            haplotype: path_name.phase(),
            fragment: path_name.fragment()
        };
        assert_eq!(from_metadata, truth, "Wrong FullPathName from metadata for path {}", path_id);
    }
}

#[test]
fn full_path_name_generic() {
    let name = "example";
    let path_name = FullPathName::generic("example");

    let string_name = name;
    assert_eq!(&path_name.to_string(), string_name, "Wrong string representation");
    check_path_name(&path_name, REF_SAMPLE, name, 0, 0, 123);
}

#[test]
fn full_path_name_reference() {
    let sample = "GRCh38";
    let contig = "chr1";
    let path_name = FullPathName::reference(sample, contig);

    let string_name = format!("{}#{}", sample, contig);
    assert_eq!(path_name.to_string(), string_name, "Wrong string representation");
    check_path_name(&path_name, sample, contig, 0, 0, 123);
}

#[test]
fn full_path_name_haplotype() {
    let sample = "NA12878";
    let contig = "chr1";
    let haplotype = 1;
    let fragment = 38000;
    let path_name = FullPathName::haplotype(sample, contig, haplotype, fragment);

    let string_name = format!("{}#{}#{}@{}", sample, haplotype, contig, fragment);
    assert_eq!(path_name.to_string(), string_name, "Wrong string representation");
    check_path_name(&path_name, sample, contig, haplotype, fragment, 4800);
}

//-----------------------------------------------------------------------------
