use super::*;

use crate::GBWT;
use crate::support;

use simple_sds::serialize;

use std::convert::TryFrom;

//-----------------------------------------------------------------------------

const SAMPLES: usize = 5;
const CONTIGS: usize = 4;
const PHASES: usize = 2;

fn create_metadata(paths: bool, samples: bool, contigs: bool) -> Metadata {
    let mut header = Header::<MetadataPayload>::new();
    header.payload_mut().sample_count = SAMPLES as u64;
    header.payload_mut().haplotype_count = (SAMPLES * PHASES) as u64;
    header.payload_mut().contig_count = CONTIGS as u64;

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

fn check_generic_path_names(metadata: &Metadata, expected_count: usize, test_case: &str) {
    let generic_sample_id = metadata.sample_id(GENERIC_SAMPLE);
    assert!(generic_sample_id.is_some(), "{}: Generic sample name {} is missing", test_case, GENERIC_SAMPLE);
    let generic_sample_id = generic_sample_id.unwrap();

    let mut generic_count = 0;
    for path_id in 0..metadata.paths() {
        let path_name = metadata.path(path_id).unwrap();
        if path_name.sample() == generic_sample_id {
            generic_count += 1;
            assert_eq!(path_name.phase(), 0, "{}: Generic path name {} has a nonzero phase", test_case, path_id);
        }
    }
    assert_eq!(generic_count, expected_count, "{}: Invalid number of generic path names", test_case);
}

#[test]
fn generic_path_names() {
    let filename = support::get_test_data("example.meta");
    let metadata: Metadata = serialize::load_from(&filename).unwrap();
    check_generic_path_names(&metadata, 2, "example.meta");

    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let metadata = index.metadata().unwrap();
    check_generic_path_names(metadata, 2, "example.gbwt");
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

        // And now in the other direction.
        let found_id = metadata.find_path(&from_metadata);
        assert_eq!(found_id, Some(path_id), "Failed to find path id for FullPathName from metadata for path {}", path_id);
    }
}

#[test]
fn full_path_name_generic() {
    let name = "example";
    let path_name = FullPathName::generic("example");

    let string_name = name;
    assert_eq!(&path_name.to_string(), string_name, "Wrong string representation");
    check_path_name(&path_name, GENERIC_SAMPLE, name, 0, 0, 123);
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

fn compare_full_metadata(metadata: &Metadata, truth: &Metadata) {
    let same_metadata = metadata == truth;
    if same_metadata {
        return;
    }

    let same_headers = metadata.header == truth.header;
    if !same_headers {
        assert_eq!(metadata.header.payload().sample_count, truth.header.payload().sample_count, "Sample counts differ");
        assert_eq!(metadata.header.payload().haplotype_count, truth.header.payload().haplotype_count, "Haplotype counts differ");
        assert_eq!(metadata.header.payload().contig_count, truth.header.payload().contig_count, "Contig counts differ");
    }
    assert!(same_headers, "Headers differ");

    let same_path_names = metadata.path_names == truth.path_names;
    if !same_path_names {
        assert_eq!(metadata.path_names.len(), truth.path_names.len(), "Path name counts differ");
        for (i, (path_name, truth_path_name)) in metadata.path_names.iter().zip(truth.path_names.iter()).enumerate() {
            assert_eq!(path_name.sample(), truth_path_name.sample(), "Sample ids differ for path name {}", i);
            assert_eq!(path_name.contig(), truth_path_name.contig(), "Contig ids differ for path name {}", i);
            assert_eq!(path_name.phase(), truth_path_name.phase(), "Haplotype numbers differ for path name {}", i);
            assert_eq!(path_name.fragment(), truth_path_name.fragment(), "Fragment numbers differ for path name {}", i);
        }
    }
    assert!(same_path_names, "Path names differ");

    let same_sample_names = metadata.sample_names == truth.sample_names;
    if !same_sample_names {
        assert_eq!(metadata.sample_names.len(), truth.sample_names.len(), "Sample name counts differ");
        for i in 0..metadata.sample_names.len() {
            let sample_name = metadata.sample_name(i);
            let truth_sample_name = truth.sample_name(i);
            assert_eq!(sample_name, truth_sample_name, "Sample names differ for sample name {}", i);
        }
    }
    assert!(same_sample_names, "Sample names differ");

    let same_contig_names = metadata.contig_names == truth.contig_names;
    if !same_contig_names {
        assert_eq!(metadata.contig_names.len(), truth.contig_names.len(), "Contig name counts differ");
        for i in 0..metadata.contig_names.len() {
            let contig_name = metadata.contig_name(i);
            let truth_contig_name = truth.contig_name(i);
            assert_eq!(contig_name, truth_contig_name, "Contig names differ for contig name {}", i);
        }
    }
    assert!(same_contig_names, "Contig names differ");
}

#[test]
fn metadata_builder_empty() {
    let empty = Metadata::from(MetadataBuilder::new());

    assert!(empty.has_path_names(), "Empty metadata should have path names");
    assert_eq!(empty.paths(), 0, "Empty metadata should have zero paths");

    assert!(empty.has_sample_names(), "Empty metadata should have sample names");
    assert_eq!(empty.samples(), 0, "Empty metadata should have zero samples");
    assert_eq!(empty.haplotypes(), 0, "Empty metadata should have zero haplotypes");

    assert!(empty.has_contig_names(), "Empty metadata should have contig names");
    assert_eq!(empty.contigs(), 0, "Empty metadata should have zero contigs");
}

#[test]
fn metadata_builder_insert() {
    let truth = create_metadata(true, true, true);

    let mut builder = MetadataBuilder::new();
    for path_id in 0..truth.paths() {
        let full_path_name = FullPathName::from_metadata(&truth, path_id).unwrap();
        let result = builder.insert(&full_path_name);
        assert!(result.is_ok(), "Failed to insert path name {} to builder: {}", path_id, result.err().unwrap());
    }

    // Check that we can get `FullPathName`s back from the builder.
    for path_id in 0..truth.paths() {
        let path_name = builder.path_name(path_id).unwrap();
        let true_name = FullPathName::from_metadata(&truth, path_id).unwrap();
        assert_eq!(path_name, true_name, "FullPathName differs for path {}", path_id);
    }
    assert!(builder.path_name(truth.paths()).is_none(), "Got a path name past the end of the builder");

    // Compare the result against the truth.
    let metadata = Metadata::from(builder);
    compare_full_metadata(&metadata, &truth);
}

#[test]
fn metadata_builder_insert_duplicates() {
    let truth = create_metadata(true, true, true);

    let mut builder = MetadataBuilder::new();
    for path_id in 0..truth.paths() {
        let full_path_name = FullPathName::from_metadata(&truth, path_id).unwrap();
        let result = builder.insert(&full_path_name);
        assert!(result.is_ok(), "Failed to insert path name {} to builder: {}", path_id, result.err().unwrap());
        let result = builder.insert(&full_path_name);
        assert!(result.is_err(), "Expected an error when inserting a duplicate path name {} to the builder", path_id);
    }

    // Compare the result against the truth.
    let metadata = Metadata::from(builder);
    compare_full_metadata(&metadata, &truth);
}

#[test]
fn metadata_builder_extend() {
    let truth = create_metadata(true, true, true);
    let path_names: Vec<FullPathName> = (0..truth.paths()).map(|path_id|
        FullPathName::from_metadata(&truth, path_id).unwrap()
    ).collect();

    let mut builder = MetadataBuilder::new();
    let result = builder.extend(&path_names);
    assert!(result.is_ok(), "Failed to extend builder with path names: {}", result.err().unwrap());

    // Compare the result against the truth.
    let metadata = Metadata::from(builder);
    compare_full_metadata(&metadata, &truth);
}

#[test]
fn metadata_builder_extend_twice() {
    let truth = create_metadata(true, true, true);
    let path_names: Vec<FullPathName> = (0..truth.paths()).map(|path_id|
        FullPathName::from_metadata(&truth, path_id).unwrap()
    ).collect();

    let mut builder = MetadataBuilder::new();
    let mid = path_names.len() / 2;
    let result = builder.extend(&path_names[..mid]);
    assert!(result.is_ok(), "Failed to extend builder with the first half of the path names: {}", result.err().unwrap());
    let result = builder.extend(&path_names[mid..]);
    assert!(result.is_ok(), "Failed to extend builder with the second half of the path names: {}", result.err().unwrap());

    // Compare the result against the truth.
    let metadata = Metadata::from(builder);
    compare_full_metadata(&metadata, &truth);
}

#[test]
fn metadata_builder_extend_duplicates() {
    let truth = create_metadata(true, true, true);
    let path_names: Vec<FullPathName> = (0..truth.paths()).map(|path_id|
        FullPathName::from_metadata(&truth, path_id).unwrap()
    ).collect();

    let mut builder = MetadataBuilder::new();
    let result = builder.extend(&path_names);
    assert!(result.is_ok(), "Failed to extend builder with path names: {}", result.err().unwrap());
    let result = builder.extend(&path_names);
    assert!(result.is_err(), "Expected an error when extending builder with duplicate path names");

    // Compare the result against the truth.
    let metadata = Metadata::from(builder);
    compare_full_metadata(&metadata, &truth);
}

//-----------------------------------------------------------------------------
