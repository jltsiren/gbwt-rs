//! Metadata for GBWT paths.
//!
//! Path metadata is stored as an [`Metadata`] structure.
//! If path names are present, they are stored as an ordered set of [`PathName`] objects.
//! Sample and contig fields in path names may have corresponding string names in the metadata.
//! [`FullPathName`] objects store sample and contig names in the path name itself.
//! [`MetadataBuilder`] can be used for creating metadata from [`FullPathName`] objects.

use crate::{GENERIC_SAMPLE, GENERIC_HAPLOTYPE};
use crate::headers::{Header, MetadataPayload};
use crate::support::{Dictionary, StringIter};

use simple_sds::serialize::{Serialize, Serializable};

use std::collections::{HashMap, HashSet};
use std::io::{Error, ErrorKind};
use std::{fmt, io, slice};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Metadata for the paths in a GBWT index.
///
/// The metadata contains some basic statistics about the paths, and it may also contain names for paths, samples, and contigs.
/// Samples correspond to biological samples.
/// Contigs usually correspond to contigs in the reference sequence as well as to connected components in the graph.
/// The number of haplotypes is an estimate for the number of full-length paths in a connected component.
///
/// Path names are structured names associated with each path in the original graph.
/// If the GBWT index is bidirectional, it contains two sequences for each path name.
/// Each path name is a unique combination of four fields: sample, contig, phase, and fragment (see [`PathName`]).
/// The first two must be in the intervals `0..self.samples()` and `0..self.contigs()`, respectively.
///
/// Metadata objects with path, sample, and contig names can be created using [`MetadataBuilder`].
///
/// # Examples
///
/// ```
/// use gbz::{Metadata, MetadataBuilder, FullPathName};
///
/// // This is the `example.meta` test case.
/// let mut builder = MetadataBuilder::new();
/// let _ = builder.insert(&FullPathName::generic("A"));
/// let _ = builder.insert(&FullPathName::generic("B"));
/// let _ = builder.insert(&FullPathName::haplotype("sample", "A", 1, 0));
/// let _ = builder.insert(&FullPathName::haplotype("sample", "A", 2, 0));
/// let _ = builder.insert(&FullPathName::haplotype("sample", "B", 1, 0));
/// let _ = builder.insert(&FullPathName::haplotype("sample", "B", 2, 0));
/// let metadata = Metadata::from(builder);
///
/// assert!(metadata.has_path_names());
/// assert_eq!(metadata.paths(), 6);
/// assert_eq!(metadata.pan_sn_path(3), Some("sample#2#A".to_string()));
/// let path = metadata.path(3).unwrap();
///
/// assert!(metadata.has_sample_names());
/// assert_eq!(metadata.sample(path.sample()), Some("sample"));
///
/// assert!(metadata.has_contig_names());
/// assert_eq!(metadata.contig(path.contig()), Some("A"));
///
/// // Find the paths over contig B.
/// let id = metadata.contig_id("B").unwrap();
/// let mut paths = Vec::<usize>::new();
/// for (i, path) in metadata.path_iter().enumerate() {
///     if path.contig() == id {
///         paths.push(i);
///     }
/// }
/// assert_eq!(paths, vec![1, 4, 5]);
///
/// // Find a path identifier by metadata.
/// let path_name = FullPathName::haplotype("sample", "A", 2, 0);
/// let path_id = metadata.find_path(&path_name);
/// assert_eq!(path_id, Some(3));
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Metadata {
    header: Header<MetadataPayload>,
    path_names: Vec<PathName>,
    sample_names: Dictionary,
    contig_names: Dictionary,
}

/// Paths.
impl Metadata {
    /// Returns `true` if the metadata contains path names.
    pub fn has_path_names(&self) -> bool {
        self.header.is_set(MetadataPayload::FLAG_PATH_NAMES)
    }

    /// Returns the number path names in the metadata.
    ///
    /// If there are path names, each name corresponds to a path in the original graph.
    /// In a bidirectional GBWT, there are twice as many sequences as paths.
    pub fn paths(&self) -> usize {
        self.path_names.len()
    }

    /// Returns the name of the path with the given identifier, or [`None`] if there is no such name.
    ///
    /// Valid identifiers are in the interval `0..self.paths()`.
    pub fn path(&self, id: usize) -> Option<PathName> {
        if id < self.paths() {
            Some(self.path_names[id])
        } else {
            None
        }
    }

    /// Returns the identifier of the path with the given name, or [`None`] if there is no such name.
    pub fn find_path(&self, name: &FullPathName) -> Option<usize> {
        if !self.has_path_names() || !self.has_sample_names() || !self.has_contig_names() {
            return None;
        }
        let path_name = PathName::from_fields(
            self.sample_id(&name.sample)?,
            self.contig_id(&name.contig)?,
            name.haplotype,
            name.fragment,
        );

        for (i, path) in self.path_names.iter().enumerate() {
            if *path == path_name {
                return Some(i);
            }
        }
        None
    }

    // TODO: We need a graph with fragmented haplotypes for testing this.
    /// Returns the identifier of the last haplotype fragment starting at or before the specified position.
    ///
    /// The fragment field is used for specifying an offset on the haplotype.
    /// The other fields are used for identifying the haplotype.
    /// Returns [`None`] if there is no such path.
    pub fn find_fragment(&self, name: &FullPathName) -> Option<usize> {
        if !self.has_path_names() || !self.has_sample_names() || !self.has_contig_names() {
            return None;
        }
        let path_name = PathName::from_fields(
            self.sample_id(&name.sample)?,
            self.contig_id(&name.contig)?,
            name.haplotype,
            name.fragment,
        );

        let mut prev_fragment = 0;
        let mut result: Option<usize> = None;
        for (i, path) in self.path_names.iter().enumerate() {
            if path.is_predecessor_of(&path_name)
            && (result.is_none() || (result.is_some() && path.fragment > prev_fragment)) {
                prev_fragment = path.fragment;
                result = Some(i);
            }
        }
        result
    }

    /// Returns the name of the path with the given identifier in the [PanSN format](https://github.com/pangenome/PanSN-spec), or [`None`] if there is no such name.
    ///
    /// Valid identifiers are in the interval `0..self.paths()`.
    /// `'#'` will be used as the separator character.
    pub fn pan_sn_path(&self, id: usize) -> Option<String> {
        let path_name = self.path(id)?;
        let result = format!("{}#{}#{}", self.sample_name(path_name.sample()), path_name.phase(), self.contig_name(path_name.contig()));
        Some(result)
    }

    /// Returns an iterator over path names.
    pub fn path_iter(&'_ self) -> slice::Iter<'_, PathName> {
        self.path_names.iter()
    }
}

/// Samples.
impl Metadata {
    /// Returns `true` if the metadata contains sample names.
    pub fn has_sample_names(&self) -> bool {
        self.header.is_set(MetadataPayload::FLAG_SAMPLE_NAMES)
    }

    /// Returns the number samples.
    pub fn samples(&self) -> usize {
        self.header.payload().sample_count as usize
    }

    /// Returns the number of haplotypes.
    ///
    /// This generally corresponds to the number of full-length paths in a graph component.
    pub fn haplotypes(&self) -> usize {
        self.header.payload().haplotype_count as usize
    }

    /// Returns the name of the sample with the given identifier, or [`None`] if there is no such sample or name.
    ///
    /// Valid identifiers are in the interval `0..self.samples()`.
    /// Also returns [`None`] if the name exists but is not valid UTF-8.
    pub fn sample(&self, id: usize) -> Option<&str> {
        if self.has_sample_names() && id < self.samples() {
            self.sample_names.str(id).ok()
        } else {
            None
        }
    }

    /// Returns the name of the sample with the given identifier.
    ///
    /// Returns a string representation of the the sample identifier when [`Metadata::sample`] would return [`None`].
    pub fn sample_name(&self, id: usize) -> String {
        if let Some(name) = self.sample(id) {
            name.to_string()
        } else {
            id.to_string()
        }
    }
 
    /// Returns the sample idenfier corresponding to the given sample name, or [`None`] if there is no such sample.
    pub fn sample_id(&self, name: &str) -> Option<usize> {
        self.sample_names.id(name)
    }

    /// Returns an iterator over sample names.
    pub fn sample_iter(&'_ self) -> StringIter<'_> {
        self.sample_names.as_ref().iter()
    }
}

/// Contigs.
impl Metadata {
    /// Returns `true` if the metadata contains contig names.
    pub fn has_contig_names(&self) -> bool {
        self.header.is_set(MetadataPayload::FLAG_CONTIG_NAMES)
    }

    /// Returns the number contigs.
    ///
    /// A contig usually corresponds to a graph component.
    pub fn contigs(&self) -> usize {
        self.header.payload().contig_count as usize
    }

    /// Returns the name of the contig with the given identifier, or [`None`] if there is no such contig or name.
    ///
    /// Valid identifiers are in the interval `0..self.contigs()`.
    /// Also returns [`None`] if the name exists but is not valid UTF-8.
    pub fn contig(&self, id: usize) -> Option<&str> {
        if self.has_contig_names() && id < self.contigs() {
            self.contig_names.str(id).ok()
        } else {
            None
        }
    }

    /// Returns the name of the contig with the given identifier.
    ///
    /// Returns a string representation of the the contig identifier when [`Metadata::contig`] would return [`None`].
    pub fn contig_name(&self, id: usize) -> String {
        if let Some(name) = self.contig(id) {
            name.to_string()
        } else {
            id.to_string()
        }
    }

    /// Returns the contig idenfier corresponding to the given contig name, or [`None`] if there is no such contig.
    pub fn contig_id(&self, name: &str) -> Option<usize> {
        self.contig_names.id(name)
    }

    /// Returns an iterator over contig names.
    pub fn contig_iter(&'_ self) -> StringIter<'_> {
        self.contig_names.as_ref().iter()
    }
}

impl Serialize for Metadata {
    fn serialize_header<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        self.header.serialize(writer)
    }

    fn serialize_body<T: io::Write>(&self, writer: &mut T) -> io::Result<()> {
        // Follow vg conventions with haplotype numbers for generic paths.
        if let Some(id) = self.sample_id(GENERIC_SAMPLE) {
            let mut path_names = self.path_names.clone();
            for path in path_names.iter_mut() {
                if path.sample() == id && path.phase == 0 {
                    path.phase = GENERIC_HAPLOTYPE;
                }
            }
            path_names.serialize(writer)?;
        } else {
            self.path_names.serialize(writer)?;
        }

        self.sample_names.serialize(writer)?;
        self.contig_names.serialize(writer)?;
        Ok(())
    }

    fn load<T: io::Read>(reader: &mut T) -> io::Result<Self> {
        let mut header = Header::<MetadataPayload>::load(reader)?;
        if let Err(msg) = header.validate() {
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }

        let mut path_names = Vec::<PathName>::load(reader)?;
        if header.is_set(MetadataPayload::FLAG_PATH_NAMES) == path_names.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Metadata: Path name flag does not match the presence of path names"));
        }

        let sample_names = Dictionary::load(reader)?;
        if header.is_set(MetadataPayload::FLAG_SAMPLE_NAMES) {
            if header.payload().sample_count as usize != sample_names.len() {
                return Err(Error::new(ErrorKind::InvalidData, "Metadata: Sample count does not match the number of sample names"));
            }
        } else if !sample_names.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Metadata: Sample names are present without the sample name flag"));
        }

        let contig_names = Dictionary::load(reader)?;
        if header.is_set(MetadataPayload::FLAG_CONTIG_NAMES) {
            if header.payload().contig_count as usize != contig_names.len() {
                return Err(Error::new(ErrorKind::InvalidData, "Metadata: Contig count does not match the number of contig names"));
            }
        } else if !contig_names.is_empty() {
            return Err(Error::new(ErrorKind::InvalidData, "Metadata: Contig names are present without the contig name flag"));
        }

        // Update the header to the latest version after we have used the
        // serialized version for loading the correct data.
        header.update();

        // Convert path names for generic paths from vg convention to ours.
        if let Some(id) = sample_names.id(GENERIC_SAMPLE) {
            for path in path_names.iter_mut() {
                if path.sample() == id && path.phase == GENERIC_HAPLOTYPE {
                    path.phase = 0;
                }
            }
        }

        Ok(Metadata {
            header, path_names, sample_names, contig_names,
        })
    }

    fn size_in_elements(&self) -> usize {
        self.header.size_in_elements() + self.path_names.size_in_elements() + self.sample_names.size_in_elements() + self.contig_names.size_in_elements()
    }
}

impl From<MetadataBuilder> for Metadata {
    fn from(builder: MetadataBuilder) -> Self {
        let mut header = Header::<MetadataPayload>::default();
        header.payload_mut().sample_count = builder.sample_names.len() as u64;
        header.payload_mut().contig_count = builder.contig_names.len() as u64;
        header.payload_mut().haplotype_count = builder.haplotypes.len() as u64;
        header.set(MetadataPayload::FLAG_PATH_NAMES);
        header.set(MetadataPayload::FLAG_SAMPLE_NAMES);
        header.set(MetadataPayload::FLAG_CONTIG_NAMES);

        // The panics in sample / contig name conversions should be unreachable.
        // The conversion can only fail if there are duplicate names, but the builder should prevent that.
        let path_names = builder.path_names;
        let sample_names = Dictionary::try_from(builder.sample_names).unwrap_or_else(|_|
            panic!("Metadata: Inconsistent sample names in MetadataBuilder")
        );
        let contig_names = Dictionary::try_from(builder.contig_names).unwrap_or_else(|_|
            panic!("Metadata: Inconsistent contig names in MetadataBuilder")
        );

        Metadata {
            header, path_names, sample_names, contig_names,
        }
    }
}

//-----------------------------------------------------------------------------

/// A structured path name relative to [`Metadata`].
///
/// Each path name in the same metadata structure must be unique.
/// A path name has four components: sample, contig, phase (haplotype), and fragment (count).
///
/// * Samples correspond to sample identifiers in the metadata.
/// * Contigs correspond to contig identifiers in the metadata.
/// * Each (sample, phase) combination corresponds to a haplotype in the metadata.
/// * Fragment field can be used as a fragment index for fragments from the sequence identified by (sample, contig, phase).
///   It can also be used as the starting offset of the fragment in the corresponding sequence.
///
/// See [`FullPathName`] for a stand-alone structure that contains the same information and documents some naming conventions.
#[repr(C)]
#[derive(Copy, Clone, Default, Debug, Hash, PartialEq, Eq)]
pub struct PathName {
    /// Sample identifier.
    pub sample: u32,

    /// Contig identifier.
    pub contig: u32,

    /// Phase / haplotype identifier.
    pub phase: u32,

    /// Fragment identifier / running count.
    pub fragment: u32,
}

impl PathName {
    /// Returns a new path name with all components set to 0.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns a path name with the given values in each field.
    pub fn from_fields(sample: usize, contig: usize, phase: usize, fragment: usize) -> Self {
        PathName {
            sample: sample as u32,
            contig: contig as u32,
            phase: phase as u32,
            fragment: fragment as u32,
        }
    }

    /// Returns the sample identifier.
    pub fn sample(&self) -> usize {
        self.sample as usize
    }

    /// Returns the contig identifier.
    pub fn contig(&self) -> usize {
        self.contig as usize
    }

    /// Returns the phase / haplotype identifier.
    pub fn phase(&self) -> usize {
        self.phase as usize
    }

    /// Returns the fragment identifier / running count.
    pub fn fragment(&self) -> usize {
        self.fragment as usize
    }

    /// Returns `true` if this is the same path or a predecessor fragment of the other path.
    pub fn is_predecessor_of(&self, other: &PathName) -> bool {
        self.sample == other.sample && self.contig == other.contig && self.phase == other.phase && self.fragment <= other.fragment
    }
}

impl Serializable for PathName {}

//-----------------------------------------------------------------------------

/// A structured path name as a stand-alone structure.
///
/// This structure stores the same information as [`PathName`]:
///
/// * Sample name.
/// * Contig name.
/// * Haplotype/phase number.
/// * Fragment number or starting offset of the fragment.
///
/// The path can be a generic path, a reference path, or a haplotype path in a similar way to vg path senses.
///
/// * Generic paths have [`GENERIC_SAMPLE`] as their sample name, and their actual name is stored as contig name.
///   Both haplotype and fragment numbers are `0`.
///   When serialized, haplotype number becomes [`GENERIC_HAPLOTYPE`] to follow vg conventions.
/// * Reference paths have `0` both as the haplotype and the fragment, and their names are of the form `sample#contig`.
/// * Haplotype paths start their haplotype numbers from `1`, and their names are of the form `sample#haplotype#contig@fragment`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FullPathName {
    /// Sample name.
    pub sample: String,

    /// Contig name.
    pub contig: String,

    /// Haplotype/phase number.
    pub haplotype: usize,

    /// Fragment number or starting offset of the fragment.
    pub fragment: usize,
}

impl FullPathName {
    /// Returns the name of the path with the given identifier in GBWT metadata, or [`None`] if the path does not exist.
    pub fn from_metadata(metadata: &Metadata, path_id: usize) -> Option<Self> {
        let name = metadata.path(path_id)?;
        let sample = metadata.sample_name(name.sample());
        let contig = metadata.contig_name(name.contig());
        Some(FullPathName {
            sample, contig,
            haplotype: name.phase(),
            fragment: name.fragment(),
        })
    }

    /// Returns a new generic path name.
    pub fn generic(name: &str) -> Self {
        FullPathName {
            sample: String::from(GENERIC_SAMPLE),
            contig: String::from(name),
            haplotype: 0,
            fragment: 0,
        }
    }

    /// Returns a new reference path name.
    pub fn reference(sample: &str, contig: &str) -> Self {
        FullPathName {
            sample: String::from(sample),
            contig: String::from(contig),
            haplotype: 0,
            fragment: 0,
        }
    }

    /// Returns a new haplotype path name.
    pub fn haplotype(sample: &str, contig: &str, haplotype: usize, fragment: usize) -> Self {
        FullPathName {
            sample: String::from(sample),
            contig: String::from(contig),
            haplotype,
            fragment,
        }
    }

    /// Returns a PanSN representation of the path name.
    ///
    /// The fragment field is assumed to be `0`.
    pub fn pan_sn_name(&self) -> String {
        format!("{}#{}#{}", self.sample, self.haplotype, self.contig)
    }

    /// Returns a path fragment name compatible with vg.
    ///
    /// The name is of the form `sample#haplotype#contig[fragment-end]`.
    pub fn path_fragment_name(&self, end: usize) -> String {
        format!("{}#{}#{}[{}-{}]", self.sample, self.haplotype, self.contig, self.fragment, end)
    }
}

impl fmt::Display for FullPathName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.sample == GENERIC_SAMPLE {
            write!(f, "{}", self.contig)
        } else if self.haplotype == 0 && self.fragment == 0 {
            write!(f, "{}#{}", self.sample, self.contig)
        } else {
            write!(f, "{}#{}#{}@{}", self.sample, self.haplotype, self.contig, self.fragment)
        }
    }
}

//-----------------------------------------------------------------------------

/// A builder that creates [`Metadata`] from [`FullPathName`] objects.
///
/// The created metadata will contain path, sample, and contig names.
/// Path names must be added in the same order as the paths in the corresponding [`crate::GBWT`] index.
/// See [`Metadata`] for an example.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MetadataBuilder {
    sample_names: Vec<String>,
    sample_set: HashMap<String, usize>,
    contig_names: Vec<String>,
    contig_set: HashMap<String, usize>,
    path_names: Vec<PathName>,
    path_set: HashSet<PathName>,
    haplotypes: HashSet<(usize, usize)>,
}

impl MetadataBuilder {
    /// Creates a new metadata builder.
    pub fn new() -> Self {
        MetadataBuilder::default()
    }

    // Returns the identifier of the given sample name.
    fn sample_id(&self, name: &FullPathName) -> usize {
        *self.sample_set.get(&name.sample).unwrap_or(&self.sample_set.len())
    }

    // Returns the identifier of the given contig name.
    fn contig_id(&self, name: &FullPathName) -> usize {
        *self.contig_set.get(&name.contig).unwrap_or(&self.contig_set.len())
    }

    // Returns the identifier of the given name in the union of two name sets.
    // Identifiers in `first` come before identifiers in `second`.
    // If the name is not present, returns the new identifier that would be assigned to the name if it were added to `second`.
    fn joint_id(first: &HashMap<String, usize>, second: &HashMap<String, usize>, name: &String) -> usize {
        if let Some(id) = first.get(name) {
            *id
        } else {
            *second.get(name).unwrap_or(&(first.len() + second.len()))
        }
    }

    /// Inserts a path name to the builder.
    ///
    /// # Errors
    ///
    /// Returns an error if the path name is a duplicate of an already inserted path name.
    /// If an error occurs, the builder is not modified.
    pub fn insert(&mut self, name: &FullPathName) -> Result<(), String> {
        let sample_id = self.sample_id(name);
        let contig_id = self.contig_id(name);
        let path_name = PathName::from_fields(sample_id, contig_id, name.haplotype, name.fragment);
        if self.path_set.contains(&path_name) {
            return Err(format!("MetadataBuilder: Duplicate path name {}", name));
        }

        if sample_id == self.sample_set.len() {
            self.sample_names.push(name.sample.clone());
            self.sample_set.insert(name.sample.clone(), sample_id);
        }
        if contig_id == self.contig_set.len() {
            self.contig_names.push(name.contig.clone());
            self.contig_set.insert(name.contig.clone(), contig_id);
        }
        self.path_names.push(path_name);
        self.path_set.insert(path_name);
        self.haplotypes.insert((sample_id, name.haplotype));

        Ok(())
    }

    /// Appends a slice of path names to the builder.
    ///
    /// # Errors
    ///
    /// Returns an error if there are duplicate path names in the input or between the input and the existing path names.
    /// If an error occurs, the builder is not modified.
    pub fn extend(&mut self, names: &[FullPathName]) -> Result<(), String> {
        // Build the extension in a separate builder.
        let mut extension = MetadataBuilder::new();
        for name in names {
            let sample_id = Self::joint_id(&self.sample_set, &extension.sample_set, &name.sample);
            if sample_id >= self.sample_set.len() + extension.sample_set.len() {
                extension.sample_names.push(name.sample.clone());
                extension.sample_set.insert(name.sample.clone(), sample_id);
            }
            let contig_id = Self::joint_id(&self.contig_set, &extension.contig_set, &name.contig);
            if contig_id >= self.contig_set.len() + extension.contig_set.len() {
                extension.contig_names.push(name.contig.clone());
                extension.contig_set.insert(name.contig.clone(), contig_id);
            }
            let path_name = PathName::from_fields(sample_id, contig_id, name.haplotype, name.fragment);
            if self.path_set.contains(&path_name) || extension.path_set.contains(&path_name) {
                return Err(format!("MetadataBuilder: Duplicate path name {}", name));
            }
            extension.path_names.push(path_name);
            extension.path_set.insert(path_name);
            extension.haplotypes.insert((sample_id, name.haplotype));
        }

        // Append the extension to the builder.
        self.sample_names.extend(extension.sample_names);
        self.sample_set.extend(extension.sample_set);
        self.contig_names.extend(extension.contig_names);
        self.contig_set.extend(extension.contig_set);
        self.path_names.extend(extension.path_names);
        self.path_set.extend(extension.path_set);
        self.haplotypes.extend(extension.haplotypes);

        Ok(())
    }

    /// Returns the name of the path with the given identifier, or [`None`] if there is no such name.
    pub fn path_name(&self, id: usize) -> Option<FullPathName> {
        if id < self.path_names.len() {
            let path_name = self.path_names[id];
            let sample = self.sample_names.get(path_name.sample())?.clone();
            let contig = self.contig_names.get(path_name.contig())?.clone();
            Some(FullPathName {
                sample, contig,
                haplotype: path_name.phase(),
                fragment: path_name.fragment(),
            })
        } else {
            None
        }
    }
}

//-----------------------------------------------------------------------------
