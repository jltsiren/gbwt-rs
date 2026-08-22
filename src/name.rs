//! Structures for storing graph name and relationship information.
//!
//! A [`GraphName`] can store a graph name along with subgraph and translation relationships between graphs.
//! This information can be imported from and exported to [`Tags`] objects as well as GFA and GAF header lines.
//!
//! These structures were originally designed for use with stable [pggname](https://github.com/jltsiren/pggname) graph names.
//! They can also be used with other naming schemes.
//! The GFA representation in optional fields of type `Z` restricts the names to printable ASCII, including space.
//! Since commas (`,`) and semicolons (`;`) are used as separators, they are not allowed in graph names.
//!
//! Graph A is a subgraph of graph B, if all nodes and edges of A are also present in B.
//! The node identifiers and sequence labels in graph A must match those in graph B.
//! This usually means that anything using graph A as a reference can also use graph B.
//!
//! For coordinate translation from graph A to graph B, we use an intermediate graph C that is a subgraph of B.
//! We require that graphs A and C are isomorphic (with matching node labels), if we break their nodes into 1 bp pieces.
//! There is therefore a one-to-one mapping between unary paths in A and C.
//! We can use this mapping to translate positions in graph A to graph C, and then use these positions in graph B.

use crate::{GBZ, Tags};

use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// A structure that stores a graph name along with subgraph and translation relationships between graphs.
///
/// Each object corresponds to a particular graph.
/// The object may store the name of the graph, if available.
/// It may contain the relationship between this graph and its parent graph, if the relationship and the parent's name are known.
/// There may also be other relationships inherited from the parent graph.
///
/// # Examples
///
/// ```
/// use gbz::GraphName;
///
/// let parent = GraphName::new(String::from("parent"));
/// let mut translated = GraphName::new(String::from("translated"));
/// translated.add_translation_to(&parent);
/// let mut subgraph = GraphName::new(String::from("subgraph"));
/// subgraph.make_subgraph_of(&translated);
///
/// assert!(subgraph.is_subgraph_of(&translated));
/// assert!(!subgraph.is_subgraph_of(&parent));
/// assert!(translated.translates_to(&parent));
/// assert!(subgraph.translates_to(&parent));
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct GraphName {
    name: Option<String>,
    subgraph: BTreeMap<String, BTreeSet<String>>,
    translation: BTreeMap<String, BTreeSet<String>>,
}

/// Constants.
impl GraphName {
    /// Name of the [`Tags`] key storing the graph name.
    const TAG_NAME: &'static str = "pggname";

    /// Name of the [`Tags`] key storing subgraph relationships.
    const TAG_SUBGRAPH: &'static str = "subgraph";

    /// Name of the [`Tags`] key storing translation relationships.
    const TAG_TRANSLATION: &'static str = "translation";

    /// GFA header tag storing the graph name.
    const GFA_HEADER_NAME: &'static str = "NM";

    /// GAF header tag storing the graph name.
    const GAF_HEADER_NAME: &'static str = "RN";

    /// GFA/GAF header tag storing subgraph relationships.
    const GFA_GAF_HEADER_SUBGRAPH: &'static str = "SG";

    /// GFA/GAF header tag storing translation relationships.
    const GFA_GAF_HEADER_TRANSLATION: &'static str = "TL";

    const GFA_HEADER_TYPE: &'static str = "H";
    const GAF_HEADER_PREFIX: &'static str = "@"; 
    const GFA_GAF_FIELD_SEPARATOR: char = '\t';
    const TAG_GFA_RELATIONSHIP_SEPARATOR: char = ',';
    const TAG_RELATIONSHIP_LIST_SEPARATOR: char = ';';
}

//-----------------------------------------------------------------------------

/// Construction.
impl GraphName {
    /// Creates a new `GraphName` with the given graph name.
    pub fn new(name: String) -> Self {
        GraphName {
            name: Some(name),
            subgraph: BTreeMap::new(),
            translation: BTreeMap::new(),
        }
    }

    /// Parses a `GraphName` from the given tags.
    ///
    /// Returns an error if tag values are malformed.
    pub fn from_tags(tags: &Tags) -> Result<Self, String> {
        let mut result = GraphName::default();

        if let Some(name_field) = tags.get(Self::TAG_NAME) {
            result.name = Some(String::from(name_field));
        }

        if let Some(subgraph_field) = tags.get(Self::TAG_SUBGRAPH) {
            let relationships: Vec<&str> = subgraph_field.split(Self::TAG_RELATIONSHIP_LIST_SEPARATOR).collect();
            for rel in relationships {
                let parts: Vec<&str> = rel.split(Self::TAG_GFA_RELATIONSHIP_SEPARATOR).collect();
                if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
                    return Err(format!("Invalid subgraph relationship: {}", rel));
                }
                result.subgraph
                    .entry(String::from(parts[0]))
                    .or_default()
                    .insert(String::from(parts[1]));
            }
        }

        if let Some(translation_field) = tags.get(Self::TAG_TRANSLATION) {
            let relationships: Vec<&str> = translation_field.split(Self::TAG_RELATIONSHIP_LIST_SEPARATOR).collect();
            for rel in relationships {
                let parts: Vec<&str> = rel.split(Self::TAG_GFA_RELATIONSHIP_SEPARATOR).collect();
                if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
                    return Err(format!("Invalid translation relationship: {}", rel));
                }
                result.translation
                    .entry(String::from(parts[0]))
                    .or_default()
                    .insert(String::from(parts[1]));
            }
        }

        Ok(result)
    }

    /// Parses a `GraphName` from the tags in the given GBZ graph.
    ///
    /// Returns an empty object if the tags cannot be parsed.
    pub fn from_gbz(gbz: &GBZ) -> Self {
        Self::from_tags(gbz.tags()).unwrap_or_default()
    }

    fn typed_field_is_string(field: &str) -> Result<bool, String> {
        let bytes = field.as_bytes();
        if field.len() < 5 || bytes[2] != b':' || bytes[4] != b':' {
            return Err(format!("Invalid GFA typed field: {}", field));
        }
        Ok(bytes[3] == b'Z')
    }

    fn parse_gfa_optional_fields(fields: &[&str], result: &mut GraphName) -> Result<(), String> {
        for &field in fields {
            if !Self::typed_field_is_string(field)? {
                continue;
            }
            let tag = &field[0..2];
            let value = &field[5..];
            match tag {
                Self::GFA_HEADER_NAME => {
                    result.name = Some(String::from(value));
                }
                Self::GFA_GAF_HEADER_SUBGRAPH => {
                    let parts: Vec<&str> = value.split(Self::TAG_GFA_RELATIONSHIP_SEPARATOR).collect();
                    if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
                        return Err(format!("Invalid subgraph field: {}", field));
                    }
                    result.subgraph
                        .entry(String::from(parts[0]))
                        .or_default()
                        .insert(String::from(parts[1]));
                }
                Self::GFA_GAF_HEADER_TRANSLATION => {
                    let parts: Vec<&str> = value.split(Self::TAG_GFA_RELATIONSHIP_SEPARATOR).collect();
                    if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
                        return Err(format!("Invalid translation field: {}", field));
                    }
                    result.translation
                        .entry(String::from(parts[0]))
                        .or_default()
                        .insert(String::from(parts[1]));
                }
                _ => {}
            }
        }
        Ok(())
    }

    fn parse_gaf_header_fields(line: &str, fields: &[&str], result: &mut GraphName) -> Result<(), String> {
        match &fields[0][1..] {
            Self::GAF_HEADER_NAME => {
                if fields.len() != 2 || fields[1].is_empty() {
                    return Err(format!("Invalid GAF name header line: {}", line));
                }
                result.name = Some(String::from(fields[1]));
            }
            Self::GFA_GAF_HEADER_SUBGRAPH => {
                if fields.len() != 3 || fields[1].is_empty() || fields[2].is_empty() {
                    return Err(format!("Invalid GAF subgraph header line: {}", line));
                }
                result.subgraph
                    .entry(String::from(fields[1]))
                    .or_default()
                    .insert(String::from(fields[2]));
            }
            Self::GFA_GAF_HEADER_TRANSLATION => {
                if fields.len() != 3 || fields[1].is_empty() || fields[2].is_empty() {
                    return Err(format!("Invalid GAF translation header line: {}", line));
                }
                result.translation
                    .entry(String::from(fields[1]))
                    .or_default()
                    .insert(String::from(fields[2]));
            }
            _ => {}
        }
        Ok(())
    }

    /// Parses a `GraphName` from the given GFA/GAF header lines.
    ///
    /// The lines must not end with a newline.
    /// Returns an error if the lines cannot be parsed.
    pub fn from_header_lines(lines: &[String]) -> Result<Self, String> {
        let mut result = GraphName::default();

        for (i, line) in lines.iter().enumerate() {
            let fields: Vec<&str> = line.split(Self::GFA_GAF_FIELD_SEPARATOR).collect();
            if fields.len() < 2 {
                return Err(format!("Error parsing header line {}: not enough fields", i + 1));
            }
            if fields[0] == Self::GFA_HEADER_TYPE {
                Self::parse_gfa_optional_fields(&fields[1..], &mut result)?;
            } else if fields[0].len() == 3 && fields[0].starts_with(Self::GAF_HEADER_PREFIX) {
                Self::parse_gaf_header_fields(line, &fields, &mut result)?;
            } else {
                return Err(format!("Error parsing header line {}: unknown first field {}", i + 1, fields[0]));
            }
        }

        Ok(result)
    }

    /// Adds a new subgraph relationship, if both names are non-empty.
    pub fn add_subgraph(&mut self, subgraph: &str, supergraph: &str) {
        if !subgraph.is_empty() && !supergraph.is_empty() {
            self.subgraph
                .entry(String::from(subgraph))
                .or_default()
                .insert(String::from(supergraph));
        }
    }

    /// Adds a new subgraph relationship between this graph and the parent graph.
    ///
    /// Also copies all relationships from the parent graph.
    /// Does nothing if either name is not available.
    pub fn make_subgraph_of(&mut self, parent: &GraphName) {
        if !self.has_name() || !parent.has_name() {
            return;
        }
        let entry = self.subgraph.entry(self.name.as_ref().unwrap().clone()).or_default();
        entry.insert(parent.name.as_ref().unwrap().clone());
        self.add_relationships(parent);
    }

    /// Adds a new translation relationship, if both names are non-empty.
    pub fn add_translation(&mut self, from: &str, to: &str) {
        if !from.is_empty() && !to.is_empty() {
            self.translation
                .entry(String::from(from))
                .or_default()
                .insert(String::from(to));
        }
    }

    /// Adds a new translation relationship between this graph and the parent graph.
    ///
    /// Also copies all relationships from the parent graph.
    /// Does nothing if either name is not available.
    pub fn add_translation_to(&mut self, parent: &GraphName) {
        if !self.has_name() || !parent.has_name() {
            return;
        }
        let entry = self.translation.entry(self.name.as_ref().unwrap().clone()).or_default();
        entry.insert(parent.name.as_ref().unwrap().clone());
        self.add_relationships(parent);
    }

    /// Adds all relationships from another `GraphName` object.
    pub fn add_relationships(&mut self, other: &GraphName) {
        for (supergraph, subgraphs) in &other.subgraph {
            let entry = self.subgraph.entry(supergraph.clone()).or_default();
            for subgraph in subgraphs {
                entry.insert(subgraph.clone());
            }
        }
        for (from, tos) in &other.translation {
            let entry = self.translation.entry(from.clone()).or_default();
            for to in tos {
                entry.insert(to.clone());
            }
        }
    }
}

//-----------------------------------------------------------------------------

/// Export to other formats.
impl GraphName {
    fn relationships_to_string(relationships: &BTreeMap<String, BTreeSet<String>>) -> String {
        let mut value = String::new();
        for (from, tos) in relationships {
            for to in tos {
                if !value.is_empty() {
                    value.push(Self::TAG_RELATIONSHIP_LIST_SEPARATOR);
                }
                value.push_str(from);
                value.push(Self::TAG_GFA_RELATIONSHIP_SEPARATOR);
                value.push_str(to);
            }
        }
        value
    }

    /// Writes the data stored in this object to the given tags.
    ///
    /// Clears existing tags if no corresponding data is available.
    pub fn set_tags(&self, tags: &mut Tags) {
        if let Some(name) = &self.name {
            tags.insert(Self::TAG_NAME, name);
        } else {
            tags.remove(Self::TAG_NAME);
        }

        if !self.subgraph.is_empty() {
            let value = Self::relationships_to_string(&self.subgraph);
            tags.insert(Self::TAG_SUBGRAPH, &value);
        } else {
            tags.remove(Self::TAG_SUBGRAPH);
        }

        if !self.translation.is_empty() {
            let value = Self::relationships_to_string(&self.translation);
            tags.insert(Self::TAG_TRANSLATION, &value);
        } else {
            tags.remove(Self::TAG_TRANSLATION);
        }
    }

    /// Returns GFA header lines representing this object.
    ///
    /// The lines do not end with a newline.
    pub fn to_gfa_header_lines(&self) -> Vec<String> {
        let mut lines = Vec::new();
        if let Some(name) = &self.name {
            lines.push(format!("H\t{}:Z:{}", Self::GFA_HEADER_NAME, name));
        }
        for (subgraph, supergraphs) in &self.subgraph {
            for supergraph in supergraphs {
                lines.push(format!("H\t{}:Z:{},{}", Self::GFA_GAF_HEADER_SUBGRAPH, subgraph, supergraph));
            }
        }
        for (from, tos) in &self.translation {
            for to in tos {
                lines.push(format!("H\t{}:Z:{},{}", Self::GFA_GAF_HEADER_TRANSLATION, from, to));
            }
        }
        lines
    }

    /// Returns GAF header lines representing this object.
    ///
    /// The lines do not end with a newline.
    pub fn to_gaf_header_lines(&self) -> Vec<String> {
        let mut lines = Vec::new();
        if let Some(name) = &self.name {
            lines.push(format!("@{}\t{}", Self::GAF_HEADER_NAME, name));
        }
        for (subgraph, supergraphs) in &self.subgraph {
            for supergraph in supergraphs {
                lines.push(format!("@{}\t{}\t{}", Self::GFA_GAF_HEADER_SUBGRAPH, subgraph, supergraph));
            }
        }
        for (from, tos) in &self.translation {
            for to in tos {
                lines.push(format!("@{}\t{}\t{}", Self::GFA_GAF_HEADER_TRANSLATION, from, to));
            }
        }
        lines
    }
}

/// Queries and operations.
impl GraphName {
    /// Returns the name of the graph, if available.
    pub fn name(&self) -> Option<&String> {
        self.name.as_ref()
    }

    /// Returns `true` if the graph has a name.
    pub fn has_name(&self) -> bool {
        self.name.is_some()
    }

    /// Returns `true` if both objects represent the same graph.
    pub fn is_same(&self, other: &GraphName) -> bool {
        match (&self.name, &other.name) {
            (Some(name1), Some(name2)) => name1 == name2,
            _ => false,
        }
    }

    /// Returns an iterator over stored subgraph relationships.
    ///
    /// The iterator yields pairs `(subgraph_name, supergraph_name)` in sorted order.
    pub fn subgraph_iter(&self) -> impl Iterator<Item = (&str, &str)> {
        self.subgraph.iter().flat_map(|(subgraph, supergraphs)| {
            supergraphs.iter().map(move |supergraph| (subgraph.as_ref(), supergraph.as_ref()))
        })
    }

    /// Returns an iterator over stored translation relationships.
    ///
    /// The iterator yields pairs `(from_name, to_name)` in sorted order.
    pub fn translation_iter(&self) -> impl Iterator<Item = (&str, &str)> {
        self.translation.iter().flat_map(|(from, tos)| {
            tos.iter().map(move |to| (from.as_ref(), to.as_ref()))
        })
    }

    // Finds a path of subgraph relationships from `from` to `to`, including both.
    // Uses relationships stored in `self`.
    fn find_subgraph_path(&self, from: &GraphName, to: &GraphName) -> Option<Vec<String>> {
        if !from.has_name() || !to.has_name() {
            return None;
        }
        let from_name = from.name().unwrap();
        let to_name = to.name().unwrap();

        // Find a shortest path using BFS.
        let mut predecessor: BTreeMap<String, String> = BTreeMap::new();
        predecessor.insert(from_name.clone(), String::new());
        let mut queue: VecDeque<String> = VecDeque::new();
        queue.push_back(from_name.clone());
        while let Some(curr) = queue.pop_front() {
            if curr == *to_name {
                break;
            }
            if let Some(supers) = self.subgraph.get(&curr) {
                for supergraph in supers {
                    if !predecessor.contains_key(supergraph) {
                        predecessor.insert(supergraph.clone(), curr.clone());
                        queue.push_back(supergraph.clone());
                    }
                }
            }
        }
        if !predecessor.contains_key(to_name) {
            return None;
        }

        // Trace back the path.
        let mut result: Vec<String> = Vec::new();
        let mut current = to_name.clone();
        while !current.is_empty() {
            result.push(current.clone());
            current = predecessor.get(&current).unwrap().clone();
        }
        result.reverse();

        Some(result)
    }

    // Finds a path of subgraph or translation relationships from `from` to `to`, including both.
    // Each step is a pair `(name, is_translation)`, where `is_translation` indicates whether the step to the next name is a translation.
    // Uses relationships stored in `self`.
    fn find_path(&self, from: &GraphName, to: &GraphName) -> Option<Vec<(String, bool)>> {
        if !from.has_name() || !to.has_name() {
            return None;
        }
        let from_name = from.name().unwrap();
        let to_name = to.name().unwrap();

        // Find a shortest path using BFS.
        let mut predecessor: BTreeMap<String, (String, bool)> = BTreeMap::new();
        predecessor.insert(from_name.clone(), (String::new(), false));
        let mut queue: VecDeque<String> = VecDeque::new();
        queue.push_back(from_name.clone());
        while let Some(curr) = queue.pop_front() {
            if curr == *to_name {
                break;
            }
            // Prioritize subgraph relationships.
            if let Some(neighbors) = self.subgraph.get(&curr) {
                for next in neighbors {
                    if !predecessor.contains_key(next) {
                        predecessor.insert(next.clone(), (curr.clone(), false));
                        queue.push_back(next.clone());
                    }
                }
            }
            // Then consider translation relationships.
            if let Some(neighbors) = self.translation.get(&curr) {
                for next in neighbors {
                    if !predecessor.contains_key(next) {
                        predecessor.insert(next.clone(), (curr.clone(), true));
                        queue.push_back(next.clone());
                    }
                }
            }
        }
        if !predecessor.contains_key(to_name) {
            return None;
        }

        // Trace back the path.
        let mut result: Vec<(String, bool)> = Vec::new();
        result.push((from_name.clone(), false));
        let (mut curr, mut is_translation) = predecessor.get(to_name).unwrap().clone();
        while !curr.is_empty() {
            result.push((curr.clone(), is_translation));
            (curr, is_translation) = predecessor.get(&curr).unwrap().clone();
        }
        result.reverse();

        Some(result)
    }

    /// Returns `true` if this graph is a subgraph of the given graph.
    ///
    /// Uses relationships stored in both graphs.
    pub fn is_subgraph_of(&self, other: &GraphName) -> bool {
        let mut merged = self.clone();
        merged.add_relationships(other);
        merged.find_subgraph_path(self, other).is_some()
    }

    /// Returns `true` if coordinates in this graph can be translated to coordinates in the given graph.
    ///
    /// Uses relationships stored in both graphs.
    pub fn translates_to(&self, other: &GraphName) -> bool {
        let mut merged = self.clone();
        merged.add_relationships(other);
        merged.find_path(self, other).is_some()
    }

    fn append_description(result: &mut String, num: usize, description: &str) {
        let line = format!("Name {} is for {}\n", num, description);
        result.push_str(&line); 
    }

    fn append_relationship(result: &mut String, step: usize, is_translation: bool) {
        let relation = if is_translation {
            "translates to"
        } else {
            "is a subgraph of"
        };
        let line = format!("Graph {} {} graph {}\n", step, relation, step + 1);
        result.push_str(&line);
    }

    fn append_graph(result: &mut String, num: usize, name: &str) {
        let line = format!("{}\t{}\n", num, name);
        result.push_str(&line);
    }

    /// Returns a description of the relationship between this graph and the given graph.
    ///
    /// Uses relationships stored in both graphs.
    /// The description consists of multiple lines and ends with a newline.
    ///
    /// # Arguments
    ///
    /// * `other`: Name of the other graph.
    /// * `self_desc`: Description of this graph to use in the output.
    /// * `other_desc`: Description of the other graph to use in the output.
    pub fn describe_relationship(&self, other: &GraphName, self_desc: &str, other_desc: &str) -> String {
        let mut merged = self.clone();
        merged.add_relationships(other);

        let no_name = String::from("(no name)");
        let mut from = (self.name.as_ref().unwrap_or(&no_name).clone(), String::from(self_desc));
        let mut to = (other.name.as_ref().unwrap_or(&no_name).clone(), String::from(other_desc));
        let mut path = merged.find_path(self, other);
        if path.is_none() {
            std::mem::swap(&mut from, &mut to);
            path = merged.find_path(other, self);
        }

        // Graph descriptions and relationships.
        let mut result = String::new();
        Self::append_description(&mut result, 1, &from.1);
        if let Some(path) = &path {
            for i in 1..path.len() {
                Self::append_relationship(&mut result, i, path[i - 1].1);
            }
            Self::append_description(&mut result, path.len(), &to.1);
        } else {
            Self::append_description(&mut result, 2, &to.1);
        }

        // Graph names.
        result.push_str("With graph names:\n");
        if let Some(path) = path {
            for (i, (name, _)) in path.iter().enumerate() {
                Self::append_graph(&mut result, i + 1, name);
            }
        } else {
            Self::append_graph(&mut result, 1, &from.0);
            Self::append_graph(&mut result, 2, &to.0);
        }

        result
    }
}

//-----------------------------------------------------------------------------
