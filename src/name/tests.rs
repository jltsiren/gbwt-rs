use super::*;

//-----------------------------------------------------------------------------

const NAME: &str = "A";
const SUBGRAPH: &[(&str, &str)] = &[
    ("A", "B"),
    ("C", "D"),
    ("D", "E"),
];
const TRANSLATION: &[(&str, &str)] = &[
    ("B", "C"),
    ("C", "F"),
];

fn manual() -> GraphName {
    let mut name = GraphName::new(String::from(NAME));
    for (subgraph, supergraph) in SUBGRAPH.iter() {
        name.add_subgraph(subgraph, supergraph);
    }
    for (from, to) in TRANSLATION.iter() {
        name.add_translation(from, to);
    }
    name
}

fn from_parents() -> GraphName {
    let f = GraphName::new(String::from("F"));
    let e = GraphName::new(String::from("E"));
    let mut d = GraphName::new(String::from("D"));
    d.make_subgraph_of(&e);
    let mut c = GraphName::new(String::from("C"));
    c.make_subgraph_of(&d);
    c.add_translation_to(&f);
    let mut b = GraphName::new(String::from("B"));
    b.add_translation_to(&c);
    let mut a = GraphName::new(String::from(NAME));
    a.make_subgraph_of(&b);
    a
}

// In the following functions, the first returned value also contains fields (e.g. tags, rows)
// unrelated to GraphName, while the second only contains those derived from GraphName.

fn tags() -> (Tags, Tags) {
    let mut all_tags = Tags::new();
    let mut name_tags = Tags::new();
    all_tags.insert(crate::SOURCE_KEY, crate::SOURCE_VALUE);
    all_tags.insert(GraphName::TAG_NAME, NAME);
    name_tags.insert(GraphName::TAG_NAME, NAME);

    let mut subgraph_value = String::new();
    for (subgraph, supergraph) in SUBGRAPH.iter() {
        if !subgraph_value.is_empty() {
            subgraph_value.push(GraphName::TAG_RELATIONSHIP_LIST_SEPARATOR);
        }
        subgraph_value.push_str(subgraph);
        subgraph_value.push(GraphName::TAG_GFA_RELATIONSHIP_SEPARATOR);
        subgraph_value.push_str(supergraph);
    }
    all_tags.insert(GraphName::TAG_SUBGRAPH, &subgraph_value);
    name_tags.insert(GraphName::TAG_SUBGRAPH, &subgraph_value);

    let mut translation_value = String::new();
    for (from, to) in TRANSLATION.iter() {
        if !translation_value.is_empty() {
            translation_value.push(GraphName::TAG_RELATIONSHIP_LIST_SEPARATOR);
        }
        translation_value.push_str(from);
        translation_value.push(GraphName::TAG_GFA_RELATIONSHIP_SEPARATOR);
        translation_value.push_str(to);
    }
    all_tags.insert(GraphName::TAG_TRANSLATION, &translation_value);
    name_tags.insert(GraphName::TAG_TRANSLATION, &translation_value);

    (all_tags, name_tags)
}

fn gfa_header_lines() -> (Vec<String>, Vec<String>) {
    let mut all_lines = Vec::new();
    all_lines.push(String::from("H\tVN:Z:1.1"));
    all_lines.push(format!("H\t{}:Z:{}", GraphName::GFA_HEADER_NAME, NAME));

    for (subgraph, supergraph) in SUBGRAPH.iter() {
        all_lines.push(format!(
            "H\t{}:Z:{}{}{}",
            GraphName::GFA_GAF_HEADER_SUBGRAPH,
            subgraph,
            GraphName::TAG_GFA_RELATIONSHIP_SEPARATOR,
            supergraph
        ));
    }

    for (from, to) in TRANSLATION.iter() {
        all_lines.push(format!(
            "H\t{}:Z:{}{}{}",
            GraphName::GFA_GAF_HEADER_TRANSLATION,
            from,
            GraphName::TAG_GFA_RELATIONSHIP_SEPARATOR,
            to
        ));
    }

    let name_lines: Vec<String> = all_lines.iter().skip(1).cloned().collect();
    (all_lines, name_lines)
}

fn gaf_header_lines() -> (Vec<String>, Vec<String>) {
    let mut all_lines = Vec::new();
    all_lines.push(String::from("@HD\tVN:Z:1.0"));
    all_lines.push(format!("@{}\t{}", GraphName::GAF_HEADER_NAME, NAME));

    for (subgraph, supergraph) in SUBGRAPH.iter() {
        all_lines.push(format!(
            "@{}\t{}\t{}",
            GraphName::GFA_GAF_HEADER_SUBGRAPH,
            subgraph,
            supergraph
        ));
    }

    for (from, to) in TRANSLATION.iter() {
        all_lines.push(format!(
            "@{}\t{}\t{}",
            GraphName::GFA_GAF_HEADER_TRANSLATION,
            from,
            to
        ));
    }

    let name_lines: Vec<String> = all_lines.iter().skip(1).cloned().collect();
    (all_lines, name_lines)
}

fn expected_description_lines(steps: usize, has_path: bool) -> usize {
    let mut lines = 2; // Descriptions of both names.
    lines += steps; // One line per step.
    lines += 1; // With graph names.
    lines += steps + 1; // One line per graph.
    if !has_path {
        lines += 1; // Final graph.
    }
    lines
}

fn test_describe_relationship(from: &GraphName, to: &GraphName, from_desc: &str, to_desc: &str, steps: usize, has_path: bool) {
    let description = from.describe_relationship(to, from_desc, to_desc);
    let expected_lines = expected_description_lines(steps, has_path);
    let description_lines: Vec<&str> = description.lines().collect();
    let from_name = from.name().cloned().unwrap_or(String::from("<no name>"));
    let to_name = to.name().cloned().unwrap_or(String::from("<no name>"));
    assert_eq!(description_lines.len(), expected_lines, "Unexpected number of description lines for {} and {}", from_name, to_name);
}

#[test]
fn graph_name_empty() {
    let default = GraphName::default();
    assert!(!default.has_name(), "Expected has_name() to be false in default GraphName");
    assert!(default.name().is_none(), "Expected no name in default GraphName");

    let empty_tags = Tags::new();
    let from_tags = GraphName::from_tags(&empty_tags);
    assert!(from_tags.is_ok(), "Failed to build GraphName from empty tags: {}", from_tags.unwrap_err());
    let from_tags = from_tags.unwrap();
    assert!(!from_tags.has_name(), "Expected has_name() to be false in GraphName built from empty tags");
    assert!(from_tags.name().is_none(), "Expected no name in GraphName built from empty tags");
    assert!(!from_tags.is_same(&default), "GraphNames with missing names should not be the same");
    let mut to_tags = Tags::new();
    default.set_tags(&mut to_tags);
    assert_eq!(to_tags, empty_tags, "Expected no tags to be written from default GraphName");
    test_describe_relationship(&default, &from_tags, "default", "from_tags", 0, false);

    let empty_header_lines = Vec::new();
    let from_headers = GraphName::from_header_lines(&empty_header_lines);
    assert!(from_headers.is_ok(), "Failed to build GraphName from empty header lines: {}", from_headers.unwrap_err());
    let from_headers = from_headers.unwrap();
    assert!(!from_headers.has_name(), "Expected has_name() to be false in GraphName built from empty header lines");
    assert!(from_headers.name().is_none(), "Expected no name in GraphName built from empty header lines");
    assert!(!from_headers.is_same(&default), "GraphNames with missing names should not be the same");
    test_describe_relationship(&default, &from_headers, "default", "from_headers", 0, false);
    let gfa_headers = default.to_gfa_header_lines();
    assert!(gfa_headers.is_empty(), "Expected no GFA header lines from default GraphName");
    let gaf_headers = default.to_gaf_header_lines();
    assert!(gaf_headers.is_empty(), "Expected no GAF header lines from default GraphName");
}

#[test]
fn graph_name_manual() {
    let manual = manual();
    let from_parents = from_parents();
    assert!(manual.is_same(&from_parents), "GraphName built manually does not match GraphName built from parents");
    assert_eq!(manual, from_parents, "GraphName built manually is not equal to GraphName built from parents");

    let subgraph: Vec<(&str, &str)> = manual.subgraph_iter().collect();
    assert_eq!(subgraph.len(), SUBGRAPH.len(), "Wrong number of subgraph relationships from iterator");
    for (i, pair) in subgraph.iter().enumerate() {
        assert_eq!(*pair, SUBGRAPH[i], "Wrong subgraph relationship {} from iterator", i);
    }

    let translation: Vec<(&str, &str)> = manual.translation_iter().collect();
    assert_eq!(translation.len(), TRANSLATION.len(), "Wrong number of translation relationships from iterator");
    for (i, pair) in translation.iter().enumerate() {
        assert_eq!(*pair, TRANSLATION[i], "Wrong translation relationship {} from iterator", i);
    }
}

#[test]
fn graph_name_tags() {
    let (all_tags, name_tags) = tags();
    let from_tags = GraphName::from_tags(&all_tags);
    assert!(from_tags.is_ok(), "Failed to build GraphName from tags: {}", from_tags.unwrap_err());
    let from_tags = from_tags.unwrap();
    assert!(from_tags.has_name(), "GraphName built from tags should have a name");
    assert_eq!(from_tags.name().unwrap(), NAME, "Wrong name in GraphName built from tags");

    let from_manual = manual();
    assert!(from_tags.is_same(&from_manual), "GraphName built from tags does not match manual GraphName");
    assert_eq!(from_tags, from_manual, "GraphName built from tags is not equal to manual GraphName");

    let mut to_tags = Tags::new();
    from_tags.set_tags(&mut to_tags);
    assert_eq!(to_tags, name_tags, "Tags written from GraphName do not match expected tags");

    let empty_name = GraphName::default();
    empty_name.set_tags(&mut to_tags);
    assert!(!to_tags.contains_key(GraphName::TAG_NAME), "Graph name tag was not cleared");
    assert!(!to_tags.contains_key(GraphName::TAG_SUBGRAPH), "Subgraph tag was not cleared");
    assert!(!to_tags.contains_key(GraphName::TAG_TRANSLATION), "Translation tag was not cleared");
}

#[test]
fn graph_name_gfa() {
    let (all_headers, name_headers) = gfa_header_lines();
    let from_headers = GraphName::from_header_lines(&all_headers);
    assert!(from_headers.is_ok(), "Failed to build GraphName from GFA header lines: {}", from_headers.unwrap_err());
    let from_headers = from_headers.unwrap();
    assert!(from_headers.has_name(), "GraphName built from GFA header lines should have a name");
    assert_eq!(from_headers.name().unwrap(), NAME, "Wrong name in GraphName built from GFA header lines");

    let from_manual = manual();
    assert!(from_headers.is_same(&from_manual), "GraphName built from GFA header lines does not match manual GraphName");
    assert_eq!(from_headers, from_manual, "GraphName built from GFA header lines is not equal to manual GraphName");

    let to_headers = from_headers.to_gfa_header_lines();
    assert_eq!(to_headers, name_headers, "GFA header lines written from GraphName do not match expected header lines");
}

#[test]
fn graph_name_gaf() {
    let (all_headers, name_headers) = gaf_header_lines();
    let from_headers = GraphName::from_header_lines(&all_headers);
    assert!(from_headers.is_ok(), "Failed to build GraphName from GAF header lines: {}", from_headers.unwrap_err());
    let from_headers = from_headers.unwrap();
    assert!(from_headers.has_name(), "GraphName built from GAF header lines should have a name");
    assert_eq!(from_headers.name().unwrap(), NAME, "Wrong name in GraphName built from GAF header lines");

    let from_manual = manual();
    assert!(from_headers.is_same(&from_manual), "GraphName built from GAF header lines does not match manual GraphName");
    assert_eq!(from_headers, from_manual, "GraphName built from GAF header lines is not equal to manual GraphName");

    let to_headers = from_headers.to_gaf_header_lines();
    assert_eq!(to_headers, name_headers, "GAF header lines written from GraphName do not match expected header lines");
}

#[test]
fn graph_name_subgraph() {
    let a = manual();
    let a_empty = GraphName::new(String::from(NAME));

    // Same graph.
    assert!(a.is_subgraph_of(&a_empty), "Graph is not a subgrap of itself (relationships in subgraph)");
    assert!(a_empty.is_subgraph_of(&a), "Graph is not a subgrap of itself (relationships in supergraph)");
    assert!(a_empty.is_subgraph_of(&a_empty), "Graph is not a subgrap of itself (no relationships)");
    test_describe_relationship(&a, &a_empty, "original", "same", 0, true);

    // Single step.
    {
        let b_empty = GraphName::new(String::from("B"));
        let mut b = b_empty.clone();
        b.add_relationships(&a);
        assert!(a.is_subgraph_of(&b_empty), "A is not a subgraph of B (relationships in subgraph)");
        assert!(a_empty.is_subgraph_of(&b), "A is not a subgraph of B (relationships in supergraph)");
        assert!(!b.is_subgraph_of(&a), "B is a subgraph of A");
        test_describe_relationship(&a, &b, "subgraph (first)", "supergraph (second)", 1, true);
        test_describe_relationship(&b, &a, "supergraph (first)", "subgraph (second)", 1, true);
    }

    // No path.
    let c_empty = GraphName::new(String::from("C"));
    let mut c = c_empty.clone();
    c.add_relationships(&a);
    assert!(!a.is_subgraph_of(&c), "A is a subgraph of C");
    assert!(!c.is_subgraph_of(&a), "C is a subgraph of A");
    // Skip description test, as there is a path with a translation.

    // Multiple steps.
    {
        let e_empty = GraphName::new(String::from("E"));
        let mut e = e_empty.clone();
        e.add_relationships(&a);
        assert!(c.is_subgraph_of(&e_empty), "C is not a subgraph of E (relationships in subgraph)");
        assert!(c_empty.is_subgraph_of(&e), "C is not a subgraph of E (relationships in supergraph)");
        assert!(!e.is_subgraph_of(&c), "E is a subgraph of C");
        test_describe_relationship(&c, &e, "subgraph", "supergraph", 2, true);
    }

    // Relationships split between graphs.
    {
        let mut from = GraphName::new(String::from("from"));
        from.add_subgraph("from", "middle");
        let mut to = GraphName::new(String::from("to"));
        to.add_subgraph("middle", "to");
        assert!(from.is_subgraph_of(&to), "from is not a subgraph of to");
        assert!(!to.is_subgraph_of(&from), "to is a subgraph of from");
        test_describe_relationship(&from, &to, "from", "to", 2, true);
    }
}

#[test]
fn graph_name_translation() {
    let a = manual();
    let a_empty = GraphName::new(String::from(NAME));

    // Same graph.
    assert!(a.translates_to(&a_empty), "Graph does not translate to itself (relationships in source)");
    assert!(a_empty.translates_to(&a), "Graph does not translate to itself (relationships in target)");
    assert!(a_empty.translates_to(&a_empty), "Graph does not translate to itself (no relationships)");
    test_describe_relationship(&a, &a_empty, "original", "same", 0, true);

    // Single subgraph step.
    let b_empty = GraphName::new(String::from("B"));
    let mut b = b_empty.clone();
    b.add_relationships(&a);
    assert!(a.translates_to(&b_empty), "A does not translate to B (relationships in source)");
    assert!(a_empty.translates_to(&b), "A does not translate to B (relationships in target)");
    assert!(!b.translates_to(&a), "B translates to A");
    test_describe_relationship(&a, &b, "source", "target", 1, true);

    // Single translation step.
    let c_empty = GraphName::new(String::from("C"));
    let mut c = c_empty.clone();
    c.add_relationships(&a);
    assert!(b.translates_to(&c_empty), "B does not translate to C (relationships in source)");
    assert!(b_empty.translates_to(&c), "B does not translate to C (relationships in target)");
    assert!(!c.translates_to(&b), "C translates to B");
    test_describe_relationship(&b, &c, "source", "target", 1, true);

    // Mixed steps.
    assert!(a.translates_to(&c_empty), "A does not translate to C (relationships in source)");
    assert!(a_empty.translates_to(&c), "A does not translate to C (relationships in target)");
    assert!(!c.translates_to(&a), "C translates to A");
    test_describe_relationship(&a, &c, "source", "target", 2, true);

    // Multiple subgraph steps.
    {
        let e_empty = GraphName::new(String::from("E"));
        let mut e = e_empty.clone();
        e.add_relationships(&a);
        assert!(c.translates_to(&e_empty), "C does not translate to E (relationships in source)");
        assert!(c_empty.translates_to(&e), "C does not translate to E (relationships in target)");
        assert!(!e.translates_to(&c), "E translates to C");
        test_describe_relationship(&c, &e, "source", "target", 2, true);
    }

    // Multiple translation steps.
    {
        let f_empty = GraphName::new(String::from("F"));
        let mut f = f_empty.clone();
        f.add_relationships(&a);
        assert!(b.translates_to(&f_empty), "B does not translate to F (relationships in source)");
        assert!(b_empty.translates_to(&f), "B does not translate to F (relationships in target)");
        assert!(!f.translates_to(&b), "F translates to B");
        test_describe_relationship(&b, &f, "source", "target", 2, true);
    }

    // Relationships split between graphs.
    {
        let mut from = GraphName::new(String::from("from"));
        from.add_translation("from", "middle");
        let mut to = GraphName::new(String::from("to"));
        to.add_translation("middle", "to");
        assert!(from.translates_to(&to), "from does not translate to to");
        assert!(!to.translates_to(&from), "to translates to from");
        test_describe_relationship(&from, &to, "source", "target", 2, true);
    }
}

//-----------------------------------------------------------------------------
