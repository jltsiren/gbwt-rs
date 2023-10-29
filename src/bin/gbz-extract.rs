use gbwt::{GBZ, Orientation, Metadata};
use gbwt::{internal, support};

use simple_sds::serialize::Serialize;
use simple_sds::serialize;

use core::slice;
use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write, Read, Seek, SeekFrom};
use std::time::Instant;
use std::{cmp, env, mem, process};

use getopts::Options;
use rayon::prelude::*;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start = Instant::now();
    let config = Config::new().map_err(|x| x.to_string())?;
    rayon::ThreadPoolBuilder::new().num_threads(config.threads).build_global().map_err(|e| e.to_string())?;

    let filename = config.input.as_ref().unwrap();
    if config.verbose {
        eprintln!("Loading GBZ graph {}", filename);
    }
    let gbz: GBZ = serialize::load_from(filename).map_err(|x| x.to_string())?;
    if !gbz.has_metadata() {
        return Err("Sequence extraction requires GBWT metadata".to_string());
    }

    match config.mode {
        Mode::Sequences => {
            extract_sequences(&gbz, &config)?;
        },
        Mode::TagArray => {
            extract_tag_array(&gbz, &config)?;
        },
    }

    if config.verbose {
        eprintln!("Finished in {:.3} seconds", start.elapsed().as_secs_f64());
        internal::report_memory_usage();
        eprintln!("");
    }
    Ok(())
}

//-----------------------------------------------------------------------------

enum Mode {
    Sequences,
    TagArray,
}

// FIXME add option for reverse complement
struct Config {
    input: Option<String>,
    output: Option<String>,
    contig: Option<String>,
    endmarker: u8,
    sa_skip: usize,
    threads: usize,
    mode: Mode,
    verbose: bool,
}

impl Config {
    const MIN_THREADS: usize = 1;
    const MAX_THREADS: usize = 64;

    const BUFFER_SIZE: usize = 1048576;

    pub fn new() -> Result<Config, String> {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        // FIXME add an option for determining if we should skip the first value in SA/BWT
        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("c", "contig", "restrict to components containing this contig", "STR");
        opts.optopt("m", "mode", "extract this data (sequences, tag-array; default sequences)", "STR");
        opts.optopt("o", "output", "base name for output", "FILE");
        opts.optopt("t", "threads", "number of parallel threads (default 1)", "INT");
        opts.optopt("", "endmarker-value", "byte value used to terminate sequences (default 0)", "INT");
        opts.optopt("", "endmarker-char", "character used to terminate sequences", "CHAR");
        opts.optflag("v", "verbose", "print progress information");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let mut config = Config {
            input: None,
            output: None,
            contig: None,
            endmarker: 0,
            sa_skip: 1,
            threads: Self::MIN_THREADS,
            mode: Mode::Sequences,
            verbose: false,
        };
        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        if let Some(s) = matches.opt_str("c") {
            config.contig = Some(s);
        }
        if let Some(s) = matches.opt_str("m") {
            match s.as_str() {
                "sequences" => config.mode = Mode::Sequences,
                "tag-array" => config.mode = Mode::TagArray,
                s => {
                    return Err(format!("Invalid mode: {}", s));
                }
            }
        }
        if let Some(s) = matches.opt_str("o") {
            config.output = Some(s);
        }
        if let Some(s) = matches.opt_str("t") {
            match s.parse::<usize>() {
                Ok(n) => {
                    if !(Self::MIN_THREADS..=Self::MAX_THREADS).contains(&n) {
                        return Err(format!("--threads: number of threads must be between {} and {}", Self::MIN_THREADS, Self::MAX_THREADS));
                    }
                    config.threads = n;
                },
                Err(f) => {
                    return Err(format!("--threads: {}", f.to_string()));
                },
            }
        }
        if let Some(s) = matches.opt_str("endmarker-value") {
            match s.parse::<u8>() {
                Ok(n) => config.endmarker = n,
                Err(f) => return Err(format!("--endmarker-value: {}", f.to_string())),
            }
        }
        if let Some(s) = matches.opt_str("endmarker-char") {
            if s.as_bytes().len() == 1 {
                config.endmarker = s.as_bytes()[0];
            } else {
                return Err(format!("Invalid endmarker character: {}", s));
            }
        }
        if matches.opt_present("v") {
            config.verbose = true;
        }

        if !matches.free.is_empty() {
            config.input = Some(matches.free[0].clone());
        } else {
            let header = format!("Usage: {} [options] -o output graph.gbz", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        }
        if config.output.is_none() {
            return Err("Option -o / --output is mandatory".to_string());
        }

        Ok(config)
    }
}

//-----------------------------------------------------------------------------

fn extract_sequence(gbz: &GBZ, path_id: usize, orientation: Orientation, config: &Config) -> Vec<u8> {
    let mut sequence: Vec<u8> = Vec::new();

    for (node_id, node_o) in gbz.path(path_id, orientation).unwrap() {
        let seq = gbz.sequence(node_id).unwrap();
        if node_o == Orientation::Reverse {
            sequence.extend(support::reverse_complement(seq));
        } else {
            sequence.extend_from_slice(seq);
        }
    }

    // Append the endmarker.
    sequence.push(config.endmarker);

    sequence
}

fn path_name_as_line(metadata: &Metadata, path_id: usize, sequence_len: usize) -> String {
    let path_name = metadata.path(path_id).unwrap();
    format!("{}\t{}\t{}\t{}\t{}\t{}\n", path_id, metadata.sample_name(path_name.sample()), metadata.contig_name(path_name.contig()), path_name.phase(), path_name.fragment(), sequence_len)
}

fn select_paths(gbz: &GBZ, metadata: &Metadata, config: &Config) -> Result<Vec<usize>, String> {
    if config.verbose {
        eprintln!("Selecting paths");
    }

    let mut selected_paths: Vec<usize> = Vec::new();
    if config.contig.is_none() {
        selected_paths.extend(0..metadata.paths());
        if config.verbose {
            eprintln!("Selected all {} paths", metadata.paths());
        }
        return Ok(selected_paths);
    }

    // Determine the paths with the given contig name.
    let contig = config.contig.as_ref().unwrap();
    if !metadata.has_contig_names() {
        return Err("Cannot select a contig without contig names".to_string());
    }
    let contig_id = metadata.contig_id(contig).ok_or(format!("The graph does not contain contig {}", contig))?;
    let mut initial_paths: Vec<usize> = Vec::new();
    for (path_id, path_name) in metadata.path_iter().enumerate() {
        if path_name.contig() == contig_id {
            initial_paths.push(path_id);
        }
    }
    if initial_paths.is_empty() {
        return Err(format!("The graph does not contain any paths for contig {}", contig));
    }
    if config.verbose {
        eprintln!("Found {} paths with contig name {}", initial_paths.len(), contig);
    }

    // Determine weakly connected components in the graph.
    let components = gbz.weakly_connected_components();
    if config.verbose {
        eprintln!("Found {} weakly connected components", components.len());
    }

    // Determine the components containing the initial paths.
    let mut node_to_component: HashMap<usize, usize> = HashMap::with_capacity(gbz.nodes());
    for component_id in 0..components.len() {
        for &node_id in components[component_id].iter() {
            node_to_component.insert(node_id, component_id);
        }
    }
    let mut selected_components: HashSet<usize> = HashSet::new();
    for path_id in initial_paths {
        if let Some((node_id, _)) = gbz.path(path_id, Orientation::Forward).unwrap().next() {
            selected_components.insert(node_to_component[&node_id]);
        }
    }
    if config.verbose {
        eprintln!("Found {} components for contig {}", selected_components.len(), contig);
    }

    // Select all paths in the selected components.
    for path_id in 0..gbz.paths() {
        if let Some((node_id, _)) = gbz.path(path_id, Orientation::Forward).unwrap().next() {
            if selected_components.contains(&node_to_component[&node_id]) {
                selected_paths.push(path_id);
            }
        }
    }
    if config.verbose {
        eprintln!("Found {} paths in the selected components", selected_paths.len());
    }

    Ok(selected_paths)
}

fn extract_sequences(gbz: &GBZ, config: &Config) -> Result<(), String> {
    let metadata = gbz.metadata().ok_or("Sequence extraction requires GBWT metadata".to_string())?;
    if !metadata.has_path_names() {
        return Err("Sequence extraction requires path names".to_string());
    }
    let selected_paths = select_paths(gbz, metadata, config)?;

    if config.verbose {
        eprintln!("Extracting paths");
    }
    let mut options = OpenOptions::new();
    options.create(true).write(true).truncate(true);
    let seq_name = config.output.as_ref().unwrap();
    let mut seq_file = options.open(seq_name).map_err(|e| e.to_string())?;
    let name_name = format!("{}.names", config.output.as_ref().unwrap());
    let mut name_file = options.open(name_name).map_err(|e| e.to_string())?;
    for path_id in selected_paths {
        let sequence = extract_sequence(&gbz, path_id, Orientation::Forward, &config);
        seq_file.write_all(&sequence).map_err(|e| e.to_string())?;
        let path_name = path_name_as_line(metadata, path_id, sequence.len() - 1);
        name_file.write_all(path_name.as_bytes()).map_err(|e| e.to_string())?;
    }
    drop(seq_file);
    drop(name_file);

    Ok(())
}

//-----------------------------------------------------------------------------

// Returns (path id, length, starting offset)
fn read_names(config: &Config) -> Result<Vec<(usize, usize, usize)>, String> {
    let filename = format!("{}.names", config.output.as_ref().unwrap());
    if config.verbose {
        eprintln!("Reading path names from {}.names", config.output.as_ref().unwrap());
    }

    let file = File::open(filename).map_err(|e| e.to_string())?;
    let mut result = Vec::new();
    let mut offset: usize = 0;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|e| e.to_string())?;
        let mut fields = line.split('\t');
        let first = fields.next().ok_or("Name line missing first field".to_string())?;
        let id = first.parse::<usize>().map_err(|e| e.to_string())?;
        let last = fields.next_back().ok_or("Name line missing last field".to_string())?;
        let len = last.parse::<usize>().map_err(|e| e.to_string())?;
        result.push((id, len, offset));
        offset += len + 1;
    }

    Ok(result)
}

// Returns (i, SA[i]).
fn read_suffix_array(expected_len: usize, config: &Config) -> Result<Vec<(usize, u64)>, String> {
    let filename = format!("{}.sa", config.output.as_ref().unwrap());
    if config.verbose {
        eprintln!("Reading suffix array from {}", &filename);
    }
    let mut file = File::open(filename).map_err(|e| e.to_string())?;
    file.seek(SeekFrom::Start((config.sa_skip * mem::size_of::<u64>()) as u64)).map_err(|e| e.to_string())?;

    let mut result: Vec<(usize, u64)> = Vec::with_capacity(expected_len);
    let mut buffer: Vec<u64> = Vec::with_capacity(Config::BUFFER_SIZE);
    let mut offset = 0;
    while result.len() < expected_len {
        let len = cmp::min(expected_len - offset, Config::BUFFER_SIZE);
        unsafe {
            let buf: &mut [u8] = slice::from_raw_parts_mut(buffer.as_mut_ptr() as *mut u8, len * mem::size_of::<u64>());
            file.read_exact(buf).map_err(|e| e.to_string())?;
            buffer.set_len(len);
        }
        result.extend((offset..offset + len).zip(buffer.iter().copied()));
        offset += len;
    }

    Ok(result)
}

fn encode_start(node_id: usize, orientation: Orientation) -> u64 {
    let o_bit: u64 = match orientation {
        Orientation::Forward => 0,
        Orientation::Reverse => 1 << 10,
    };
    (node_id as u64) << 11 + o_bit
}

fn extract_path(gbz: &GBZ, path_id: usize, orientation: Orientation) -> Vec<u64> {
    let mut result: Vec<u64> = Vec::new();

    if let Some(iter) = gbz.path(path_id, orientation) {
        for (node_id, node_o) in iter {
            // FIXME warn if the node is too long
            let len = gbz.sequence_len(node_id).unwrap();
            let mut pos = encode_start(node_id, node_o);
            for _ in 0..len {
                result.push(pos);
                pos += 1;
            }
        }
    }
    result.push(0); // Append the endmarker.

    result
}

fn count_bwt_runs(expected_len: usize, config: &Config) -> Result<(), String> {
    let filename = format!("{}.bwt", config.output.as_ref().unwrap());
    if config.verbose {
        eprintln!("Counting BWT runs in {}", &filename);
    }
    let mut file = File::open(filename).map_err(|e| e.to_string())?;
    file.seek(SeekFrom::Start(config.sa_skip as u64)).map_err(|e| e.to_string())?;

    let mut runs: usize = 0;
    let mut prev: usize = 256;
    let mut buffer: Vec<u8> = Vec::with_capacity(Config::BUFFER_SIZE);
    let mut offset = 0;
    while offset < expected_len {
        let len = cmp::min(expected_len - offset, Config::BUFFER_SIZE);
        unsafe {
            let buf: &mut [u8] = slice::from_raw_parts_mut(buffer.as_mut_ptr() as *mut u8, len);
            file.read_exact(buf).map_err(|e| e.to_string())?;
            buffer.set_len(len);
        }
        for byte in buffer.iter().copied() {
            let val = byte as usize;
            if val != prev {
                runs += 1;
            }
            prev = val;
        }
        offset += len;
    }
    if config.verbose {
        eprintln!("BWT runs: {}", runs);
    }

    Ok(())
}

fn extract_tag_array(gbz: &GBZ, config: &Config) -> Result<(), String> {
    let names = read_names(config)?;
    let last = names.last().ok_or("No path names found".to_string())?;
    let expected_len = last.1 + last.2 + 1;
    if config.verbose {
        eprintln!("Found {} paths; expected file size {}", names.len(), expected_len);
    }

    // Read the suffix array into (i, SA[i]).
    let mut values = read_suffix_array(expected_len, config)?;

    // Sort by second field to get (ISA[i], i).
    if config.verbose {
        eprintln!("Converting suffix array to inverse suffix array");
    }
    values.par_sort_unstable_by(|a, b| a.1.cmp(&b.1));

    // FIXME in parallel?
    // Extract paths and replace the second field with tags.
    if config.verbose {
        eprintln!("Extracting paths");
    }
    for (path_id, len, offset) in names {
        let path = extract_path(gbz, path_id, Orientation::Forward);
        if path.len() != len + 1 {
            return Err(format!("Invalid length for path {}: expected {}, got {}", path_id, len, path.len() - 1));
        }
        for i in 0..path.len() {
            values[offset + i].1 = path[i];
        }
    }

    // Sort by first field to get (i, TAG[i]).
    if config.verbose {
        eprintln!("Sorting the tag array");
    }
    values.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));

    // Write the tag array to file.
    let tag_name = format!("{}.tags", config.output.as_ref().unwrap());
    if config.verbose {
        eprintln!("Writing the tag array to {}", tag_name);
    }
    let mut options = OpenOptions::new();
    options.create(true).write(true).truncate(true);
    let mut tag_file = options.open(tag_name).map_err(|e| e.to_string())?;
    let mut buffer: Vec<u64> = Vec::with_capacity(Config::BUFFER_SIZE);
    let mut prev = u64::MAX;
    let mut runs: usize = 0;
    for (_, tag) in values {
        if tag != prev {
            runs += 1;
        }
        prev = tag;
        if buffer.len() >= Config::BUFFER_SIZE {
            buffer.serialize_body(&mut tag_file).map_err(|e| e.to_string())?;
            buffer.clear();
        }
        buffer.push(tag);
    }
    if !buffer.is_empty() {
        buffer.serialize_body(&mut tag_file).map_err(|e| e.to_string())?;
        buffer.clear();
    }
    if config.verbose {
        eprintln!("Tag array runs: {}", runs);
    }

    // Finally count BWT runs.
    if config.verbose {
        count_bwt_runs(expected_len, config)?;
    }

    Ok(())
}

//-----------------------------------------------------------------------------

