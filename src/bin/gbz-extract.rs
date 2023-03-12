use gbwt::{GBZ, Orientation, Metadata};
use gbwt::internal;

use simple_sds::serialize::Serialize;
use simple_sds::serialize;

use std::collections::{HashMap, HashSet};
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;
use std::{env, process};

use getopts::Options;

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
    let metadata = gbz.metadata().unwrap();
    if !metadata.has_path_names() {
        return Err("Sequence extraction requires path names".to_string());
    }

    // Select the paths to extract.
    // FIXME move the restricted case to a function
    let mut selected_paths: Vec<usize> = Vec::new();
    if let Some(contig) = config.contig.as_ref() {
        if config.verbose {
            eprintln!("Finding paths for contig {}", contig);
        }
        // Determine the paths with the given contig name.
        if !metadata.has_contig_names() {
            return Err("Cannot select a contig without contig names".to_string());
        }
        let contig_id = metadata.contig_id(contig);
        if contig_id.is_none() {
            return Err(format!("The graph does not contain contig {}", contig));            
        }
        let contig_id = contig_id.unwrap();
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
    } else {
        selected_paths.extend(0..gbz.paths());
        eprintln!("Selected all {} paths", gbz.paths());
    }

    // FIXME in parallel
    // Extract the paths.
    if config.verbose {
        eprintln!("Extracting paths");
    }
    let mut options = OpenOptions::new();
    options.create(true).write(true).truncate(true);
    let seq_name = config.output.as_ref().unwrap();
    let mut seq_file = options.open(seq_name).map_err(|e| e.to_string())?;
    let pos_name = format!("{}.pos", config.output.as_ref().unwrap());
    let mut pos_file = options.open(pos_name).map_err(|e| e.to_string())?;
    let name_name = format!("{}.names", config.output.as_ref().unwrap());
    let mut name_file = options.open(name_name).map_err(|e| e.to_string())?;
    for path_id in selected_paths {
        let (sequence, positions) = extract_sequence(&gbz, path_id, Orientation::Forward);
        seq_file.write_all(&sequence).map_err(|e| e.to_string())?;
        positions.serialize_body(&mut pos_file).map_err(|e| e.to_string())?;
        let path_name = path_name_as_line(metadata, path_id);
        name_file.write_all(path_name.as_bytes()).map_err(|e| e.to_string())?;
    }
    drop(seq_file);
    drop(pos_file);
    drop(name_file);

    if config.verbose {
        eprintln!("Extracted the sequences in {:.3} seconds", start.elapsed().as_secs_f64());
        internal::report_memory_usage();
        eprintln!("");
    }
    Ok(())
}

//-----------------------------------------------------------------------------

// FIXME add option for reverse complement
// FIXME add terminator character
// FIXME add option for positions
struct Config {
    input: Option<String>,
    output: Option<String>,
    contig: Option<String>,
    threads: usize,
    verbose: bool,
}

impl Config {
    const MIN_THREADS: usize = 1;
    const MAX_THREADS: usize = 64;

    pub fn new() -> Result<Config, String> {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("c", "contig", "restrict to components containing this contig", "STR");
        opts.optopt("o", "output", "base name for output", "FILE");
        opts.optopt("t", "threads", "number of threads for extracting paths (default 1)", "INT");
        opts.optflag("v", "verbose", "print progress information");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let mut config = Config {
            input: None,
            output: None,
            contig: None,
            threads: Self::MIN_THREADS,
            verbose: false,
        };
        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz > graph.gfa", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        if let Some(s) = matches.opt_str("c") {
            config.contig = Some(s);
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

        Ok(config)
    }
}

//-----------------------------------------------------------------------------

fn encode_start(node_id: usize, orientation: Orientation) -> u64 {
    let o_bit: u64 = match orientation {
        Orientation::Forward => 0,
        Orientation::Reverse => 1 << 10,
    };
    (node_id as u64) << 11 + o_bit
}

// FIXME optimize and move to support
fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::with_capacity(sequence.len());
    for &c in sequence.iter().rev() {
        result.push(match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            c => c,
        });
    }
    result
}

fn extract_sequence(gbz: &GBZ, path_id: usize, orientation: Orientation) -> (Vec<u8>, Vec<u64>) {
    let mut sequence: Vec<u8> = Vec::new();
    let mut positions: Vec<u64> = Vec::new();

    for (node_id, node_o) in gbz.path(path_id, orientation).unwrap() {
        let seq = gbz.sequence(node_id).unwrap();
        if node_o == Orientation::Reverse {
            sequence.extend(reverse_complement(seq));
        } else {
            sequence.extend_from_slice(seq);
        }
        let mut pos = encode_start(node_id, node_o);
        for _ in 0..seq.len() {
            positions.push(pos);
            pos += 1;
        }
    }

    // Append terminators.
    sequence.push(0);
    positions.push(0);

    (sequence, positions)
}

fn path_name_as_line(metadata: &Metadata, path_id: usize) -> String {
    let path_name = metadata.path(path_id).unwrap();
    format!("{}\t{}\t{}\t{}\n", metadata.sample_name(path_name.sample()), metadata.contig_name(path_name.contig()), path_name.phase(), path_name.fragment())
}

//-----------------------------------------------------------------------------

