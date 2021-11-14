use gbwt::{GBZ, Orientation, PathName};
use gbwt::REF_SAMPLE;
use gbwt::internal;

use simple_sds::serialize::Serialize;
use simple_sds::serialize;

use std::io::{Write, BufWriter};
use std::time::Instant;
use std::{env, io, process};

use getopts::Options;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start = Instant::now();
    let config = Config::new().map_err(|x| x.to_string())?;

    let filename = config.filename.as_ref().unwrap();
    if config.verbose {
        eprintln!("Loading GBZ graph {}", filename);
    }
    let gbz: GBZ = serialize::load_from(filename).map_err(|x| x.to_string())?;
    if !gbz.has_metadata() {
        return Err("GFA decompression requires GBWT metadata".to_string());
    }
    if let Some(metadata) = gbz.metadata() {
        if !metadata.has_sample_names() || !metadata.has_contig_names() || !metadata.has_path_names() {
            return Err("GFA decompression requires sample / contig / path names".to_string());
        }
    }
    if config.verbose {
        let (size, units) = internal::readable_size(gbz.size_in_bytes());
        eprintln!("GBZ size: {:.3} {}", size, units);
        eprintln!("");
    }

    let stdout = io::stdout();
    let mut buffer = BufWriter::with_capacity(config.buffer_size, stdout.lock());
    write_gfa(&gbz, &mut buffer, &config).map_err(|x| x.to_string())?;

    if config.verbose {
        eprintln!("");
        eprintln!("GFA decompressed in {:.3} seconds", start.elapsed().as_secs_f64());
        internal::report_memory_usage();
        eprintln!("");
    }
    Ok(())
}

//-----------------------------------------------------------------------------

pub struct Config {
    pub filename: Option<String>,
    pub buffer_size: usize,
    pub verbose: bool,
}

impl Config {
    const BUFFER_SIZE: usize = 8 * 1048576;

    pub fn new() -> Result<Config, String> {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optopt("b", "buffer-size", "output buffer size in megabytes (default 8)", "INT");
        opts.optflag("h", "help", "print this help");
        opts.optflag("v", "verbose", "write progress information");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let mut config = Config {
            filename: None,
            buffer_size: Self::BUFFER_SIZE,
            verbose: false,
        };
        if let Some(s) = matches.opt_str("b") {
            match s.parse::<usize>() {
                Ok(n) => {
                    if n == 0 {
                        return Err("--buffer-size: buffer size must be > 0".to_string());
                    }
                    config.buffer_size = n * 1048576;
                },
                Err(f) => {
                    return Err(format!("--buffer-size: {}", f.to_string()));
                },
            }
        }
        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz > graph.gfa", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        if matches.opt_present("v") {
            config.verbose = true;
        }

        if !matches.free.is_empty() {
            config.filename = Some(matches.free[0].clone());
        } else {
            let header = format!("Usage: {} [options] graph.gbz > graph.gfa", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        }

        Ok(config)
    }
}

//-----------------------------------------------------------------------------

fn write_gfa<T: Write>(gbz: &GBZ, output: &mut T, config: &Config) -> io::Result<()> {
    output.write_all(b"H\tVN:Z:1.0\n")?;
    write_segments(gbz, output, config)?;
    write_links(gbz, output, config)?;
    write_paths(gbz, output, config)?;
    write_walks(gbz, output, config)?;
    Ok(())
}

//-----------------------------------------------------------------------------

fn write_segments<T: Write>(gbz: &GBZ, output: &mut T, config: &Config) -> io::Result<()> {
    let start = Instant::now();
    let mut segments = 0;
    if config.verbose {
        eprintln!("Writing segments");
    }

    match gbz.segment_iter() {
        Some(iter) => {
            for segment in iter {
                write_segment(segment.name, segment.sequence, output)?;
                segments += 1;
            }
        },
        None => {
            for node_id in gbz.node_iter() {
                write_segment(node_id.to_string().as_bytes(), gbz.sequence(node_id).unwrap(), output)?;
                segments += 1;
            }
        },
    }

    if config.verbose {
        eprintln!("Wrote {} segments in {:.3} seconds", segments, start.elapsed().as_secs_f64());
    }
    Ok(())
}

fn write_segment<T: Write>(name: &[u8], sequence: &[u8], output: &mut T) -> io::Result<()> {
    output.write_all(b"S\t")?;
    output.write_all(name)?;
    output.write_all(b"\t")?;
    output.write_all(sequence)?;
    output.write_all(b"\n")
}

//-----------------------------------------------------------------------------

fn write_links<T: Write>(gbz: &GBZ, output: &mut T, config: &Config) -> io::Result<()> {
    let start = Instant::now();
    let mut links = 0;
    if config.verbose {
        eprintln!("Writing links");
    }

    match gbz.segment_iter() {
        Some(iter) => {
            for segment in iter {
                for (successor, orientation) in gbz.segment_successors(&segment, Orientation::Forward).unwrap() {
                    // A link from forward orientation is canonical if it is a self-loop or the successor
                    // has a greater id.
                    if successor.id >= segment.id {
                        write_link((segment.name, Orientation::Forward), (successor.name, orientation), output)?;
                        links += 1;
                    }
                }
                for (successor, orientation) in gbz.segment_successors(&segment, Orientation::Reverse).unwrap() {
                    // A link from reverse orientation is canonical if it is a self-loop that to forward
                    // orientation or the successor has a greater id.
                    if successor.id > segment.id || (successor.id == segment.id && orientation == Orientation::Forward) {
                        write_link((segment.name, Orientation::Reverse), (successor.name, orientation), output)?;
                        links += 1;
                    }
                }
            }
        },
        None => {
            for node_id in gbz.node_iter() {
                for (successor, orientation) in gbz.successors(node_id, Orientation::Forward).unwrap() {
                    if successor >= node_id {
                        write_link((node_id.to_string().as_bytes(), Orientation::Forward), (successor.to_string().as_bytes(), orientation), output)?;
                        links += 1;
                    }
                }
                for (successor, orientation) in gbz.successors(node_id, Orientation::Reverse).unwrap() {
                    if successor > node_id || (successor == node_id && orientation == Orientation::Forward) {
                        write_link((node_id.to_string().as_bytes(), Orientation::Reverse), (successor.to_string().as_bytes(), orientation), output)?;
                        links += 1;
                    }
                }
            }
        },
    }

    if config.verbose {
        eprintln!("Wrote {} links in {:.3} seconds", links, start.elapsed().as_secs_f64());
    }
    Ok(())
}

fn write_link<T: Write>(from: (&[u8], Orientation), to: (&[u8], Orientation), output: &mut T) -> io::Result<()> {
    output.write_all(b"L\t")?;
    output.write_all(from.0)?;
    match from.1 {
        Orientation::Forward => output.write_all(b"\t+\t")?,
        Orientation::Reverse => output.write_all(b"\t-\t")?,
    }
    output.write_all(to.0)?;
    match to.1 {
        Orientation::Forward => output.write_all(b"\t+\t*\n"),
        Orientation::Reverse => output.write_all(b"\t-\t*\n"),
    }
}

//-----------------------------------------------------------------------------

fn write_paths<T: Write>(gbz: &GBZ, output: &mut T, config: &Config) -> io::Result<()> {
    let start = Instant::now();
    let mut paths = 0;
    let metadata = gbz.metadata().unwrap();
    let ref_sample = metadata.sample_id(REF_SAMPLE);
    if ref_sample.is_none() {
        eprintln!("No reference paths in the graph");
        return Ok(());
    }
    let ref_sample = ref_sample.unwrap();
    if config.verbose {
        eprintln!("Writing paths");
    }

    for (path_id, path_name) in metadata.path_iter().enumerate() {
        if path_name.sample() == ref_sample {
            write_path(gbz, path_id, metadata.contig(path_name.contig()).unwrap(), output)?;
            paths += 1;
        }
    }

    if config.verbose {
        eprintln!("Wrote {} paths in {:.3} seconds", paths, start.elapsed().as_secs_f64());
    }
    Ok(())
}

fn write_path<T: Write>(gbz: &GBZ, path_id: usize, path_name: &str, output: &mut T) -> io::Result<()> {
    let mut len = 0;
    output.write_all(b"P\t")?;
    output.write_all(path_name.as_bytes())?;
    output.write_all(b"\t")?;

    match gbz.segment_path(path_id, Orientation::Forward) {
        Some(iter) => {
            for (segment, orientation) in iter {
                if len > 0 {
                    output.write_all(b",")?;
                }
                output.write_all(segment.name)?;
                match orientation {
                    Orientation::Forward => output.write_all(b"+")?,
                    Orientation::Reverse => output.write_all(b"-")?,
                }
                len += 1;
            }
        },
        None => {
            for (node_id, orientation) in gbz.path(path_id, Orientation::Forward).unwrap() {
                if len > 0 {
                    output.write_all(b",")?;
                }
                output.write_all(node_id.to_string().as_bytes())?;
                match orientation {
                    Orientation::Forward => output.write_all(b"+")?,
                    Orientation::Reverse => output.write_all(b"-")?,
                }
                len += 1;
            }
        },
    }

    output.write_all(b"\t")?;
    if len > 1 {
        output.write_all(b"*")?;
        for _ in 2..len {
            output.write_all(b",*")?;
        }
    }
    output.write_all(b"\n")?;
    Ok(())
}

//-----------------------------------------------------------------------------

fn write_walks<T: Write>(gbz: &GBZ, output: &mut T, config: &Config) -> io::Result<()> {
    let start = Instant::now();
    let mut walks = 0;
    let metadata = gbz.metadata().unwrap();
    let ref_sample = metadata.sample_id(REF_SAMPLE).unwrap_or(metadata.samples());
    if config.verbose {
        eprintln!("Writing walks");
    }

    for (path_id, path_name) in metadata.path_iter().enumerate() {
        if path_name.sample() != ref_sample {
            write_walk(gbz, path_id, *path_name, output)?;
            walks += 1;
        }
    }

    if config.verbose {
        eprintln!("Wrote {} walks in {:.3} seconds", walks, start.elapsed().as_secs_f64());
    }
    Ok(())
}

fn write_walk<T: Write>(gbz: &GBZ, path_id: usize, path_name: PathName, output: &mut T) -> io::Result<()> {
    let metadata = gbz.metadata().unwrap();
    output.write_all(b"W\t")?;
    output.write_all(metadata.sample(path_name.sample()).unwrap().as_bytes())?;
    output.write_all(b"\t")?;
    output.write_all(path_name.phase().to_string().as_bytes())?;
    output.write_all(b"\t")?;
    output.write_all(metadata.contig(path_name.contig()).unwrap().as_bytes())?;
    output.write_all(b"\t")?;
    output.write_all(path_name.fragment().to_string().as_bytes())?;
    output.write_all(b"\t")?;

    match gbz.segment_path(path_id, Orientation::Forward) {
        Some(iter) => {
            let mut path: Vec<(&[u8], Orientation)> = Vec::new();
            let mut len: usize = 0;
            for (segment, orientation) in iter {
                path.push((segment.name, orientation));
                len += segment.sequence.len();
            }
            output.write_all((path_name.fragment() + len).to_string().as_bytes())?;
            output.write_all(b"\t")?;
            for (name, orientation) in path.iter() {
                match orientation {
                    Orientation::Forward => output.write_all(b">")?,
                    Orientation::Reverse => output.write_all(b"<")?,
                }
                output.write_all(name)?;
            }
        },
        None => {
            let mut path: Vec<(usize, Orientation)> = Vec::new();
            let mut len: usize = 0;
            for (node_id, orientation) in gbz.path(path_id, Orientation::Forward).unwrap() {
                path.push((node_id, orientation));
                len += gbz.sequence_len(node_id).unwrap_or(0);
            }
            output.write_all((path_name.fragment() + len).to_string().as_bytes())?;
            output.write_all(b"\t")?;
            for (node_id, orientation) in path.iter() {
                match orientation {
                    Orientation::Forward => output.write_all(b">")?,
                    Orientation::Reverse => output.write_all(b"<")?,
                }
                output.write_all(node_id.to_string().as_bytes())?;
            }
        },
    }

    output.write_all(b"\n")?;
    Ok(())
}

//-----------------------------------------------------------------------------