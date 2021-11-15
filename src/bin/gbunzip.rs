use gbwt::{GBZ, Orientation};
use gbwt::REF_SAMPLE;
use gbwt::internal;

use simple_sds::serialize::Serialize;
use simple_sds::serialize;

use std::io::{Write, BufWriter};
use std::sync::Arc;
use std::thread::JoinHandle;
use std::time::Instant;
use std::{env, io, process, thread};

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
    let gbz = Arc::new(gbz);
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

    write_gfa(&gbz, &config).map_err(|x| x.to_string())?;

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
    pub threads: usize,
    pub buffer_size: usize,
    pub verbose: bool,
}

impl Config {
    const MIN_THREADS: usize = 1;
    const MAX_THREADS: usize = 31;
    const BUFFER_SIZE: usize = 8 * 1048576;

    pub fn new() -> Result<Config, String> {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optopt("b", "buffer-size", "output buffer size in megabytes (default 8)", "INT");
        opts.optflag("h", "help", "print this help");
        opts.optopt("t", "threads", "number of threads for extracting paths (default 1)", "INT");
        opts.optflag("v", "verbose", "write progress information");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let mut config = Config {
            filename: None,
            threads: Self::MIN_THREADS,
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

fn write_gfa(gbz: &Arc<GBZ>, config: &Config) -> io::Result<()> {
    let stdout = io::stdout();
    let mut buffer = BufWriter::with_capacity(config.buffer_size, stdout.lock());
    buffer.write_all(b"H\tVN:Z:1.0\n")?;
    write_segments(gbz, &mut buffer, config)?;
    write_links(gbz, &mut buffer, config)?;
    buffer.flush()?;

    let stdout = io::stdout();
    let mut handle = stdout.lock();
    write_paths(gbz, &mut handle, config)?;
    write_walks(gbz, &mut handle, config)?;

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

fn write_paths<T: Write>(gbz: &Arc<GBZ>, output: &mut T, config: &Config) -> io::Result<()> {
    let start = Instant::now();
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

    // Determine the identifiers of reference paths.
    let mut paths: Vec<usize> = Vec::new();
    for (path_id, path_name) in metadata.path_iter().enumerate() {
        if path_name.sample() == ref_sample {
            paths.push(path_id);
        }
    }

    // Extract the paths in order using multiple threads.
    let mut thread_pool = create_thread_pool(config.threads);
    let mut next: usize = 0;
    for &path_id in paths.iter() {
        if thread_pool[next].is_some() {
            join_thread(&mut thread_pool, next, output)?;
        }
        let private_ref = Arc::clone(&gbz);
        let handle = thread::spawn(move || {
            path_to_p_line(&private_ref, path_id)
        });
        thread_pool[next] = Some(handle);
        next = (next + 1) % thread_pool.len();
    }
    join_thread_pool(&mut thread_pool, next, output)?;

    if config.verbose {
        eprintln!("Wrote {} paths in {:.3} seconds", paths.len(), start.elapsed().as_secs_f64());
    }
    Ok(())
}

fn path_to_p_line(gbz: &Arc<GBZ>, path_id: usize) -> Vec<u8> {
    let mut len = 0;
    let mut buffer: Vec<u8> = Vec::new();

    buffer.push(b'P');
    buffer.push(b'\t');
    let metadata = gbz.metadata().unwrap();
    let path_name = metadata.path(path_id).unwrap();
    let contig_name = metadata.contig(path_name.contig()).unwrap();
    buffer.extend_from_slice(contig_name.as_bytes());
    buffer.push(b'\t');

    match gbz.segment_path(path_id, Orientation::Forward) {
        Some(iter) => {
            for (segment, orientation) in iter {
                if len > 0 {
                    buffer.push(b',');
                }
                buffer.extend_from_slice(segment.name);
                match orientation {
                    Orientation::Forward => buffer.push(b'+'),
                    Orientation::Reverse => buffer.push(b'-'),
                }
                len += 1;
            }
        },
        None => {
            for (node_id, orientation) in gbz.path(path_id, Orientation::Forward).unwrap() {
                if len > 0 {
                    buffer.push(b',');
                }
                buffer.extend_from_slice(node_id.to_string().as_bytes());
                match orientation {
                    Orientation::Forward => buffer.push(b'+'),
                    Orientation::Reverse => buffer.push(b'-'),
                }
                len += 1;
            }
        },
    }

    buffer.push(b'\t');
    if len > 1 {
        buffer.push(b'*');
        for _ in 2..len {
            buffer.push(b',');
            buffer.push(b'*');
        }
    }
    buffer.push(b'\n');
    buffer
}

//-----------------------------------------------------------------------------

fn write_walks<T: Write>(gbz: &Arc<GBZ>, output: &mut T, config: &Config) -> io::Result<()> {
    let start = Instant::now();
    let metadata = gbz.metadata().unwrap();
    let ref_sample = metadata.sample_id(REF_SAMPLE).unwrap_or(metadata.samples());
    if config.verbose {
        eprintln!("Writing walks");
    }

    // Determine the identifiers of non-reference paths.
    let mut paths: Vec<usize> = Vec::new();
    for (path_id, path_name) in metadata.path_iter().enumerate() {
        if path_name.sample() != ref_sample {
            paths.push(path_id);
        }
    }

    // Extract the paths in order using multiple threads.
    let mut thread_pool = create_thread_pool(config.threads);
    let mut next: usize = 0;
    for &path_id in paths.iter() {
        if thread_pool[next].is_some() {
            join_thread(&mut thread_pool, next, output)?;
        }
        let private_ref = Arc::clone(&gbz);
        let handle = thread::spawn(move || {
            path_to_w_line(&private_ref, path_id)
        });
        thread_pool[next] = Some(handle);
        next = (next + 1) % thread_pool.len();
    }
    join_thread_pool(&mut thread_pool, next, output)?;

    if config.verbose {
        eprintln!("Wrote {} walks in {:.3} seconds", paths.len(), start.elapsed().as_secs_f64());
    }
    Ok(())
}

fn path_to_w_line(gbz: &Arc<GBZ>, path_id: usize) -> Vec<u8> {
    let mut buffer: Vec<u8> = Vec::new();

    let metadata = gbz.metadata().unwrap();
    let path_name = metadata.path(path_id).unwrap();
    buffer.push(b'W');
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.sample(path_name.sample()).unwrap().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(path_name.phase().to_string().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.contig(path_name.contig()).unwrap().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(path_name.fragment().to_string().as_bytes());
    buffer.push(b'\t');

    match gbz.segment_path(path_id, Orientation::Forward) {
        Some(iter) => {
            let mut path: Vec<(&[u8], Orientation)> = Vec::new();
            let mut len: usize = 0;
            for (segment, orientation) in iter {
                path.push((segment.name, orientation));
                len += segment.sequence.len();
            }
            buffer.extend_from_slice((path_name.fragment() + len).to_string().as_bytes());
            buffer.push(b'\t');
            for (name, orientation) in path.iter() {
                match orientation {
                    Orientation::Forward => buffer.push(b'>'),
                    Orientation::Reverse => buffer.push(b'<'),
                }
                buffer.extend_from_slice(name);
            }
        },
        None => {
            let mut path: Vec<(usize, Orientation)> = Vec::new();
            let mut len: usize = 0;
            for (node_id, orientation) in gbz.path(path_id, Orientation::Forward).unwrap() {
                path.push((node_id, orientation));
                len += gbz.sequence_len(node_id).unwrap_or(0);
            }
            buffer.extend_from_slice((path_name.fragment() + len).to_string().as_bytes());
            buffer.push(b'\t');
            for (node_id, orientation) in path.iter() {
                match orientation {
                    Orientation::Forward => buffer.push(b'>'),
                    Orientation::Reverse => buffer.push(b'<'),
                }
                buffer.extend_from_slice(node_id.to_string().as_bytes());
            }
        },
    }

    buffer.push(b'\n');
    buffer
}

//-----------------------------------------------------------------------------

fn create_thread_pool(threads: usize) -> Vec<Option<JoinHandle<Vec<u8>>>> {
    let mut thread_pool = Vec::new();
    for _ in 0..threads {
        thread_pool.push(None);
    }
    thread_pool
}

fn join_thread<T: Write>(thread_pool: &mut Vec<Option<JoinHandle<Vec<u8>>>>, next: usize, output: &mut T) -> io::Result<()> {
    thread_pool.push(None);
    let line = thread_pool.swap_remove(next).unwrap().join().unwrap();
    output.write_all(&line)
}

fn join_thread_pool<T: Write>(thread_pool: &mut Vec<Option<JoinHandle<Vec<u8>>>>, next: usize, output: &mut T) -> io::Result<()> {
    let mut next = next;
    while thread_pool[next].is_some() {
        join_thread(thread_pool, next, output)?;
        next = (next + 1) % thread_pool.len();
    }
    Ok(())
}

//-----------------------------------------------------------------------------
