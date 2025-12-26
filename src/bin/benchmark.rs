#![allow(
    clippy::uninlined_format_args,
    clippy::new_without_default
)]

use gbwt::{GBWT, Pos};
use gbwt::internal;

use simple_sds::serialize::Serialize;
use simple_sds::serialize;

use std::time::Instant;
use std::{env, process};

use getopts::Options;
use rand::Rng;

//-----------------------------------------------------------------------------

fn main() {
    let config = Config::new();

    let filename = config.filename.as_ref().unwrap();
    eprintln!("Loading GBWT index {}", filename);
    let index: GBWT = serialize::load_from(filename).unwrap();
    let (size, units) = internal::readable_size(index.size_in_bytes());
    eprintln!("Index size: {:.3} {}", size, units);
    if index.is_empty() {
        eprintln!("Cannot perform benchmarks with an empty index");
        process::exit(1);
    }
    eprintln!();

    let queries = generate_queries(&index, &config);
    unidirectional_search(&index, &queries);

    internal::report_memory_usage();
    eprintln!();
}

//-----------------------------------------------------------------------------

pub struct Config {
    pub filename: Option<String>,
    pub queries: usize,
    pub query_len: usize,
}

impl Config {
    const QUERIES: usize = 1000000;
    const QUERY_LEN: usize = 10;

    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("n", "queries", "number of queries (default 1000000)", "INT");
        opts.optopt("l", "query-len", "query length (default 10)", "INT");
        let matches = match opts.parse(&args[1..]) {
            Ok(m) => m,
            Err(f) => {
                eprintln!("{}", f);
                process::exit(1);
            }
        };

        let mut config = Config {
            filename: None,
            queries: Self::QUERIES,
            query_len: Self::QUERY_LEN,
        };
        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbwt", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        if let Some(s) = matches.opt_str("n") {
            match s.parse::<usize>() {
                Ok(n) => {
                    if n == 0 {
                        eprintln!("--queries: number of queries must be non-zero");
                        process::exit(1);
                    }
                    config.queries = n;
                },
                Err(f) => {
                    eprintln!("--queries: {}", f);
                    process::exit(1);
                },
            }
        }
        if let Some(s) = matches.opt_str("l") {
            match s.parse::<usize>() {
                Ok(n) => {
                    if n == 0 {
                        eprintln!("--query-len: query length must be non-zero");
                        process::exit(1);
                    }
                    config.query_len = n;
                },
                Err(f) => {
                    eprintln!("--query-len: {}", f);
                    process::exit(1);
                },
            }
        }

        if !matches.free.is_empty() {
            config.filename = Some(matches.free[0].clone());
        } else {
            let header = format!("Usage: {} [options] graph.gbwt", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        }

        config
    }
}

//-----------------------------------------------------------------------------

fn generate_queries(index: &GBWT, config: &Config) -> Vec<Vec<usize>> {
    println!("Generating {} queries of length {}", config.queries, config.query_len);
    let mut queries: Vec<Vec<usize>> = Vec::new();
    let mut rng = rand::rng();

    while queries.len() < config.queries {
        let mut query: Vec<usize> = Vec::new();
        let mut curr = Pos::new(rng.random_range(index.first_node()..index.alphabet_size()), 0);
        if let Some(state) = index.find(curr.node) {
            curr.offset = rng.random_range(0..state.len());
        } else {
            continue;
        }
        query.push(curr.node);
        while query.len() < config.query_len {
            if let Some(next) = index.forward(curr) {
                query.push(next.node);
                curr = next;
            } else {
                break;
            }
        }
        if query.len() == config.query_len {
            queries.push(query);
        }
    }

    println!();
    queries
}

fn unidirectional_search(index: &GBWT, queries: &[Vec<usize>]) {
    println!("Running {} unidirectional queries", queries.len());
    let now = Instant::now();
    let mut total_len = 0;
    let mut total_occs = 0;
    for query in queries {
        let mut state = index.find(query[0]).unwrap();
        for node in query.iter().skip(1) {
            state = index.extend(&state, *node).unwrap();
        }
        total_len += query.len();
        total_occs += state.len();
    }
    internal::report_results(queries.len(), total_len, total_occs, now.elapsed());
}

//-----------------------------------------------------------------------------
