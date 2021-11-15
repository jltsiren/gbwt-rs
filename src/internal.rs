use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread::JoinHandle;
use std::time::Duration;

//-----------------------------------------------------------------------------

pub fn readable_size(bytes: usize) -> (f64, &'static str) {
    let units: Vec<(f64, &'static str)> = vec![
        (1.0, "B"),
        (1024.0, "KiB"),
        (1024.0 * 1024.0, "MiB"),
        (1024.0 * 1024.0 * 1024.0, "GiB"),
        (1024.0 * 1024.0 * 1024.0 * 1024.0, "TiB"),
    ];

    let value = bytes as f64;
    let mut unit = 0;
    for i in 1..units.len() {
        if value >= units[i].0 {
            unit = i;
        } else {
            break;
        }
    }

    (value / units[unit].0, units[unit].1)
}

#[cfg(target_os = "linux")]
pub fn peak_memory_usage() -> Result<usize, &'static str> {
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
        match retval {
            0 => Ok(rusage.ru_maxrss as usize * 1024),
            _ => Err("libc::getrusage call failed"),
        }
    }
}

#[cfg(target_os = "macos")]
pub fn peak_memory_usage() -> Result<usize, &'static str> {
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
        match retval {
            0 => Ok(rusage.ru_maxrss as usize),
            _ => Err("libc::getrusage call failed"),
        }
    }
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
pub fn peak_memory_usage() -> Result<usize, &'static str> {
    Err("No peak_memory_usage implementation for this OS")
}

//-----------------------------------------------------------------------------

pub fn report_results(queries: usize, total_len: usize, total_occs: usize, duration: Duration) {
    let us = (duration.as_micros() as f64) / (queries as f64);
    let ns = (duration.as_nanos() as f64) / (total_len as f64);
    let occs = (total_occs as f64) / (queries as f64);
    eprintln!("Time:        {:.3} seconds ({:.3} us/query, {:.1} ns/node)", duration.as_secs_f64(), us, ns);
    eprintln!("Occurrences: {} total ({:.3} per query)", total_occs, occs);
    eprintln!("");
}

pub fn report_memory_usage() {
    match peak_memory_usage() {
        Ok(bytes) => {
            let (size, unit) = readable_size(bytes);
            eprintln!("Peak memory usage: {:.3} {}", size, unit);
        },
        Err(f) => {
            eprintln!("{}", f);
        },
    }
}

//-----------------------------------------------------------------------------

pub struct ThreadPool<T> {
    threads: Vec<Option<JoinHandle<T>>>,
    // A thread can set its signal to `true` to indicate that it is about to finish.
    signals: Vec<Arc<AtomicBool>>,
}

impl<T> ThreadPool<T> {
    pub fn new(num_threads: usize) -> Self {
        let mut threads = Vec::new();
        let mut signals = Vec::new();
        for _ in 0..num_threads {
            threads.push(None);
            signals.push(Arc::new(AtomicBool::new(false)));
        }
        ThreadPool {
            threads: threads,
            signals: signals,
        }
    }

    pub fn len(&self) -> usize {
        self.threads.len()
    }

    // Returns a copy of the signal for the given thread.
    pub fn signal(&self, thread_id: usize) -> Arc<AtomicBool> {
        self.signals[thread_id].clone()
    }

    // Joins the given thread if it is running and returns the result.
    // Also resets the signal for the given thread.
    pub fn join(&mut self, thread_id: usize) -> Option<T> {
        if self.threads[thread_id].is_some() {
            self.threads.push(None);
            self.signals[thread_id].store(false, Ordering::SeqCst);
            self.threads.swap_remove(thread_id).unwrap().join().ok()
        } else {
            None
        }
    }

    // Returns an unused thread id if one is available. May join an existing
    // thread that has signaled it is about to finish. In that case, also
    // returns the result from that thread.
    pub fn try_join(&mut self) -> Option<(usize, Option<T>)> {
        for thread_id in 0..self.len() {
            if self.threads[thread_id].is_none() || self.signals[thread_id].load(Ordering::SeqCst) {
                return Some((thread_id, self.join(thread_id)));
            }
        }
        None
    }

    // Joins all threads and returns the results, starting from the given id.
    pub fn join_all(&mut self, thread_id: usize) -> Vec<T> {
        let mut results = Vec::new();
        for offset in 0..self.len() {
            let candidate = (thread_id + offset) % self.len();
            if let Some(result) = self.join(candidate) {
                results.push(result);
            }
        }
        results
    }

    // Inserts the given thread into the pool with the given id.
    // Requires that the identifier is unused.
    pub fn insert(&mut self, handle: JoinHandle<T>, thread_id: usize) {
        assert!(self.threads[thread_id].is_none(), "Thread id {} is already in use", thread_id);
        self.threads[thread_id] = Some(handle);
    }
}

//-----------------------------------------------------------------------------
