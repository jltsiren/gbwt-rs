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
