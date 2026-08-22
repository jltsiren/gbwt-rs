# Releases

## Current version

* Supports `GBZ` version 3 and `GBWT` version 6 with Zstandard compression for the BWT.
* Implements `SerializeVersion` trait for serializing different versions of `GBZ`, `GBWT`, and `Sequences`.
* `GraphName` structure for representing graph names and relationships between graphs (imported from [pggname](https://github.com/jltsiren/pggname)).

## GBZ 0.6.2 (2026-08-11)

* Chain-finding algorithm also records the number of trivial chains with no snarls.

## GBZ 0.6.1 (2026-04-17)

* Command line arguments can use suffixes (e.g. `k`, `MiB`) for sizes and counts that can plausibly be large.
* Top-level chains:
  * `Chains` storing a sketch of the high-level structure of the graph (imported from GBZ-base).
  * Algorithm for finding the chains in a graph.

## GBZ 0.6.0 (2026-03-30)

* GBWT construction:
  * Insert batches of sequences into `MutableGBWT` and convert it to `GBWT`.
  * Insert paths one by one into `GBWTBuilder` to build the `GBWT` in a background thread.
* `Graph` has been renamed to `Sequences` to better reflect the purpose of the structure.
* `gbunzip` extracts actual reference paths in addition to generic paths with option `--paths ref-only`.

## GBZ 0.5.1 (2026-02-13)

* Convert between internal and vg conventions for generic path names when (de)serializing metadata.

## GBZ 0.5.0 (2026-02-12)

* Support for GBZ version 2 with Zstandard compression for sequences.
* Renamed the crate from `gbwt` to `gbz` due to a name conflict in crates.io.

## GBWT-rs 0.4.0 (2025-12-26)

* New functionality:
  * `GBZ::reference_positions` for finding positions on reference paths for indexing by sequence offsets.
  * GBZ methods for setting reference samples and iterating over them.
  * Longest common subsequence algorithms:
    * `lcs()` for two integer sequences.
    * `path_lcs()` for two paths, weighted by sequence lengths.
    * `fast_weighted_lcs()` for two integer sequences, weighted by an arbitrary function.
  * `GraphPosition` for storing a position in the graph.
  * `GBWT` passes through document array samples from the C++ implementation.
    * The C++ implementation no longer has to resample the GBWT if a GBZ file has been modified by this implementation.
* Error type is always `String`, `Utf8Error`, or `io::Error`.

## GBWT-rs 0.3.0 (2024-01-29)

* `gbunzip`:
  * Outputs GFA version 1.1.
  * Option for using [PanSN](https://github.com/pangenome/PanSN-spec) path names.
  * Suppors the reference samples header tag `RS` used by vg.
* Terminology change for compatibility with vg:
  * Paths with sample name `_gbwt_ref` are now generic named paths.
  * Paths can be promoted to reference paths by specifying sample names in GBWT tag `reference_samples`.
* New functionality:
  * `GBZ::weakly_connected_components` for finding weakly connected components in the graph.
  * Functions for determining if an edge or a path is in canonical orientation.

## GBWT-rs 0.2.2 (2022-02-22)

Another patch release for the GBZ paper.

* Increase the maximum number of decompression threads in `gbunzip` from 31 to 64.

## GBWT-rs 0.2.1 (2022-02-17)

Minor patch release for the GBZ paper.

* Support for empty paths.
* `gbunzip` uses Rayon for multithreading.

## GBWT-rs 0.2.0 (2021-11-17)

* `Graph`: GBWTGraph implementation storing node sequences and node-to-segment translation.
* `GBZ`: Graph interface.
  * Iterators in the original graph: nodes, predecessors/successors of a node, paths over nodes, extensions of a set of paths.
  * Iterators in the GFA graph: segments, predecessors/successors of a segment, paths over segments.
* Interface changes:
  * Runs and GBWT positions use `Run` and `Pos` types instead of `(usize, usize)` pairs.
  * Node / path orientation uses `enum Orientation` instead of `bool`.
* Proof-of-concept `gbunzip` tool for decompressing GFA from GBZ.
  * Somewhat slower than the C++ implementation but supports multi-threaded path extraction for faster decompression.
  * Better parallelization when allowed to reorder the paths opportunistically.

## GBWT-rs 0.1.0 (2021-09-16)

This is a reimplementation of some parts of the GBWT and the GBWTGraph in Rust. At the moment, the main purpose is ensuring that the Simple-SDS file format specifications are complete enough. I may use this in future projects, but at the moment everything can change without warning.

The first pre-release includes supports the GBWT Simple-SDS file format as well as path navigation, unidirectional/bidirectional search, and GBWT metadata. The next version will probably support GBWTGraph, the GBZ format, and GFA extraction.

## Release process

* Clean up with `cargo clean`.
* Update version in `Cargo.toml`.
* Switch to crates.io versions of dependencies, if necessary.
* Update `RELEASES.md`.
* Run `cargo clippy --features binaries`.
* Run tests with `cargo test` and `cargo test -- --ignored`.
* Build documentation with `cargo doc`.
* Build the optimized version with `cargo build --release --features binaries`.
* Commit the final changes for the release.
* Publish in crates.io with `cargo publish`.
* Push to GitHub.
* Draft a new release in GitHub.
