# GBWT (in Rust)

This is a Rust reimplementation of parts of the [GBWT](https://github.com/jltsiren/gbwt) and the [GBWTGraph](https://github.com/jltsiren/gbwtgraph).
It is based on the [Simple-SDS](https://github.com/jltsiren/simple-sds) library.

## Scope

The initial goal is to support the following functionality:

* Simple-SDS file formats for [GBWT](https://github.com/jltsiren/gbwt/SERIALIZATION.md) and [GBZ](https://github.com/jltsiren/gbwtgraph/SERIALIZATION.md).
* Unidirectional and bidirectional search in GBWT.
* Iteration over paths.
* GFA extraction from a GBZ file.

Various construction/merging algorithms and extensions such as the caching layer and document listing structures may be implemented later.

## Notes

* The included `.cargo/config.toml` sets the target CPU to `native`.
* See [Simple-SDS](https://github.com/jltsiren/simple-sds) for assumptions on the environment.
