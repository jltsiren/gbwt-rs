# GBWT (in Rust)

This is a Rust reimplementation of parts of the [GBWT](https://github.com/jltsiren/gbwt) and the [GBWTGraph](https://github.com/jltsiren/gbwtgraph).
It is based on the [Simple-SDS](https://github.com/jltsiren/simple-sds) library.

## Scope

### GBWT

- [x] Simple-SDS [file format](https://github.com/jltsiren/gbwt/blob/master/SERIALIZATION.md)
- [x] Iteration over paths
- [x] Unidirectional search
- [x] Bidirectional search
- [x] Metadata
- [ ] Locate queries

### GBWTGraph / GBZ

- [x] Simple-SDS [file format](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md)
- [x] Iteration over nodes and edges
- [x] Iteration over segments and links
- [ ] Iteration over paths and path extensions
- [ ] GFA extraction

### Possible future extensions

* GBWT construction
* GBWT merging
* Cached GBWT

## Notes

* The included `.cargo/config.toml` sets the target CPU to `native`.
* See [Simple-SDS](https://github.com/jltsiren/simple-sds) for assumptions on the environment.
