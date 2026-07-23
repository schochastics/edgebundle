# Changelog

## edgebundle (development version)

- switched the test suite to testthat edition 3 and added correctness
  tests for every exported function
- modernised all `igraph` calls to the 2.x API, removing deprecation
  warnings
- [`edge_bundle_force()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_force.md)
  now also accepts a two-column edgelist matrix
- fixed
  [`edge_bundle_stub()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_stub.md):
  the angular grouping was broken (wrong circular-gap handling and a
  bundle-size cap that summed cluster ids instead of counting edges) and
  vertical edges could produce `NaN`. Bundling output changes as a
  result
- fixed
  [`edge_bundle_path()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_path.md):
  the path-length used in the distortion test was computed along the
  wrong vertices, so edges were essentially never routed. Bundling
  output changes as a result
- added divided edge bundling for directed graphs via
  `edge_bundle_force(directed = TRUE)` (Selassie et al. 2011): edges
  running in opposite directions are kept in separate lanes
- [`edge_bundle_hammer()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_hammer.md)
  is now a native C++ implementation of KDE edge bundling (Hurter et
  al. 2012) and no longer depends on Python/`reticulate`/datashader.
  `reticulate` was dropped from Imports and `install_bundle_py()` was
  removed. Bundling output differs from the datashader-based version
- added
  [`edge_bundle_mingle()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_mingle.md),
  multilevel agglomerative edge bundling (Gansner et al. 2011)
- [`metro_multicriteria()`](https://schochastics.github.io/edgebundle/reference/metro_multicriteria.md)
  is deprecated in favour of `graphlayouts::layout_as_metromap()` and
  will be removed in a future release

## edgebundle 0.4.2

CRAN release: 2023-12-16

- clean up codebase

## edgebundle 0.4.1

CRAN release: 2022-11-22

- added package level docs
- added citation file
- some minor bug fixes

## edgebundle 0.4.0

CRAN release: 2022-07-04

- changed line straightening in
  [`tnss_tree()`](https://schochastics.github.io/edgebundle/reference/tnss_tree.md)
  to Visvalingam algorithm (NOTE: The meaning of the epsilon parameter
  is now reversed!)
- fixed a bug in
  [`tnss_tree()`](https://schochastics.github.io/edgebundle/reference/tnss_tree.md)
  which created duplicated mesh points

## edgebundle 0.3.1

CRAN release: 2022-01-23

- minor bug fixes

## edgebundle 0.3.0

CRAN release: 2021-10-30

- added
  [`edge_bundle_path()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_path.md)

## edgebundle 0.2.1

CRAN release: 2021-08-09

- Fixed an error that occurred with valgrind

## edgebundle 0.2.0

CRAN release: 2021-08-02

- added more documentation
- cleanup of code
- fixed readme

## edgebundle 0.1.0

- added smoothing for
  [`tnss_tree()`](https://schochastics.github.io/edgebundle/reference/tnss_tree.md)

## edgebundle 0.0.1

- Added a `NEWS.md` file to track changes to the package.
