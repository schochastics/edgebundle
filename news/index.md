# Changelog

## edgebundle (development version)

- switched the test suite to testthat edition 3 and added correctness
  tests for every exported function
- modernised all `igraph` calls to the 2.x API, removing deprecation
  warnings
- [`edge_bundle_force()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_force.md)
  now also accepts a two-column edgelist matrix

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
