# edgebundle (development version)

* switched the test suite to testthat edition 3 and added correctness tests for every exported function
* modernised all `igraph` calls to the 2.x API, removing deprecation warnings
* `edge_bundle_force()` now also accepts a two-column edgelist matrix
* fixed `edge_bundle_stub()`: the angular grouping was broken (wrong circular-gap handling and a bundle-size cap that summed cluster ids instead of counting edges) and vertical edges could produce `NaN`. Bundling output changes as a result
* fixed `edge_bundle_path()`: the path-length used in the distortion test was computed along the wrong vertices, so edges were essentially never routed. Bundling output changes as a result
* added divided edge bundling for directed graphs via `edge_bundle_force(directed = TRUE)` (Selassie et al. 2011): edges running in opposite directions are kept in separate lanes

# edgebundle 0.4.2

* clean up codebase

# edgebundle 0.4.1

* added package level docs
* added citation file
* some minor bug fixes

# edgebundle 0.4.0

* changed line straightening in `tnss_tree()` to Visvalingam algorithm (NOTE: The meaning of the epsilon parameter is now reversed!)
* fixed a bug in `tnss_tree()` which created duplicated mesh points

# edgebundle 0.3.1

* minor bug fixes

# edgebundle 0.3.0

* added `edge_bundle_path()`

# edgebundle 0.2.1

* Fixed an error that occurred with valgrind

# edgebundle 0.2.0

* added more documentation
* cleanup of code
* fixed readme

# edgebundle 0.1.0

* added smoothing for `tnss_tree()`

# edgebundle 0.0.1

* Added a `NEWS.md` file to track changes to the package.
