# edgebundle (development version)

* switched the test suite to testthat edition 3 and added correctness tests for every exported function
* modernised all `igraph` calls to the 2.x API, removing deprecation warnings
* `edge_bundle_force()` now also accepts a two-column edgelist matrix

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
