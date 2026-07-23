# hammer edge bundling

Implements the hammer edge bundling by Ian Calvert.

## Usage

``` r
edge_bundle_hammer(object, xy, bw = 0.05, decay = 0.7)
```

## Arguments

- object:

  a graph object (igraph/network/tbl_graph)

- xy:

  coordinates of vertices

- bw:

  bandwidth parameter

- decay:

  decay parameter

## Value

data.frame containing the bundled edges

## Details

This function only wraps existing python code from the datashader
library. Original code can be found at
https://gitlab.com/ianjcalvert/edgehammer. Datashader is a huge library
with a lot of dependencies, so think twice if you want to install it
just for edge bundling. Check
https://datashader.org/user_guide/Networks.html for help concerning
parameters bw and decay. To install all dependencies, use
[install_bundle_py](https://schochastics.github.io/edgebundle/reference/install_bundle_py.md).

see [online](https://github.com/schochastics/edgebundle) for plotting
tips

## See also

[edge_bundle_force](https://schochastics.github.io/edgebundle/reference/edge_bundle_force.md),[edge_bundle_stub](https://schochastics.github.io/edgebundle/reference/edge_bundle_stub.md),
[edge_bundle_path](https://schochastics.github.io/edgebundle/reference/edge_bundle_path.md)

## Author

David Schoch
