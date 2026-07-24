# Spiral flow tree

Computes a one-to-many flow map as an angle-restricted spiral tree
(Verbeek, Buchin and Speckmann 2011). This is the recommended flow map
layout; unlike
[tnss_tree](https://schochastics.github.io/edgebundle/reference/tnss_tree.md)
it needs no dummy nodes or triangulation, keeps the real node positions
fixed, produces a planar (crossing-free) tree with smooth, shallow-angle
merges, and is tuned by a single parameter `alpha`.

## Usage

``` r
flow_tree(object, xy, root, alpha = 40, n = 20, optimize = FALSE)
```

## Arguments

- object:

  a one-to-many flow network (star graph) as igraph/tbl_graph, with a
  `weight` edge attribute giving the flow

- xy:

  coordinates of vertices

- root:

  root node id of the flow

- alpha:

  restricting angle in degrees (0-90). Smaller values bundle more
  tightly toward the root; typical values are 20-45.

- n:

  number of points sampled per tree edge

- optimize:

  logical. If `TRUE`, refine the tree with an approximate FLOWTREE
  optimization (Verbeek et al. 2011, section 5): join points stay fixed
  while edge interiors are relaxed for smoothness and to keep clear of
  node obstacles. This can slightly relax the strict `alpha` bound.

## Value

data.frame with columns `x`, `y`, `flow` and `edge` (one smooth arc per
tree edge)

## References

Verbeek, Kevin, Kevin Buchin, and Bettina Speckmann. "Flow map layout
via spiral trees." IEEE Transactions on Visualization and Computer
Graphics 17, no. 12 (2011): 2536-2545.

## See also

[tnss_tree](https://schochastics.github.io/edgebundle/reference/tnss_tree.md)

## Author

David Schoch

## Examples

``` r
xy <- cbind(state.center$x, state.center$y)[!state.name %in% c("Alaska", "Hawaii"), ]
flow <- flow_tree(cali2010, xy, root = 4, alpha = 40)
```
