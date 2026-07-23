# Create Steiner tree from real and dummy points

creates an approximated Steiner tree for a flow map visualization

## Usage

``` r
tnss_tree(
  g,
  xy,
  xydummy,
  root,
  gamma = 0.9,
  epsilon = 0.3,
  elen = Inf,
  order = "random"
)
```

## Arguments

- g:

  original flow network (must be a one-to-many flow network, i.e star
  graph). Must have a weight attribute indicating the flow

- xy:

  coordinates of "real" nodes

- xydummy:

  coordinates of "dummy" nodes

- root:

  root node id of the flow

- gamma:

  edge length decay parameter

- epsilon:

  percentage of points keept on a line after straightening with
  Visvalingam Algorithm

- elen:

  maximal length of edges in triangulation

- order:

  in which order shortest paths are calculated
  ("random","weight","near","far")

## Value

approximated Steiner tree from dummy and real nodes as igraph object

## Details

Use
[tnss_smooth](https://schochastics.github.io/edgebundle/reference/tnss_smooth.md)
to smooth the edges of the tree

## References

Sun, Shipeng. "An automated spatial flow layout algorithm using
triangulation, approximate Steiner tree, and path smoothing." AutoCarto,
2016.

## Author

David Schoch

## Examples

``` r
xy <- cbind(state.center$x, state.center$y)[!state.name %in% c("Alaska", "Hawaii"), ]
xy_dummy <- tnss_dummies(xy, root = 4)
gtree <- tnss_tree(cali2010, xy, xy_dummy, root = 4, gamma = 0.9)
```
