# Sample points for triangulated networks

uses various sampling strategies to create dummy nodes for the
[tnss_tree](https://schochastics.github.io/edgebundle/reference/tnss_tree.md)

## Usage

``` r
tnss_dummies(
  xy,
  root,
  circ = TRUE,
  line = TRUE,
  diag = TRUE,
  grid = FALSE,
  rand = FALSE,
  ncirc = 9,
  rcirc = 2,
  nline = 10,
  ndiag = 50,
  ngrid = 50,
  nrand = 50
)
```

## Arguments

- xy:

  coordinates of "real" nodes

- root:

  root node id

- circ:

  logical. create circular dummy nodes around leafs.

- line:

  logical. create dummy nodes on a straight line between root and leafs.

- diag:

  logical. create dummy nodes diagonally through space.

- grid:

  logical. create dummy nodes on a grid.

- rand:

  logical. create random dummy nodes.

- ncirc:

  numeric. number of circular dummy nodes per leaf.

- rcirc:

  numeric. radius of circles around leaf nodes.

- nline:

  numeric. number of straight line nodes per leaf.

- ndiag:

  numeric. number of dummy nodes on diagonals.

- ngrid:

  numeric. number of dummy nodes per dim on grid.

- nrand:

  numeric. number of random nodes to create.

## Value

coordinates of dummy nodes

## Author

David Schoch

## Examples

``` r
# dummy nodes for tree rooted in California
xy <- cbind(state.center$x, state.center$y)
xy_dummy <- tnss_dummies(xy, 4)
```
