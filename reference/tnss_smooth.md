# Smooth a Steiner tree

Converts the Steiner tree to smooth paths

## Usage

``` r
tnss_smooth(g, bw = 3, n = 10)
```

## Arguments

- g:

  Steiner tree computed with
  [tnss_tree](https://schochastics.github.io/edgebundle/reference/tnss_tree.md)

- bw:

  bandwidth of Gaussian Kernel

- n:

  number of extra nodes to include per edge

## Value

data.frame containing the smoothed paths

## Details

see see [online](https://github.com/schochastics/edgebundle) for tips on
plotting the result

## Author

David Schoch

## Examples

``` r
xy <- cbind(state.center$x, state.center$y)[!state.name %in% c("Alaska", "Hawaii"), ]
xy_dummy <- tnss_dummies(xy, root = 4)
gtree <- tnss_tree(cali2010, xy, xy_dummy, root = 4, gamma = 0.9)
tree_smooth <- tnss_smooth(gtree, bw = 10, n = 10)
```
