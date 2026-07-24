# Choosing a bundling algorithm

`edgebundle` implements several edge bundling algorithms. They all take
the same input — a graph plus a node layout — and return the same thing:
a data frame of points along each edge (`x`, `y`, `index`, `group`) that
you plot yourself with, for example, `ggplot2::geom_path()`. Which one
should you use? This article summarises the trade-offs.

``` r

library(edgebundle)
library(igraph)
```

## At a glance

| Method | Function | Directed? | Cost | Good for |
|----|----|----|----|----|
| Force-directed (FDEB) | [`edge_bundle_force()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_force.md) | no | O(E²) | general graphs; the classic look |
| Divided | `edge_bundle_force(directed = TRUE)` | **yes** | O(E²) | directed flow; opposite directions kept in separate lanes |
| Stub | [`edge_bundle_stub()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_stub.md) | no | fast | geographic networks; keeps edges readable near nodes |
| Edge-path | [`edge_bundle_path()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_path.md) | no | O(E · sp) | less ambiguous bundling; routes edges along existing short paths |
| Hammer (KDE) | [`edge_bundle_hammer()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_hammer.md) | no | ~O(E) | large, dense graphs; smooth organic bundles |
| MINGLE | [`edge_bundle_mingle()`](https://schochastics.github.io/edgebundle/reference/edge_bundle_mingle.md) | no | O(E log E) | very large graphs; scales best |

All of them are also reachable through a single dispatcher:

``` r

edge_bundle(g, xy, type = "force")   # or "divided", "stub", "path", "hammer", "mingle"
```

## When to use which

- **Force-directed (`"force"`)** is the classic Holten & van Wijk
  algorithm and a sensible default for small to medium graphs. It is
  O(E²) in the number of edges, so it becomes slow for large graphs.
- **Divided (`"divided"`)** is the directed variant (Selassie et
  al. 2011). Edges running in opposite directions are pulled into
  *separate lanes* instead of merging, so the drawing shows the
  asymmetry of flow. Use it when edge direction matters and you pass a
  directed graph.
- **Stub (`"stub"`)** only bundles the parts of edges near their
  endpoints, which keeps individual edges distinguishable — a good fit
  for geographic networks. Plot the result with
  `ggforce::geom_bezier()`.
- **Edge-path (`"path"`)** routes each edge along a short path in the
  graph if doing so is not too much of a detour (`max_distortion`). This
  produces *less ambiguous* bundles than force-directed bundling.
- **Hammer (`"hammer"`)** bundles by kernel density estimation. It is
  now a native C++ implementation (no Python needed) and scales close to
  linearly, making it a good choice for large, dense graphs.
- **MINGLE (`"mingle"`)** groups edges bottom-up whenever routing them
  together reduces total drawn length. It scales best (O(E log E)) and
  is the right choice for very large graphs.

## Flow maps

For one-to-many *flow* data (a weighted star graph), use a flow map
rather than a bundler.
[`flow_tree()`](https://schochastics.github.io/edgebundle/reference/flow_tree.md)
is the recommended layout: an angle-restricted spiral tree that is
planar, keeps node positions fixed, and needs a single parameter
`alpha`.

``` r

xy <- cbind(state.center$x, state.center$y)[!state.name %in% c("Alaska", "Hawaii"), ]
flow <- flow_tree(cali2010, xy, root = 4, alpha = 40)
```

The older
[`tnss_tree()`](https://schochastics.github.io/edgebundle/reference/tnss_tree.md)
(triangulation + approximate Steiner tree) is kept as an alternative but
requires the `interp` package and a set of dummy nodes.

## A note on ggraph

Since version 2.2.0, `ggraph` supports force-directed edge bundling
natively via `geom_edge_bundle_*()`. `edgebundle` stays useful as a
ggplot-agnostic toolkit (it returns plain data frames) and provides
methods ggraph does not: the divided/directed variant, stub, edge-path,
MINGLE, a Python-free hammer, and the flow map layouts.

## Disclaimer

Edge bundling produces neat-looking visualizations but does not
necessarily improve readability, and the algorithms are sensitive to
their parameters. Consult the original papers and experiment — do not
expect miracles.
