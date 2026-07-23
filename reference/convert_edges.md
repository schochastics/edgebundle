# Convert edges

converts edges of an igraph/network/tidygraph object into format useable
for edge bundling

## Usage

``` r
convert_edges(object, coords)

# Default S3 method
convert_edges(object, coords)

# S3 method for class 'igraph'
convert_edges(object, coords)

# S3 method for class 'network'
convert_edges(object, coords)

# S3 method for class 'tbl_graph'
convert_edges(object, coords)
```

## Arguments

- object:

  graph object

- coords:

  coordinates of vertices

## Value

data frame of edges with coordinates

## Author

David Schoch
