test_that("convert_edges.igraph returns an m x 4 endpoint matrix", {
    g <- igraph::graph_from_edgelist(matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE), FALSE)
    xy <- cbind(c(0, 1, 2), c(0, 1, 0))
    ce <- convert_edges(g, xy)
    expect_true(is.matrix(ce))
    expect_equal(dim(ce), c(2L, 4L))
    expect_equal(ce[1, ], c(0, 0, 1, 1))
})

test_that("convert_edges errors when coords do not match vertex count", {
    g <- igraph::graph_from_edgelist(matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE), FALSE)
    expect_error(convert_edges(g, cbind(0:1, 0:1)), "number of rows")
})

test_that("convert_edges default method rejects unknown classes", {
    expect_error(convert_edges(list(), cbind(0, 0)), "don't know how to handle")
})
