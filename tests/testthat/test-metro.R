test_that("metro_multicriteria returns new station coordinates", {
    g <- igraph::simplify(metro_berlin)
    xy <- cbind(igraph::V(g)$lon, igraph::V(g)$lat) * 100
    xy_new <- metro_multicriteria(g, xy, l = 2, gr = 0.5, w = c(100, 100, 1, 1, 100), bsize = 5)
    expect_true(is.matrix(xy_new))
    expect_equal(dim(xy_new), dim(xy))
    expect_false(anyNA(xy_new))
})
