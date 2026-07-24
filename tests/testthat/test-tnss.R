test_that("tnss_dummies returns a two-column coordinate matrix", {
    fx <- tnss_fixture()
    d <- tnss_dummies(fx$xy, root = fx$root)
    expect_true(is.matrix(d))
    expect_equal(ncol(d), 2L)
    expect_gt(nrow(d), 0L)
})

test_that("tnss_tree builds a connected Steiner tree spanning all leaves", {
    skip_if_not_installed("interp")
    fx <- tnss_fixture()
    d <- tnss_dummies(fx$xy, root = fx$root)
    g <- tnss_tree(cali2010, fx$xy, d, root = fx$root, gamma = 0.9, order = "near")
    expect_s3_class(g, "steiner_tree")
    expect_true(igraph::is_connected(g))
    leafs <- which(igraph::V(g)$tnss == "leaf")
    expect_true(all(igraph::distances(g, fx$root, leafs) < Inf))
})

test_that("tnss_tree order = 'weight' processes leaves by flow", {
    skip_if_not_installed("interp")
    fx <- tnss_fixture()
    d <- tnss_dummies(fx$xy, root = fx$root)
    g <- tnss_tree(cali2010, fx$xy, d, root = fx$root, gamma = 0.9, order = "weight")
    expect_s3_class(g, "steiner_tree")
    expect_true(igraph::is_connected(g))
})

test_that("tnss_smooth returns smoothed paths with flow", {
    skip_if_not_installed("interp")
    fx <- tnss_fixture()
    d <- tnss_dummies(fx$xy, root = fx$root)
    g <- tnss_tree(cali2010, fx$xy, d, root = fx$root, gamma = 0.9, order = "near")
    sm <- tnss_smooth(g, bw = 10, n = 10)
    expect_s3_class(sm, "data.frame")
    expect_named(sm, c("x", "y", "flow", "destination"))
    expect_gt(nrow(sm), 0L)
})

test_that("tnss_smooth rejects non-steiner-tree input", {
    expect_error(tnss_smooth(igraph::make_ring(3)), "steiner tree")
})
