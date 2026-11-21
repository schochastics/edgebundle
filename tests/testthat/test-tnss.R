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

test_that("tnss_tree flow uses only root-leaf edges (#32)", {
    skip_if_not_installed("interp")
    xy <- cbind(c(0, 1, 1, -1, -1, 2, -2), c(0, 1, -1, 1, -1, 0, 0))
    leaf_flow <- function(g) {
        set.seed(1)
        tr <- tnss_tree(g, xy, tnss_dummies(xy, root = 1), root = 1)
        el <- igraph::as_edgelist(tr, names = FALSE)
        vapply(2:7, function(lf) igraph::E(tr)$flow[el[, 1] == lf | el[, 2] == lf][1], numeric(1))
    }
    star <- igraph::make_star(7, mode = "out")
    igraph::E(star)$weight <- c(10, 20, 30, 40, 50, 60)

    expect_equal(leaf_flow(star), c(10, 20, 30, 40, 50, 60))
    expect_no_warning(
        f <- leaf_flow(igraph::add_edges(star, c(2, 3, 4, 5), weight = c(5, 5)))
    )
    expect_equal(f, c(10, 20, 30, 40, 50, 60))
    expect_no_warning(
        f <- leaf_flow(igraph::add_edges(star, c(1, 2), weight = 15))
    )
    expect_equal(f, c(25, 20, 30, 40, 50, 60))
})

test_that("tnss_tree errors clearly without a weight attribute", {
    xy <- cbind(c(0, 1, -1), c(0, 1, -1))
    g <- igraph::make_star(3, mode = "out")
    expect_error(tnss_tree(g, xy, tnss_dummies(xy, root = 1), root = 1), "weight")
})
