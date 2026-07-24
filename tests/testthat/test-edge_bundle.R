test_that("edge_bundle dispatches to the individual bundlers", {
    fx <- force_fixture()
    expect_equal(edge_bundle(fx$g, fx$xy, type = "force"), edge_bundle_force(fx$g, fx$xy))

    px <- path_fixture()
    expect_equal(edge_bundle(px$g, px$xy, type = "path"), edge_bundle_path(px$g, px$xy))

    sx <- stub_fixture()
    expect_equal(edge_bundle(sx$g, sx$xy, type = "stub"), edge_bundle_stub(sx$g, sx$xy))

    expect_equal(edge_bundle(fx$g, fx$xy, type = "mingle"), edge_bundle_mingle(fx$g, fx$xy))
})

test_that("edge_bundle type = 'divided' routes to directed force bundling", {
    g <- igraph::graph_from_edgelist(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE), directed = TRUE)
    xy <- rbind(c(0, 0), c(10, 0), c(0, 1), c(10, 1))
    expect_equal(
        edge_bundle(g, xy, type = "divided", use_connectivity = FALSE),
        edge_bundle_force(g, xy, directed = TRUE, use_connectivity = FALSE)
    )
})

test_that("edge_bundle forwards tuning parameters", {
    px <- path_fixture()
    expect_equal(
        edge_bundle(px$g, px$xy, type = "path", max_distortion = 3),
        edge_bundle_path(px$g, px$xy, max_distortion = 3)
    )
})

test_that("edge_bundle rejects an unknown type", {
    fx <- force_fixture()
    expect_error(edge_bundle(fx$g, fx$xy, type = "nope"))
})
