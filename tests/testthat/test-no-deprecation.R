# Guards the igraph API modernisation: exported functions must not emit
# deprecation (or any) warnings.

test_that("edge bundlers run without warnings", {
    fx <- force_fixture()
    expect_no_warning(edge_bundle_force(fx$g, fx$xy))
    sx <- stub_fixture()
    expect_no_warning(edge_bundle_stub(sx$g, sx$xy))
    px <- path_fixture()
    expect_no_warning(edge_bundle_path(px$g, px$xy))
})

test_that("flow map and metro run without warnings", {
    fx <- tnss_fixture()
    d <- tnss_dummies(fx$xy, root = fx$root)
    g <- expect_no_warning(tnss_tree(cali2010, fx$xy, d, root = fx$root, order = "near"))
    expect_no_warning(tnss_smooth(g, bw = 10, n = 10))

    mg <- igraph::simplify(metro_berlin)
    xy <- cbind(igraph::V(mg)$lon, igraph::V(mg)$lat) * 100
    expect_no_warning(metro_multicriteria(mg, xy, l = 2, gr = 0.5, w = c(100, 100, 1, 1, 100), bsize = 5))
})
