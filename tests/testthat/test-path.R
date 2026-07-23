test_that("edge_bundle_path returns the canonical bundle data frame", {
    fx <- path_fixture()
    res <- edge_bundle_path(fx$g, fx$xy)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("x", "y", "index", "group"))
    expect_equal(sort(unique(res$group)), seq_len(igraph::ecount(fx$g)))
    expect_equal(unname(unique(table(res$group))), 20L)
    expect_true(all(res$index >= 0 & res$index <= 1))
})

test_that("edge_bundle_path requires an igraph object", {
    expect_error(edge_bundle_path(list(), cbind(0, 0)), "input graph")
})
