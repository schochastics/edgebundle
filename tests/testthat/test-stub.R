test_that("edge_bundle_stub returns the canonical bundle data frame", {
    fx <- stub_fixture()
    res <- edge_bundle_stub(fx$g, fx$xy)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("x", "y", "index", "group"))
    # two stubs per edge, four control points each
    expect_equal(length(unique(res$group)), 2L * igraph::ecount(fx$g))
    expect_equal(unname(unique(table(res$group))), 4L)
    expect_true(all(res$index >= 0 & res$index <= 1))
})

test_that("edge_bundle_stub rejects unsupported input", {
    expect_error(edge_bundle_stub(list(), cbind(0, 0)), "igraph")
})
