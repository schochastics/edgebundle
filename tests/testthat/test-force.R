test_that("edge_bundle_force returns the canonical bundle data frame", {
    fx <- force_fixture()
    res <- edge_bundle_force(fx$g, fx$xy)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("x", "y", "index", "group"))
    expect_equal(sort(unique(res$group)), 1:6)
    expect_equal(unname(unique(table(res$group))), 34L)
    expect_true(all(res$index >= 0 & res$index <= 1))
})

test_that("edge_bundle_force is deterministic", {
    fx <- force_fixture()
    expect_equal(edge_bundle_force(fx$g, fx$xy), edge_bundle_force(fx$g, fx$xy))
})

# Characterisation pin: the Phase 1 refactor must not change force output.
test_that("edge_bundle_force output is stable (invariance pin for the refactor)", {
    fx <- force_fixture()
    res <- edge_bundle_force(fx$g, fx$xy)
    expect_equal(dim(res), c(204L, 4L))
    expect_equal(res$x[1], 0)
    expect_equal(res$y[1], 1)
    expect_equal(sum(res$x), 101.261, tolerance = 1e-3)
    expect_equal(sum(res$y), 714, tolerance = 1e-6)
})
