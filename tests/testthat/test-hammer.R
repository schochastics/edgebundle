test_that("edge_bundle_hammer returns the canonical data frame with no Python", {
    fx <- force_fixture()
    res <- edge_bundle_hammer(fx$g, fx$xy)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("x", "y", "index", "group"))
    expect_equal(sort(unique(res$group)), 1:6)
    expect_equal(unname(unique(table(res$group))), 50L)
    expect_true(all(res$index >= 0 & res$index <= 1))
    expect_true(all(is.finite(res$x)) && all(is.finite(res$y)))
})

test_that("hammer bundling pulls nearby parallel edges together", {
    xy <- rbind(c(0, -1.5), c(10, -1.5), c(0, -0.5), c(10, -0.5),
                c(0, 0.5), c(10, 0.5), c(0, 1.5), c(10, 1.5))
    g <- igraph::graph_from_edgelist(
        matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2, byrow = TRUE), FALSE
    )
    res <- edge_bundle_hammer(g, xy)
    mids <- tapply(seq_len(nrow(res)), res$group, function(ix) {
        e <- res[ix, ]
        e$y[round(nrow(e) / 2)]
    })
    expect_lt(diff(range(mids)), 3) # initial spread is 3
    expect_true(all(is.finite(res$y)))
})
