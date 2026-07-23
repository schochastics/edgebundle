mingle_parallel <- function(n = 6) {
    xy <- rbind(cbind(0, seq(-2.5, 2.5, length.out = n)),
                cbind(10, seq(-2.5, 2.5, length.out = n)))
    el <- cbind(1:n, (n + 1):(2 * n))
    g <- igraph::graph_from_edgelist(el, FALSE)
    list(g = g, xy = xy)
}

mid_spread <- function(res) {
    mids <- tapply(seq_len(nrow(res)), res$group, function(ix) {
        e <- res[ix, ]
        e$y[round(nrow(e) / 2)]
    })
    diff(range(mids))
}

test_that("edge_bundle_mingle returns the canonical data frame", {
    fx <- mingle_parallel()
    res <- edge_bundle_mingle(fx$g, fx$xy)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("x", "y", "index", "group"))
    expect_equal(sort(unique(res$group)), 1:6)
    expect_equal(unname(unique(table(res$group))), 50L)
    expect_true(all(res$index >= 0 & res$index <= 1))
    expect_true(all(is.finite(res$x)) && all(is.finite(res$y)))
})

test_that("mingle keeps edge endpoints fixed", {
    fx <- mingle_parallel()
    res <- edge_bundle_mingle(fx$g, fx$xy)
    e1 <- res[res$group == 1, ]
    expect_equal(c(e1$x[1], e1$y[1]), c(0, -2.5))
    expect_equal(c(e1$x[nrow(e1)], e1$y[nrow(e1)]), c(10, -2.5))
})

test_that("mingle bundles parallel edges, strength 0 leaves them straight", {
    fx <- mingle_parallel()
    expect_equal(mid_spread(edge_bundle_mingle(fx$g, fx$xy, bundle_strength = 0)), 5)
    expect_lt(mid_spread(edge_bundle_mingle(fx$g, fx$xy, bundle_strength = 0.9)), 2)
})

test_that("mingle leaves a single edge straight", {
    g <- igraph::graph_from_edgelist(matrix(c(1, 2), ncol = 2), FALSE)
    res <- edge_bundle_mingle(g, rbind(c(0, 0), c(10, 10)))
    expect_true(all(abs(res$y - res$x) < 1e-6))
})
