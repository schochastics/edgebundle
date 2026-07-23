# Two parallel edges between the left and right sides, offset in y by 1.
divided_fixture <- function(second = c(3, 4)) {
    xy <- rbind(c(0, 0), c(10, 0), c(0, 1), c(10, 1))
    g <- igraph::graph_from_edgelist(
        matrix(c(1, 2, second), ncol = 2, byrow = TRUE), directed = TRUE
    )
    list(g = g, xy = xy)
}

mid_gap <- function(res) {
    e1 <- res[res$group == 1, ]
    e2 <- res[res$group == 2, ]
    abs(e1$y[round(nrow(e1) / 2)] - e2$y[round(nrow(e2) / 2)])
}

test_that("directed bundling returns the canonical data frame", {
    fx <- divided_fixture()
    res <- edge_bundle_force(fx$g, fx$xy, directed = TRUE)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("x", "y", "index", "group"))
    expect_equal(sort(unique(res$group)), 1:2)
    expect_equal(unname(unique(table(res$group))), 33L)
    expect_true(all(res$index >= 0 & res$index <= 1))
    expect_true(all(is.finite(res$x)) && all(is.finite(res$y)))
})

test_that("directed bundling keeps edge endpoints fixed", {
    fx <- divided_fixture()
    res <- edge_bundle_force(fx$g, fx$xy, directed = TRUE)
    e1 <- res[res$group == 1, ]
    expect_equal(c(e1$x[1], e1$y[1]), c(0, 0))
    expect_equal(c(e1$x[nrow(e1)], e1$y[nrow(e1)]), c(10, 0))
})

test_that("same-direction edges merge, opposite-direction edges stay in lanes", {
    same <- edge_bundle_force(divided_fixture(c(3, 4))$g, divided_fixture()$xy,
        directed = TRUE, use_connectivity = FALSE
    )
    opp <- edge_bundle_force(divided_fixture(c(4, 3))$g, divided_fixture()$xy,
        directed = TRUE, use_connectivity = FALSE
    )
    expect_lt(mid_gap(same), 0.05)
    expect_gt(mid_gap(opp), 0.1)
})

test_that("directed bundling needs a graph object", {
    expect_error(
        edge_bundle_force(list(), cbind(0, 0), directed = TRUE),
        "requires an .igraph"
    )
})
