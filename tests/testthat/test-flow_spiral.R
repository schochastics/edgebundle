cali_xy <- function() {
    cbind(state.center$x, state.center$y)[!state.name %in% c("Alaska", "Hawaii"), ]
}

# max angle (deg) between any rendered segment and the radial direction to root
max_radial_angle <- function(paths, rootxy) {
    mx <- 0
    for (k in unique(paths$edge)) {
        p <- paths[paths$edge == k, ]
        for (i in seq_len(nrow(p) - 1)) {
            mid <- c((p$x[i] + p$x[i + 1]) / 2, (p$y[i] + p$y[i + 1]) / 2) - rootxy
            seg <- c(p$x[i + 1] - p$x[i], p$y[i + 1] - p$y[i])
            if (sqrt(sum(seg^2)) < 1e-9 || sqrt(sum(mid^2)) < 1e-9) next
            cosang <- abs(sum(mid * seg)) / (sqrt(sum(mid^2)) * sqrt(sum(seg^2)))
            mx <- max(mx, acos(pmin(1, cosang)) * 180 / pi)
        }
    }
    mx
}

# proper crossings between edges that do not share an endpoint (planarity)
count_crossings <- function(paths) {
    segint <- function(p1, p2, p3, p4) {
        d <- function(a, b, c) (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
        d1 <- d(p3, p4, p1); d2 <- d(p3, p4, p2); d3 <- d(p1, p2, p3); d4 <- d(p1, p2, p4)
        ((d1 > 0 & d2 < 0) | (d1 < 0 & d2 > 0)) & ((d3 > 0 & d4 < 0) | (d3 < 0 & d4 > 0))
    }
    ek <- unique(paths$edge)
    polys <- lapply(ek, function(k) as.matrix(paths[paths$edge == k, c("x", "y")]))
    shares_end <- function(A, B) {
        ends <- rbind(A[1, ], A[nrow(A), ], B[1, ], B[nrow(B), ])
        any(as.matrix(dist(ends)) < 1e-6 & upper.tri(matrix(0, 4, 4)))
    }
    cr <- 0
    for (a in seq_along(polys)) {
        for (b in seq_along(polys)) {
            if (a >= b) next
            A <- polys[[a]]; B <- polys[[b]]
            if (shares_end(A, B)) next
            for (i in seq_len(nrow(A) - 1)) {
                for (j in seq_len(nrow(B) - 1)) {
                    if (segint(A[i, ], A[i + 1, ], B[j, ], B[j + 1, ])) cr <- cr + 1
                }
            }
        }
    }
    cr
}

test_that("flow_tree returns smooth arcs with flow", {
    xy <- cali_xy()
    res <- flow_tree(cali2010, xy, root = 4, alpha = 40, n = 20)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("x", "y", "flow", "edge"))
    expect_equal(unname(unique(table(res$edge))), 20L)
    expect_true(all(is.finite(res$x)) && all(is.finite(res$y)))
    expect_true(all(res$flow > 0))
})

test_that("flow_tree respects the restricting angle alpha", {
    xy <- cali_xy()
    root <- xy[4, ]
    for (a in c(20, 40, 60)) {
        res <- flow_tree(cali2010, xy, root = 4, alpha = a, n = 30)
        expect_lte(max_radial_angle(res, root), a + 1)
    }
})

test_that("flow_tree is planar (no crossings between non-adjacent edges)", {
    xy <- cali_xy()
    res <- flow_tree(cali2010, xy, root = 4, alpha = 40, n = 30)
    expect_equal(count_crossings(res), 0)
})

test_that("flow_tree keeps leaf endpoints fixed at their coordinates", {
    xy <- cali_xy()
    res <- flow_tree(cali2010, xy, root = 4, alpha = 40, n = 20)
    # every leaf coordinate must appear as the start of some edge arc
    starts <- do.call(rbind, lapply(unique(res$edge), function(k) {
        res[res$edge == k, c("x", "y")][1, ]
    }))
    leaves <- which(seq_len(nrow(xy)) != 4)
    hit <- vapply(leaves, function(v) {
        any(abs(starts$x - xy[v, 1]) < 1e-6 & abs(starts$y - xy[v, 2]) < 1e-6)
    }, logical(1))
    expect_true(all(hit))
})

test_that("flow_tree requires a graph object", {
    expect_error(flow_tree(list(), cbind(0, 0), root = 1), "requires an .igraph")
})
