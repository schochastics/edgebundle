# Correctness tests for the bugs fixed in Phase 2.

test_that("path_length measures distance along the given vertex path", {
    xy <- matrix(c(0, 0, 5, 0, 100, 100, 200, 200), ncol = 2, byrow = TRUE)
    # path visiting vertices 3 -> 4
    expect_equal(edgebundle:::path_length(c(3, 4), xy), 100 * sqrt(2))
})

test_that("stub bundle_edges groups by angle without spurious splits", {
    deg <- function(a) a * pi / 180
    mk <- function(angles) t(sapply(angles, function(a) c(0, 0, cos(deg(a)), sin(deg(a)))))
    # {10,20} one bundle; gap; {40,50,60,70} one bundle. Small gamma must not split B.
    b <- edgebundle:::bundle_edges(mk(c(10, 20, 40, 50, 60, 70)), gamma = 5, alpha = 11)
    expect_equal(length(unique(b)), 2L)
    expect_equal(b[1], b[2])
    expect_true(all(b[3] == b[c(4, 5, 6)]))
    expect_false(b[1] == b[3])
})
