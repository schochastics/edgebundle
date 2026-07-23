# Shared fixtures for the test suite.

# Bipartite crossing graph used in the edge_bundle_force() example: deterministic.
force_fixture <- function() {
    g <- igraph::graph_from_edgelist(
        matrix(c(1, 12, 2, 11, 3, 10, 4, 9, 5, 8, 6, 7), ncol = 2, byrow = TRUE),
        FALSE
    )
    xy <- cbind(c(rep(0, 6), rep(1, 6)), c(1:6, 1:6))
    list(g = g, xy = xy)
}

# Star graph used in the edge_bundle_stub() example.
stub_fixture <- function() {
    g <- igraph::make_star(10, "undirected")
    ang <- c(90, 80, 70, 330, 320, 310, 210, 200, 190)
    xy <- rbind(c(0, 0), cbind(cos(ang * pi / 180), sin(ang * pi / 180)))
    list(g = g, xy = xy)
}

# Small graph with one long edge that a detour can bundle: edge_bundle_path() example.
path_fixture <- function() {
    g <- igraph::graph_from_edgelist(
        matrix(c(1, 2, 1, 6, 1, 4, 2, 3, 3, 4, 4, 5, 5, 6), ncol = 2, byrow = TRUE),
        FALSE
    )
    xy <- cbind(c(0, 10, 25, 40, 50, 50), c(0, 15, 25, 15, 0, -10))
    list(g = g, xy = xy)
}

# US states (minus Alaska/Hawaii) with California as root, for the flow map.
tnss_fixture <- function() {
    xy <- cbind(state.center$x, state.center$y)[!state.name %in% c("Alaska", "Hawaii"), ]
    list(xy = xy, root = 4)
}
