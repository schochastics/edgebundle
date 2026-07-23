# Divided edge bundling (Selassie et al. 2011); reached via
# edge_bundle_force(directed = TRUE). Coordinates are scaled into a 1000x1000
# box for the run (the force constants are tuned for that scale) and unscaled
# for the returned paths.
.edge_bundle_divided <- function(object, xy, lane_width = 25, k_spring = 5e-4,
                                 k_charge = 2e4, decay = 35, friction = 0.8,
                                 step = 20, passes = 5, iterations = 30,
                                 use_connectivity = TRUE) {
    g <- .as_igraph(object)
    if (is.null(g)) {
        stop("directed edge bundling requires an `igraph` or `tbl_graph` object")
    }
    inp <- .bundle_inputs(g, xy)
    el <- inp$el
    xy <- inp$xy

    rngx <- range(xy[, 1])
    rngy <- range(xy[, 2])
    fac <- 1000 / max(diff(rngx), diff(rngy))
    sxy <- cbind((xy[, 1] - rngx[1]) * fac, (xy[, 2] - rngy[1]) * fac)
    exy <- cbind(sxy[el[, 1], 1], sxy[el[, 1], 2], sxy[el[, 2], 1], sxy[el[, 2], 2])

    if (use_connectivity) {
        node_dist <- igraph::distances(igraph::as_undirected(g, mode = "collapse"))
    } else {
        node_dist <- matrix(0, 1, 1)
    }

    w <- igraph::edge_attr(g, "weight")
    if (is.null(w)) {
        w <- rep(1, nrow(el))
    } else {
        w <- as.numeric(w) / max(as.numeric(w))
    }

    elist <- divided_bundle_iter(
        exy, el, node_dist, w, k_spring, k_charge, lane_width, decay,
        friction, step, passes, iterations, use_connectivity
    )

    m <- nrow(el)
    coords <- do.call(rbind, elist)
    coords[, 1] <- coords[, 1] / fac + rngx[1]
    coords[, 2] <- coords[, 2] / fac + rngy[1]

    .as_bundle_df(coords, m, nrow(elist[[1]]))
}
