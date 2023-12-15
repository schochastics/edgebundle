#' @title Edge-Path Bundling
#' @description Implements edge-path bundling.
#' @details This is a re-implementation of https://github.com/mwallinger-tu/edge-path-bundling
#' @param g an igraph object
#' @param xy coordinates of vertices
#' @param max_distortion maximum distortion
#' @param weight_fac edge weight factor
#' @param segments number of subdivisions of edges
#' @return data.frame containing the bundled edges
#' @author David Schoch
#' @details see [online](https://github.com/schochastics/edgebundle) for plotting tips
#' @seealso [edge_bundle_hammer],[edge_bundle_stub],[edge_bundle_force]
#' @references
#' Wallinger, M., Archambault, D., Auber, D., Nollenburg, M., & Peltonen, J. (2021). Edge-Path Bundling: A Less Ambiguous Edge Bundling Approach. IEEE Transactions on Visualization and Computer Graphics.
#' @examples
#' library(igraph)
#' g <- graph_from_edgelist(matrix(c(
#'     1, 2, 1, 6,
#'     1, 4, 2, 3, 3, 4, 4, 5, 5, 6
#' ), ncol = 2, byrow = TRUE), FALSE)
#' xy <- cbind(c(0, 10, 25, 40, 50, 50), c(0, 15, 25, 15, 0, -10))
#' edge_bundle_path(g, xy)
#' @export

edge_bundle_path <- function(g, xy, max_distortion = 2, weight_fac = 2, segments = 20) {
    # preprocess
    if (!igraph::is.igraph(g)) {
        stop("edge_bundle_path requires the input graph to be an ingraph object")
    }
    m <- igraph::ecount(g)
    lock <- rep(FALSE, m)
    skip <- rep(FALSE, m)

    el <- igraph::get.edgelist(g, names = FALSE)
    exy <- cbind(
        xy[el[, 1], 1], xy[el[, 1], 2],
        xy[el[, 2], 1], xy[el[, 2], 2]
    )
    elen <- sqrt((exy[, 1] - exy[, 3])^2 + (exy[, 2] - exy[, 4])^2)
    weights <- elen^weight_fac
    sortedEdges <- order(weights, decreasing = TRUE)
    cpoints <- vector("list", m)
    # iterate
    for (e in sortedEdges) {
        s <- el[e, 1]
        t <- el[e, 2]
        cpoints[[e]] <- xy[c(s, t), ]
        if (lock[e]) {
            next()
        }
        skip[e] <- TRUE
        g1 <- igraph::delete.edges(g, which(skip))
        sp_verts <- suppressWarnings(igraph::shortest_paths(g1, s, t, weights = weights[!skip])$vpath[[1]])
        if (length(sp_verts) < 2) {
            skip[e] <- FALSE
            next
        }
        sp_len <- path_length(sp_verts, xy)
        if (sp_len >= max_distortion * elen[e]) {
            skip[e] <- FALSE
            next
        }
        lock[igraph::get.edge.ids(g, rep(as.integer(sp_verts), each = 2)[-c(1, 2 * length(sp_verts))])] <- TRUE
        cpoints[[e]] <- xy[sp_verts, ]
    }
    cpoints_bezier <- lapply(cpoints, approximateBezier, n = segments)

    idx <- seq(0, 1, length.out = segments)
    data_bundle <- as.data.frame(cbind(
        do.call("rbind", cpoints_bezier),
        rep(idx, m),
        rep(1:m, each = segments)
    ))

    names(data_bundle) <- c("x", "y", "index", "group")

    data_bundle
}

path_length <- function(verts, xy) {
    plen <- 0
    for (i in 1:(length(verts) - 1)) {
        plen <- plen + sqrt((xy[i, 1] - xy[i + 1, 1])^2 + (xy[i, 2] - xy[i + 1, 2])^2)
    }
    plen
}

approximateBezier <- function(points, n) {
    pnrow <- nrow(points) - 1
    tseq <- seq(0, 1, length.out = n)
    if (pnrow == 1) {
        bezier <- cbind(
            tseq * points[1, 1] + (1 - tseq) * points[2, 1],
            tseq * points[1, 2] + (1 - tseq) * points[2, 2]
        )
    }
    binoms <- choose(pnrow, seq(0, pnrow))
    bezier <- matrix(0, length(tseq), 2)
    b <- 1
    for (t in tseq) {
        p <- c(0, 0)
        for (i in 0:pnrow) {
            tpi <- (1 - t)^(pnrow - i)
            coeff <- tpi * t^i
            p[1] <- p[1] + binoms[i + 1] * coeff * points[i + 1, 1]
            p[2] <- p[2] + binoms[i + 1] * coeff * points[i + 1, 2]
        }
        bezier[b, ] <- p
        b <- b + 1
    }
    bezier
}
