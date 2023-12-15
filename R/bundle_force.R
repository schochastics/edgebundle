#' @title force directed edge bundling
#' @description Implements the classic edge bundling by Holten.
#' @details This is a re-implementation of https://github.com/upphiminn/d3.ForceBundle. Force directed edge bundling is slow (O(E^2)).
#' @param object a graph object (igraph/network/tbl_graph)
#' @param xy coordinates of vertices
#' @param K spring constant
#' @param C number of iteration cycles
#' @param P number of initial edge divisions
#' @param S initial step size
#' @param P_rate rate of edge divisions
#' @param I number of initial iterations
#' @param I_rate rate of iteration decrease per cycle
#' @param compatibility_threshold threshold for when edges are considered compatible
#' @param eps accuracy
#' @return data.frame containing the bundled edges
#' @author David Schoch
#' @details see [online](https://github.com/schochastics/edgebundle) for plotting tips
#' @seealso [edge_bundle_hammer],[edge_bundle_stub],[edge_bundle_path]
#' @references
#' Holten, Danny, and Jarke J. Van Wijk. "Force-Directed Edge Bundling for Graph Visualization." Computer Graphics Forum (Blackwell Publishing Ltd) 28, no. 3 (2009): 983-990.
#' @examples
#' library(igraph)
#' g <- graph_from_edgelist(
#'     matrix(c(
#'         1, 12, 2, 11, 3, 10,
#'         4, 9, 5, 8, 6, 7
#'     ), ncol = 2, byrow = TRUE), FALSE
#' )
#' xy <- cbind(c(rep(0, 6), rep(1, 6)), c(1:6, 1:6))
#' edge_bundle_force(g, xy)
#' @export

edge_bundle_force <- function(object, xy, K = 1, C = 6, P = 1, S = 0.04,
                              P_rate = 2, I = 50, I_rate = 2 / 3,
                              compatibility_threshold = 0.6,
                              eps = 1e-8) {
    # initialize matrix with coordinates
    edges_xy <- convert_edges(object, xy)
    m <- nrow(edges_xy)

    # initialize edge subdivision list
    elist <- unname(lapply(
        split(edges_xy, rep(seq_len(nrow(edges_xy)), ncol(edges_xy))),
        function(y) matrix(y, 2, 2, byrow = TRUE)
    ))

    # main force bundling routine
    elist <- force_bundle_iter(
        edges_xy, elist, K, C, P, P_rate,
        S, I, I_rate, compatibility_threshold, eps
    )

    # assemble data frame
    segments <- nrow(elist[[1]])

    idx <- seq(0, 1, length.out = segments)
    data_bundle <- as.data.frame(cbind(
        do.call("rbind", elist),
        rep(idx, m),
        rep(1:m, each = segments)
    ))

    names(data_bundle) <- c("x", "y", "index", "group")

    data_bundle
}

#' @importFrom Rcpp sourceCpp
NULL

#' @useDynLib edgebundle, .registration = TRUE
NULL
