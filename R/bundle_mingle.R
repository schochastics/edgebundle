#' @title MINGLE edge bundling
#' @description Multilevel agglomerative edge bundling (MINGLE).
#' @details Groups edges bottom-up: two bundles are merged whenever routing them
#' through shared meeting points reduces the total drawn length ("ink"). Meeting
#' points are the geometric medians of the source-side and target-side endpoints.
#' A kNN proximity graph over the edges is built once with a kd-tree and then
#' coarsened level by level (O(E log E)).
#' @param object a graph object (igraph/network/tbl_graph)
#' @param xy coordinates of vertices
#' @param k number of nearest neighbours considered as merge candidates per edge
#' @param segments number of points sampled per bundled edge
#' @param bundle_strength strength of bundling between 0 (straight edges) and 1
#'   (route fully through the meeting points)
#' @return data.frame containing the bundled edges
#' @author David Schoch
#' @details see [online](https://github.com/schochastics/edgebundle) for plotting tips
#' @seealso [edge_bundle_force],[edge_bundle_stub],[edge_bundle_path],[edge_bundle_hammer]
#' @references
#' Gansner, Emden R., Yifan Hu, Stephen North, and Carlos Scheidegger. "Multilevel agglomerative edge bundling for visualizing large graphs." 2011 IEEE Pacific Visualization Symposium (2011): 187-194.
#' @examples
#' library(igraph)
#' g <- graph_from_edgelist(
#'     matrix(c(1, 12, 2, 11, 3, 10, 4, 9, 5, 8, 6, 7), ncol = 2, byrow = TRUE), FALSE
#' )
#' xy <- cbind(c(rep(0, 6), rep(1, 6)), c(1:6, 1:6))
#' edge_bundle_mingle(g, xy)
#' @export
edge_bundle_mingle <- function(object, xy, k = 10, segments = 50, bundle_strength = 0.9) {
    edges_xy <- .bundle_inputs(object, xy)$exy
    m <- nrow(edges_xy)

    elist <- mingle_iter(edges_xy, as.integer(k), segments, bundle_strength)

    .as_bundle_df(do.call("rbind", elist), m, segments)
}
