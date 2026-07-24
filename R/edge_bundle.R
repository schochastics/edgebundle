#' @title Edge bundling
#' @description Dispatches to one of the edge bundling algorithms by name. This
#' is a thin wrapper over the individual `edge_bundle_*()` functions and mirrors
#' the way frontends (e.g. vellumplot) select a bundling method.
#' @param object a graph object (igraph/network/tbl_graph)
#' @param xy coordinates of vertices
#' @param type bundling algorithm: `"force"` (force-directed, Holten),
#'   `"divided"` (force-directed for directed graphs, Selassie et al.),
#'   `"stub"` (Nocaj & Brandes), `"path"` (edge-path, Wallinger et al.),
#'   `"hammer"` (KDE, Hurter et al.), or `"mingle"` (Gansner et al.)
#' @param ... arguments passed on to the selected `edge_bundle_*()` function
#' @return data.frame containing the bundled edges (`x`, `y`, `index`, `group`)
#' @author David Schoch
#' @seealso [edge_bundle_force], [edge_bundle_stub], [edge_bundle_path],
#'   [edge_bundle_hammer], [edge_bundle_mingle]
#' @examples
#' library(igraph)
#' g <- graph_from_edgelist(
#'     matrix(c(1, 12, 2, 11, 3, 10, 4, 9, 5, 8, 6, 7), ncol = 2, byrow = TRUE), FALSE
#' )
#' xy <- cbind(c(rep(0, 6), rep(1, 6)), c(1:6, 1:6))
#' edge_bundle(g, xy, type = "force")
#' @export
edge_bundle <- function(object, xy, type = c("force", "divided", "stub", "path", "hammer", "mingle"), ...) {
    type <- match.arg(type)
    switch(type,
        force = edge_bundle_force(object, xy, ...),
        divided = edge_bundle_force(object, xy, directed = TRUE, ...),
        stub = edge_bundle_stub(object, xy, ...),
        path = edge_bundle_path(object, xy, ...),
        hammer = edge_bundle_hammer(object, xy, ...),
        mingle = edge_bundle_mingle(object, xy, ...)
    )
}
