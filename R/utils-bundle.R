# Coerce a graph-ish object to a plain igraph, or NULL if it needs the
# convert_edges() S3 fallback (e.g. a `network` object).
.as_igraph <- function(object) {
    if (inherits(object, "tbl_graph")) {
        if (!requireNamespace("tidygraph", quietly = TRUE)) {
            stop("The `tidygraph` package is required for this functionality")
        }
        return(tidygraph::as.igraph(object))
    }
    if (igraph::is_igraph(object)) {
        return(object)
    }
    if (is.matrix(object) && ncol(object) == 2) {
        return(igraph::graph_from_edgelist(object, directed = FALSE))
    }
    NULL
}

# Unified input coercion for the bundlers: endpoint coordinates (`exy`, one row
# x0,y0,x1,y1 per edge), the edgelist, an optional adjacency list, and sizes.
# `network` objects route through the exported convert_edges() S3.
.bundle_inputs <- function(object, xy, need_adj = FALSE) {
    xy <- as.matrix(xy)
    g <- .as_igraph(object)
    if (is.null(g)) {
        exy <- convert_edges(object, xy)
        return(list(g = NULL, el = NULL, exy = exy, adj = NULL, n = nrow(xy), xy = xy))
    }
    if (igraph::vcount(g) != nrow(xy)) {
        stop("number of rows in `xy` does not match number of vertices")
    }
    el <- igraph::as_edgelist(g, names = FALSE)
    exy <- cbind(xy[el[, 1], 1], xy[el[, 1], 2], xy[el[, 2], 1], xy[el[, 2], 2])
    adj <- if (need_adj) igraph::as_adj_list(g, mode = "all") else NULL
    list(g = g, el = el, exy = exy, adj = adj, n = nrow(xy), xy = xy)
}

# Assemble the canonical bundle data frame from `m` edges, each sampled into
# `segments` points, with `index` running 0..1 along every edge.
.as_bundle_df <- function(coords, m, segments) {
    idx <- seq(0, 1, length.out = segments)
    out <- as.data.frame(cbind(coords, rep(idx, m), rep(seq_len(m), each = segments)))
    names(out) <- c("x", "y", "index", "group")
    out
}
