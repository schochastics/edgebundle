#' @title stub edge bundling
#' @description Implements the stub edge bundling by Nocaj and Brandes
#' @param object a graph object (igraph/tbl_graph). Does not support network objects
#' @param xy coordinates of vertices
#' @param alpha maximal angle (in degree) between consecutive edges in a bundle
#' @param beta angle (in degree) at which to connect two stubs
#' @param gamma maximal overall angle (in degree) of an edge bundle
#' @param tshift numeric between 0 and 1. The closer to one, the longer the bigger bundle
#' @param t numeric between 0 and 1. control point location
#' @return data.frame containing the bundled edges
#' @author David Schoch
#' @details see [online](https://github.com/schochastics/edgebundle) for plotting tips
#' @seealso [edge_bundle_hammer],[edge_bundle_force], [edge_bundle_path]
#' @references
#' Nocaj, Arlind, and Ulrik Brandes. "Stub bundling and confluent spirals for geographic networks." International Symposium on Graph Drawing. Springer, Cham, 2013.
#' @examples
#' library(igraph)
#' g <- make_star(10, "undirected")
#'
#' xy <- matrix(c(
#'     0, 0,
#'     cos(90 * pi / 180), sin(90 * pi / 180),
#'     cos(80 * pi / 180), sin(80 * pi / 180),
#'     cos(70 * pi / 180), sin(70 * pi / 180),
#'     cos(330 * pi / 180), sin(330 * pi / 180),
#'     cos(320 * pi / 180), sin(320 * pi / 180),
#'     cos(310 * pi / 180), sin(310 * pi / 180),
#'     cos(210 * pi / 180), sin(210 * pi / 180),
#'     cos(200 * pi / 180), sin(200 * pi / 180),
#'     cos(190 * pi / 180), sin(190 * pi / 180)
#' ), ncol = 2, byrow = TRUE)
#'
#' edge_bundle_stub(g, xy)
#' # use ggforce::geom_bezier for plotting
#' @export

edge_bundle_stub <- function(object, xy, alpha = 11, beta = 75, gamma = 40, t = 0.5, tshift = 0.5) {
    if (any(class(object) == "igraph")) {
        if (!requireNamespace("igraph", quietly = TRUE)) {
            stop("The `igraph` package is required for this functionality")
        }
        el <- igraph::as_edgelist(object, FALSE)
        adj <- igraph::as_adj_list(object, "all")
    } else if (any(class(object) == "tbl_graph")) {
        if (!requireNamespace("tidygraph", quietly = TRUE)) {
            stop("The `tidygraph` package is required for this functionality")
        }
        object <- tidygraph::as.igraph(object)
        el <- igraph::as_edgelist(object, FALSE)
        adj <- igraph::as_adj_list(object, "all")
    } else if (any(class(object) == "network")) {
        stop("`network` objects not supported. Convert object to `igraph` or `tbl_graph` first.")
    } else {
        stop("only `igraph` or `tbl_graph` objects supported.")
    }

    idx <- seq(0, 1, length.out = 8) # 4 control points per stub/2 stubs

    B <- compute_bundle_list(xy, adj, gamma, alpha)

    res <- vector("list", nrow(el))
    for (e in seq_len(nrow(el))) {
        v <- el[e, 1]
        w <- el[e, 2]
        idw <- which(adj[[w]] == v)
        Bw <- adj[[w]][which(B[[w]] == B[[w]][idw])]
        idv <- which(adj[[v]] == w)
        Bv <- adj[[v]][which(B[[v]] == B[[v]][idv])]
        pv <- xy[v, ]
        pw <- xy[w, ]
        tst <- as.data.frame(control_points(xy, pv, pw, Bv, Bw, tshift, beta, t))
        tst$idx <- idx
        tst$grp <- c(rep(paste0(e, ".", 1), 4), rep(paste0(e, ".", 2), 4))
        rownames(tst) <- NULL
        res[[e]] <- tst
    }
    data_bundle <- do.call(rbind, res)
    names(data_bundle) <- c("x", "y", "index", "group")

    data_bundle
}


euclidean_dist <- function(p, q) {
    sqrt(sum((p - q)^2))
}

# Acute angle (radians) between the two edges, robust to vertical edges.
angle_edge <- function(P, Q) {
    ap <- atan2(P[4] - P[2], P[3] - P[1])
    aq <- atan2(Q[4] - Q[2], Q[3] - Q[1])
    d <- abs(ap - aq) %% pi
    if (d > pi / 2) d <- pi - d
    d
}

angle_edge_vec <- Vectorize(angle_edge)

bundle_edges <- function(edges_xy, gamma, alpha) {
    n <- nrow(edges_xy)
    if (n == 1) {
        return(1)
    }
    vec_x <- edges_xy[, 3] - edges_xy[, 1]
    vec_y <- edges_xy[, 4] - edges_xy[, 2]
    ang <- atan2(vec_y, vec_x) * 180 / pi
    ang <- ifelse(ang < 0, 360 + ang, ang)
    aord <- order(ang)
    av <- ang[aord]

    # forward angular gap to the next edge around the circle; gaps sum to 360
    gap <- c(av[-1], av[1] + 360) - av

    # walk the edges as a ring starting just after the widest gap, so the arcs
    # between wide gaps stay contiguous; split when a gap exceeds alpha or the
    # bundle reaches gamma edges
    start <- which.max(gap)
    walk <- ((start + seq_len(n) - 1) %% n) + 1
    bundles <- integer(n)
    cur <- 1L
    count <- 0L
    prev <- NA_integer_
    for (p in walk) {
        if (!is.na(prev) && (gap[prev] >= alpha || count >= gamma)) {
            cur <- cur + 1L
            count <- 0L
        }
        bundles[aord[p]] <- cur
        count <- count + 1L
        prev <- p
    }
    bundles
}

compute_bundle_list <- function(xy, adj, gamma, alpha) {
    bundles_lst <- lapply(seq_len(length(adj)), function(i) {
        if (length(adj[[i]]) > 1) {
            edges_xy <- (cbind(matrix(xy[i, ], nrow = length(adj[[i]]), ncol = 2, byrow = TRUE), xy[adj[[i]], ]))
        } else {
            edges_xy <- t((c(xy[i, ], xy[adj[[i]], ])))
        }
        bundle_edges(edges_xy, gamma, alpha)
    })
    bundles_lst
}

weighted_midpoint <- function(pv, pw, Bv, Bw, tshift) {
    0.5 * (pv + pw) + (Bv / (Bv + Bw) - 0.5) * tshift * (pw - pv)
}

centroid <- function(Bv, xy) {
    if (length(Bv) > 1) {
        colMeans(xy[Bv, ])
    } else {
        xy[Bv, ]
    }
}

control_points <- function(xy, pv, pw, Bv, Bw, tshift, beta, t) {
    # pv2 is the point on c such that angle(pv,pv2,pm)=beta
    beta <- beta * pi / 180

    # weighted midpoints
    pm <- weighted_midpoint(pv, pw, length(Bv), length(Bw), tshift)
    cv <- centroid(Bv, xy)
    cw <- centroid(Bw, xy)

    # line between pv and cv
    slope_v <- (cv[2] - pv[2]) / (cv[1] - pv[1])
    incep_v <- cv[2] - slope_v * cv[1]

    # line between pw and cw
    slope_w <- (cw[2] - pw[2]) / (cw[1] - pw[1])
    incep_w <- cw[2] - slope_w * cw[1]

    # control points v -----------------------------------------------------------
    # angle pm-pv-cv
    alpha <- angle_edge(c(pv, pm), c(pv, cv)) * pi / 180
    # remaining angle in triangle
    gamma <- pi - alpha - beta

    # length of side opposite of beta
    b <- euclidean_dist(pv, pm)

    l_pvcv <- sin(gamma) * b / sin(beta)

    pv2x <- c(pv[1] + l_pvcv / (sqrt(1 + slope_v^2)), pv[1] - l_pvcv / (sqrt(1 + slope_v^2)))
    pv2y <- slope_v * pv2x + incep_v

    idx <- which.min(c(euclidean_dist(pm, c(pv2x[1], pv2y[1])), euclidean_dist(pm, c(pv2x[2], pv2y[2]))))
    pv2 <- c(pv2x[idx], pv2y[idx])

    pv1 <- pv + t * (pv2 - pv)

    # control points w -----------------------------------------------------------
    # angle pm-pv-cv
    alpha <- angle_edge(c(pw, pm), c(pw, cw)) * pi / 180
    # remaining angle in triangle
    gamma <- pi - alpha - beta

    # length of side opposite of beta
    b <- euclidean_dist(pw, pm)

    l_pwcw <- sin(gamma) * b / sin(beta)

    pw2x <- c(pw[1] + l_pwcw / (sqrt(1 + slope_w^2)), pw[1] - l_pwcw / (sqrt(1 + slope_w^2)))
    pw2y <- slope_w * pw2x + incep_w

    idx <- which.min(c(euclidean_dist(pm, c(pw2x[1], pw2y[1])), euclidean_dist(pm, c(pw2x[2], pw2y[2]))))
    pw2 <- c(pw2x[idx], pw2y[idx])

    pw1 <- pw + t * (pw2 - pw)

    # control point both
    pvwm <- 0.5 * (pv2 + pw2)

    rbind(
        rbind(pv, pv1, pv2, pvwm),
        # rbind(pvwm,pw2,pw1,pw)
        rbind(pw, pw1, pw2, pvwm)
    )
}
