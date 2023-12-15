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
#' g <- graph.star(10, "undirected")
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
        el <- igraph::get.edgelist(object, FALSE)
        adj <- igraph::get.adjlist(object, "all")
    } else if (any(class(object) == "tbl_graph")) {
        if (!requireNamespace("tidygraph", quietly = TRUE)) {
            stop("The `tidygraph` package is required for this functionality")
        }
        object <- tidygraph::as.igraph(object)
        el <- igraph::get.edgelist(object, FALSE)
        adj <- igraph::get.adjlist(object, "all")
    } else if (any(class(object) == "network")) {
        el <- network::as.edgelist(object)

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

angle_edge <- function(P, Q) {
    pvec <- c(P[3] - P[1], P[4] - P[2])
    qvec <- c(Q[3] - Q[1], Q[4] - Q[2])
    # dot_pq <- sum(pvec*qvec)
    # mag_pq <- sum(pvec^2)*sum(qvec^2)
    # acos(dot_pq/mag_pq)*180/pi

    m1 <- pvec[2] / pvec[1]
    m2 <- qvec[2] / qvec[1]
    a <- atan((m1 - m2) / (1 + m1 * m2))
    aa <- c(a * 180 / pi, -a * 180 / pi)
    aa[aa < 0] <- 360 + aa[aa < 0]
    aa <- aa * pi / 180
    min(aa)
}

angle_edge_vec <- Vectorize(angle_edge)

bundle_edges <- function(edges_xy, gamma, alpha) {
    if (nrow(edges_xy) == 1) {
        return(1)
        # return(t(c(1,1)))
    }
    dat <- as.data.frame(edges_xy)
    dat[["V3"]] <- dat[["V3"]] - dat[["V1"]]
    dat[["V4"]] <- dat[["V4"]] - dat[["V2"]]
    dat[["V1"]] <- dat[["V2"]] <- 0
    elen <- sqrt(dat[["V3"]]^2 + dat[["V4"]]^2)
    dat[["x"]] <- dat[["V3"]] / elen
    dat[["y"]] <- dat[["V4"]] / elen
    dat[["angle"]] <- atan2(dat[["y"]], dat[["x"]]) * 180 / pi
    dat[["angle"]] <- ifelse(dat[["angle"]] < 0, 360 + dat[["angle"]], dat[["angle"]])
    aord <- order(dat[["angle"]])
    angle_vec <- dat[["angle"]][aord]
    w <- pmin(
        abs(angle_vec - angle_vec[c(2:length(angle_vec), 1)]),
        360 + angle_vec[c(2:length(angle_vec), 1)] - angle_vec
    )

    if (length(w) == 1) {
        if (w > alpha) {
            return(c(1, 2))
        } else {
            return(c(1, 1))
        }
    }

    start <- which.min(w)
    bundles <- rep(0, length(aord))
    if (start != 1) {
        w <- w[c(start:length(w), 1:(start - 1))]
        aord <- aord[c(start:length(w), 1:(start - 1))]
    }
    bundles[aord[1]] <- 1
    bundles[aord[2]] <- 1
    cur <- 1
    for (i in 2:length(w) - 1) {
        if (w[i] < alpha && sum(bundles[bundles == cur]) < gamma) {
            bundles[aord[i + 1]] <- cur
        } else {
            cur <- cur + 1
            bundles[aord[i + 1]] <- cur
        }
    }

    # pa <- graph.ring(nrow(edges_xy),directed = FALSE)
    # E(pa)$weight <- 360-w
    # bundles <- cluster_louvain(pa,weights = E(pa)$weight)$membership[order(aord)]

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
