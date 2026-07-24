#' @title Spiral flow tree
#' @description Computes a one-to-many flow map as an angle-restricted spiral
#' tree (Verbeek, Buchin and Speckmann 2011). This is the recommended flow map
#' layout; unlike [tnss_tree] it needs no dummy nodes or triangulation, keeps the
#' real node positions fixed, produces a planar (crossing-free) tree with smooth,
#' shallow-angle merges, and is tuned by a single parameter `alpha`.
#' @param object a one-to-many flow network (star graph) as igraph/tbl_graph,
#' with a `weight` edge attribute giving the flow
#' @param xy coordinates of vertices
#' @param root root node id of the flow
#' @param alpha restricting angle in degrees (0-90). Smaller values bundle more
#' tightly toward the root; typical values are 20-45.
#' @param n number of points sampled per tree edge
#' @param optimize logical. If `TRUE`, refine the tree with an approximate
#' FLOWTREE optimization (Verbeek et al. 2011, section 5): join points stay
#' fixed while edge interiors are relaxed for smoothness and to keep clear of
#' node obstacles. This can slightly relax the strict `alpha` bound.
#' @return data.frame with columns `x`, `y`, `flow` and `edge` (one smooth arc
#' per tree edge)
#' @author David Schoch
#' @references
#' Verbeek, Kevin, Kevin Buchin, and Bettina Speckmann. "Flow map layout via spiral trees." IEEE Transactions on Visualization and Computer Graphics 17, no. 12 (2011): 2536-2545.
#' @seealso [tnss_tree]
#' @examples
#' xy <- cbind(state.center$x, state.center$y)[!state.name %in% c("Alaska", "Hawaii"), ]
#' flow <- flow_tree(cali2010, xy, root = 4, alpha = 40)
#' @export
flow_tree <- function(object, xy, root, alpha = 40, n = 20, optimize = FALSE) {
    g <- .as_igraph(object)
    if (is.null(g)) {
        stop("flow_tree requires an `igraph` or `tbl_graph` object")
    }
    xy <- as.matrix(xy)
    nv <- igraph::vcount(g)
    if (nv != nrow(xy)) {
        stop("number of rows in `xy` does not match number of vertices")
    }
    el <- igraph::as_edgelist(g, names = FALSE)
    w <- igraph::edge_attr(g, "weight")
    if (is.null(w)) w <- rep(1, nrow(el))

    # flow per non-root vertex = summed weight of its edges to the root side
    wvec <- numeric(nv)
    for (e in seq_len(nrow(el))) {
        leaf <- if (el[e, 1] == root) el[e, 2] else el[e, 1]
        wvec[leaf] <- wvec[leaf] + w[e]
    }
    leaves <- which(wvec > 0 & seq_len(nv) != root)
    if (length(leaves) == 0) {
        stop("no leaves with flow found; is `root` correct?")
    }

    st <- .spiral_build(xy, root, leaves, wvec[leaves], alpha)
    if (optimize) {
        .spiral_optimize(st, n)
    } else {
        .spiral_render(st, n)
    }
}

# Greedy angle-restricted spiral tree. Returns a mutable-env tree description.
.spiral_build <- function(xy, root, leaves, weight, alpha) {
    tan_a <- tan(alpha * pi / 180)
    rootxy <- xy[root, ]
    rel <- sweep(xy[leaves, , drop = FALSE], 2, rootxy)
    tt <- sqrt(rel[, 1]^2 + rel[, 2]^2)
    phi <- atan2(rel[, 2], rel[, 1])

    e <- new.env()
    e$t <- tt
    e$phi <- phi
    e$x <- xy[leaves, 1]
    e$y <- xy[leaves, 2]
    e$w <- weight
    e$isleaf <- rep(TRUE, length(leaves))
    edges <- list()
    # leaves sitting on the root connect straight away
    active <- which(tt > .Machine$double.eps)
    for (i in which(tt <= .Machine$double.eps)) edges[[length(edges) + 1]] <- c(i, -1L)

    add_node <- function(t, phi) {
        e$t <- c(e$t, t)
        e$phi <- c(e$phi, phi)
        e$x <- c(e$x, rootxy[1] + t * cos(phi))
        e$y <- c(e$y, rootxy[2] + t * sin(phi))
        e$w <- c(e$w, NA_real_)
        e$isleaf <- c(e$isleaf, FALSE)
        length(e$t)
    }
    meeting <- function(i, j) {
        t1 <- e$t[i]
        t2 <- e$t[j]
        gap <- (e$phi[j] - e$phi[i]) %% (2 * pi)
        if (gap <= 0) gap <- gap + 2 * pi
        lnm <- (log(t1) + log(t2)) / 2 - gap / (2 * tan_a)
        tm <- exp(lnm)
        if (tm <= min(t1, t2)) {
            return(list(ok = TRUE, t = tm, phi = e$phi[i] + tan_a * (log(t1) - lnm), kind = "steiner"))
        }
        if (gap <= tan_a * abs(log(t2) - log(t1))) {
            if (t1 < t2) {
                return(list(ok = TRUE, t = t1, kind = "absorb", keep = i, drop = j))
            }
            return(list(ok = TRUE, t = t2, kind = "absorb", keep = j, drop = i))
        }
        list(ok = FALSE)
    }

    while (length(active) > 1) {
        ord <- active[order(e$phi[active])]
        na <- length(ord)
        best <- NULL
        for (a in seq_len(na)) {
            i <- ord[a]
            j <- ord[if (a < na) a + 1 else 1]
            mm <- meeting(i, j)
            if (mm$ok && (is.null(best) || mm$t > best$t)) best <- c(mm, list(i = i, j = j))
        }
        if (is.null(best)) break
        if (best$kind == "steiner") {
            s <- add_node(best$t, best$phi)
            e$w[s] <- e$w[best$i] + e$w[best$j]
            edges[[length(edges) + 1]] <- c(best$i, s)
            edges[[length(edges) + 1]] <- c(best$j, s)
            active <- c(setdiff(active, c(best$i, best$j)), s)
        } else {
            e$w[best$keep] <- e$w[best$keep] + e$w[best$drop]
            edges[[length(edges) + 1]] <- c(best$drop, best$keep)
            active <- setdiff(active, best$drop)
        }
    }
    for (i in active) edges[[length(edges) + 1]] <- c(i, -1L)

    list(env = e, edges = edges, rootxy = rootxy)
}

# Sample one tree edge (child -> parent) as n points along its log-spiral arc
# (straight to the root). Returns an n x 2 matrix.
.spiral_arc <- function(st, ed, n) {
    e <- st$env
    child <- ed[1]
    parent <- ed[2]
    u <- seq(0, 1, length.out = n)
    if (parent == -1L) {
        xs <- e$x[child] + u * (st$rootxy[1] - e$x[child])
        ys <- e$y[child] + u * (st$rootxy[2] - e$y[child])
    } else {
        dphi <- ((e$phi[parent] - e$phi[child] + pi) %% (2 * pi)) - pi
        lnt <- log(e$t[child]) + u * (log(e$t[parent]) - log(e$t[child]))
        ph <- e$phi[child] + u * dphi
        xs <- st$rootxy[1] + exp(lnt) * cos(ph)
        ys <- st$rootxy[2] + exp(lnt) * sin(ph)
    }
    cbind(xs, ys)
}

# Render each tree edge as a logarithmic-spiral arc.
.spiral_render <- function(st, n) {
    e <- st$env
    res <- vector("list", length(st$edges))
    for (k in seq_along(st$edges)) {
        arc <- .spiral_arc(st, st$edges[[k]], n)
        res[[k]] <- data.frame(x = arc[, 1], y = arc[, 2], flow = e$w[st$edges[[k]][1]], edge = k)
    }
    do.call(rbind, res)
}

# Approximate FLOWTREE refinement (Verbeek et al. 2011, section 5): tree nodes
# (leaves, root, join points) are held fixed so the planar skeleton is
# preserved; only the interior points of each edge are relaxed to (a) be smooth
# (FS) and (b) keep clear of node obstacles (Fobs), anchored to the spiral arc
# so the angle restriction is not badly violated.
.spiral_optimize <- function(st, n, iterations = 60, buffer = NULL,
                             w_smooth = 1, w_obs = 0.6, w_anchor = 0.15, lr = 0.25) {
    e <- st$env
    ntn <- length(e$t)
    rootid <- ntn + 1L
    X <- c(e$x, st$rootxy[1])
    Y <- c(e$y, st$rootxy[2])
    fixed <- rep(TRUE, ntn + 1L)

    obs_ids <- c(which(e$isleaf), rootid)
    ox <- X[obs_ids]
    oy <- Y[obs_ids]

    chains <- vector("list", length(st$edges))
    own_ends <- vector("list", 0)
    for (k in seq_along(st$edges)) {
        ed <- st$edges[[k]]
        parent <- if (ed[2] == -1L) rootid else ed[2]
        arc <- .spiral_arc(st, ed, n)
        if (n > 2) {
            ids <- length(X) + seq_len(n - 2)
            X <- c(X, arc[2:(n - 1), 1])
            Y <- c(Y, arc[2:(n - 1), 2])
            fixed <- c(fixed, rep(FALSE, n - 2))
            for (id in ids) own_ends[[id]] <- c(ed[1], parent)
            chains[[k]] <- c(ed[1], ids, parent)
        } else {
            chains[[k]] <- c(ed[1], parent)
        }
    }
    anchorX <- X
    anchorY <- Y
    Np <- length(X)
    nbr <- vector("list", Np)
    for (ch in chains) {
        for (i in seq_along(ch)) {
            if (i > 1) nbr[[ch[i]]] <- c(nbr[[ch[i]]], ch[i - 1])
            if (i < length(ch)) nbr[[ch[i]]] <- c(nbr[[ch[i]]], ch[i + 1])
        }
    }

    if (is.null(buffer)) {
        lx <- e$x[e$isleaf]
        ly <- e$y[e$isleaf]
        nn <- vapply(seq_along(lx), function(i) {
            d <- sqrt((lx[i] - lx[-i])^2 + (ly[i] - ly[-i])^2)
            if (length(d)) min(d) else 0
        }, numeric(1))
        buffer <- 0.5 * stats::median(nn)
    }

    movable <- which(!fixed)
    for (it in seq_len(iterations)) {
        nX <- X
        nY <- Y
        for (p in movable) {
            nb <- nbr[[p]]
            mvx <- w_smooth * (mean(X[nb]) - X[p]) + w_anchor * (anchorX[p] - X[p])
            mvy <- w_smooth * (mean(Y[nb]) - Y[p]) + w_anchor * (anchorY[p] - Y[p])
            dx <- X[p] - ox
            dy <- Y[p] - oy
            d <- sqrt(dx * dx + dy * dy)
            keep <- d < buffer & d > 1e-9 & !(obs_ids %in% own_ends[[p]])
            if (any(keep)) {
                f <- w_obs * (buffer - d[keep]) / buffer
                mvx <- mvx + sum(f * dx[keep] / d[keep])
                mvy <- mvy + sum(f * dy[keep] / d[keep])
            }
            nX[p] <- X[p] + lr * mvx
            nY[p] <- Y[p] + lr * mvy
        }
        X <- nX
        Y <- nY
    }

    res <- vector("list", length(chains))
    for (k in seq_along(chains)) {
        ch <- chains[[k]]
        res[[k]] <- data.frame(x = X[ch], y = Y[ch], flow = e$w[st$edges[[k]][1]], edge = k)
    }
    do.call(rbind, res)
}
