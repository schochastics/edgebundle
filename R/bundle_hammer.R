#' @title hammer edge bundling
#' @description Implements the hammer edge bundling by Ian Calvert.
#' @details This function only wraps existing python code from the datashader library. Original code can be found at https://gitlab.com/ianjcalvert/edgehammer.
#' Datashader is a huge library with a lot of dependencies, so think twice if you want to install it just for edge bundling.
#' Check https://datashader.org/user_guide/Networks.html for help concerning parameters bw and decay.
#' To install all dependencies, use [install_bundle_py].
#' @param object a graph object (igraph/network/tbl_graph)
#' @param xy coordinates of vertices
#' @param bw bandwidth parameter
#' @param decay decay parameter
#' @return data.frame containing the bundled edges
#' @author David Schoch
#' @details see [online](https://github.com/schochastics/edgebundle) for plotting tips
#' @seealso [edge_bundle_force],[edge_bundle_stub], [edge_bundle_path]
#' @export

edge_bundle_hammer <- function(object,xy,bw=0.05,decay=0.7){
  if (!requireNamespace('reticulate', quietly = TRUE)) {
    stop('The `reticulate` package is required for this functionality')
  }
  if(any(class(object)=="igraph")){
    if (!requireNamespace('igraph', quietly = TRUE)) {
      stop('The `igraph` package is required for this functionality')
    }
    nodes <- data.frame(name=paste0("node",0:(igraph::vcount(object)-1)),x=xy[,1],y=xy[,2])
    el <- igraph::get.edgelist(object,names = FALSE)
    el1 <- data.frame(source=el[,1]-1,target=el[,2]-1)

  } else if(any(class(object)=="tbl_graph")){
    if (!requireNamespace('tidygraph', quietly = TRUE)) {
      stop('The `tidygraph` package is required for this functionality')
    }
    object <- tidygraph::as.igraph(object)
    nodes <- data.frame(name=paste0("node",0:(igraph::vcount(object)-1)),x=xy[,1],y=xy[,2])
    el <- igraph::get.edgelist(object,names = FALSE)
    el1 <- data.frame(source=el[,1]-1,target=el[,2]-1)

  } else if(any(class(object)=="network")){
    nodes <- data.frame(name=paste0("node",0:(network::get.network.attribute(object,"n")-1)),x=xy[,1],y=xy[,2])
    el <- network::as.edgelist(object)
    el1 <- data.frame(source=el[,1]-1,target=el[,2]-1)
  } else{
    stop("only `igraph`, `network` or `tbl_graph` objects supported.")
  }
  data_bundle <- shader_env$datashader_bundling$hammer_bundle(nodes,el1,initial_bandwidth = bw,decay = decay)
  data_bundle$group <- is.na(data_bundle$y)+0
  data_bundle$group <- cumsum(data_bundle$group)+1
  data_bundle <- data_bundle[!is.na(data_bundle$y),]
  data_bundle$index <- unlist(sapply(table(data_bundle$group),function(x) seq(0,1,length.out=x)))
  data_bundle[,c("x","y","index","group")]
}

#' @title install python dependencies for hammer bundling
#' @description install datashader and scikit-image
#' @param method Installation method (by default, "auto" automatically finds a
#' method that will work in the local environment, but note that the
#' "virtualenv" method is not available on Windows)
#' @param conda Path to conda executable (or "auto" to find conda using the PATH
#' and other conventional install locations)
#' @export
#'
install_bundle_py <- function(method = "auto", conda = "auto") {
  if (!requireNamespace('reticulate', quietly = TRUE)) {
    stop('The `reticulate` package is required for this functionality')
  }
  reticulate::py_install("datashader", method = method, conda = conda, pip = TRUE)
  reticulate::py_install("scikit-image", method = method, conda = conda, pip = TRUE)
}

# Environment for globals
shader_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname)
  assign("datashader_bundling", reticulate::import("datashader.bundling", delay_load = TRUE), shader_env)
}
