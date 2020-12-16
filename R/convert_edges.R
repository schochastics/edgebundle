#' @title Convert edges
#' @description converts edges of an igraph/network/tidygraph object into format useable for edge bundling
#' @param object graph object
#' @param coords coordinates of vertices
#' @return data frame of edges with coordinates
#' @author David Schoch
#' @export
#'
convert_edges <- function(object,coords) UseMethod("convert_edges")

#' @rdname convert_edges
#' @method convert_edges default
#' @export
convert_edges.default <- function(object,coords){
  stop("don't know how to handle class ", dQuote(data.class(object)))
}

#' @rdname convert_edges
#' @method convert_edges igraph
#' @export
convert_edges.igraph <- function(object,coords){
  if (!requireNamespace('igraph', quietly = TRUE)) {
    stop('The `igraph` package is required for this functionality')
  }
  el <- igraph::as_edgelist(object,names = FALSE)
  if(igraph::vcount(object)!=nrow(coords)){
    stop('number of rows in `coords` does not match number of vertices')
  }
  edges_xy <- cbind(coords[el[,1],1],coords[el[,1],2],coords[el[,2],1],coords[el[,2],2])
  edges_xy
}

#' @rdname convert_edges
#' @method convert_edges network
#' @export
convert_edges.network <- function(object,coords){
  if (!requireNamespace('network', quietly = TRUE)) {
    stop('The `network` package is required for this functionality')
  }
  el <- network::as.edgelist(object)
  if(network::get.network.attribute(object,"n")!=nrow(coords)){
    stop('number of rows in `coords` does not match number of vertices')
  }
  edges_xy <- cbind(coords[el[,1],1],coords[el[,1],2],coords[el[,2],1],coords[el[,2],2])
  edges_xy
}


#' @rdname convert_edges
#' @method convert_edges tbl_graph
#' @export
convert_edges.tbl_graph <- function(object,coords){
  if (!requireNamespace('tidygraph', quietly = TRUE)) {
    stop('The `tidygraph` package is required for this functionality')
  }
  el <- as.matrix(tidygraph::as_tibble(object,"edges"))
  edges_xy <- cbind(coords[el[,1],1],coords[el[,1],2],coords[el[,2],1],coords[el[,2],2])
  edges_xy
}
