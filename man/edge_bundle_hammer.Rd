% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bundle_hammer.R
\name{edge_bundle_hammer}
\alias{edge_bundle_hammer}
\title{hammer edge bundling}
\usage{
edge_bundle_hammer(object, xy, bw = 0.05, decay = 0.7)
}
\arguments{
\item{object}{a graph object (igraph/network/tbl_graph)}

\item{xy}{coordinates of vertices}

\item{bw}{bandwidth parameter}

\item{decay}{decay parameter}
}
\value{
data.frame containing the bundled edges
}
\description{
Implements the hammer edge bundling by Ian Calvert.
}
\details{
This function only wraps existing python code from the datashader library. Original code can be found at https://gitlab.com/ianjcalvert/edgehammer.
Datashader is a huge library with a lot of dependencies, so think twice if you want to install it just for edge bundling.
Check https://datashader.org/user_guide/Networks.html for help concerning parameters bw and decay.
To install all dependencies, use \link{install_bundle_py}.

see \href{https://github.com/schochastics/edgebundle}{online} for plotting tips
}
\seealso{
\link{edge_bundle_force},\link{edge_bundle_stub}, \link{edge_bundle_path}
}
\author{
David Schoch
}
