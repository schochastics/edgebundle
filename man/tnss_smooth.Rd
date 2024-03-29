% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flow_tnss.R
\name{tnss_smooth}
\alias{tnss_smooth}
\title{Smooth a Steiner tree}
\usage{
tnss_smooth(g, bw = 3, n = 10)
}
\arguments{
\item{g}{Steiner tree computed with \link{tnss_tree}}

\item{bw}{bandwidth of Gaussian Kernel}

\item{n}{number of extra nodes to include per edge}
}
\value{
data.frame containing the smoothed paths
}
\description{
Converts the Steiner tree to smooth paths
}
\details{
see see \href{https://github.com/schochastics/edgebundle}{online} for tips on plotting the result
}
\examples{
xy <- cbind(state.center$x, state.center$y)[!state.name \%in\% c("Alaska", "Hawaii"), ]
xy_dummy <- tnss_dummies(xy, root = 4)
gtree <- tnss_tree(cali2010, xy, xy_dummy, root = 4, gamma = 0.9)
tree_smooth <- tnss_smooth(gtree, bw = 10, n = 10)
}
\author{
David Schoch
}
