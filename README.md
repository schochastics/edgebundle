
<!-- README.md is generated from README.Rmd. Please edit that file -->

# edgebundle

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/edgebundle)](https://CRAN.R-project.org/package=edgebundle)
<!-- badges: end -->

An R package that implements several edge bundling algorithms. So far it
includes

  - force directed edge bundling
    ([paper](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.212.7989&rep=rep1&type=pdf))
  - stub bundling
    ([paper](https://www.uni-konstanz.de/mmsp/pubsys/publishedFiles/NoBr13.pdf))
  - hammer bundling ([python
    code](https://datashader.org/_modules/datashader/bundling.html))

(*part of this package will eventually migrate to ggraph*)

## Installation

You can install the dev version of edgebundle with:

``` r
# install.packages("remotes")
remotes::install_github("schochastics/edgebundle")
```

Note that `edgebundle` imports `reticulate` and uses a pretty big python
library (datashader).

## Usage

All functions return a data frame of points along the edges of the
network.

``` r
library(igraph)
g <- graph_from_edgelist(
  matrix(c(1,12,2,11,3,10,4,9,5,8,6,7),ncol=2,byrow = T),F)
xy <- cbind(c(rep(0,6),rep(1,6)),c(1:6,1:6))

fbundle <- edge_bundle_force(g,xy,compatibility_threshold = 0.1)
head(fbundle)
#>            x       y     index group
#> 1 0.00000000 1.00000 0.0000000     1
#> 2 0.00611816 1.19977 0.0303030     1
#> 3 0.00987237 1.29767 0.0606061     1
#> 4 0.01929293 1.52427 0.0909091     1
#> 5 0.02790686 1.68643 0.1212121     1
#> 6 0.03440142 1.81285 0.1515152     1
```

The result can be visualized with ggplot.

``` r
library(ggplot2)

ggplot(fbundle)+
  geom_path(aes(x,y,group=group,col=as.factor(group)),size = 2,show.legend = FALSE)+
  geom_point(data=as.data.frame(xy),aes(V1,V2),size=5)+
  theme_void()
```

<img src="man/figures/README-plot-1.png" width="100%" />

For `edge_bundle_stub()`, you need `geom_bezier()` from the package
`ggforce`.

``` r
library(ggforce)
g <- graph.star(10,"undirected")

xy <- matrix(c(
  0,0,
  cos(90*pi/180),sin(90*pi/180),
  cos(80*pi/180),sin(80*pi/180),
  cos(70*pi/180),sin(70*pi/180),
  cos(330*pi/180),sin(330*pi/180),
  cos(320*pi/180),sin(320*pi/180),
  cos(310*pi/180),sin(310*pi/180),
  cos(210*pi/180),sin(210*pi/180),
  cos(200*pi/180),sin(200*pi/180),
  cos(190*pi/180),sin(190*pi/180)
),ncol=2,byrow=TRUE)

sbundle <- edge_bundle_stub(g,xy,beta = 90)

ggplot(sbundle)+
  geom_bezier(aes(x,y,group=group),size=2,col="grey66")+
  geom_point(data=as.data.frame(xy),aes(V1,V2),size=5)+
  theme_void()
```

<img src="man/figures/README-bezier-1.png" width="100%" />

## Showcase (US flight dataset)

*(The dataset is included in the package)*

![](man/figures/flights_fdeb.png)

![](man/figures/flights_heb.png)

![](man/figures/flights_seb.png) Code:

``` r
g <- us_flights
xy <- cbind(V(g)$longitude,V(g)$latitude)
verts <- data.frame(x=V(g)$longitude,y=V(g)$latitude)
fbundle <- edge_bundle_force(g,xy,compatibility_threshold = 0.6)
sbundle <- edge_bundle_stub(g,xy)
hbundle <- edge_bundle_hammer(g,xy,bw = 0.7,decay = 0.5)

states <- map_data("state")


p1 <- ggplot()+
  geom_polygon(data=states,aes(long,lat,group=group),col="white",size=0.1,fill=NA)+
  geom_path(data = fbundle,aes(x,y,group=group),col="#9d0191",size=0.05)+
  geom_path(data = fbundle,aes(x,y,group=group),col="white",size=0.005)+
  geom_point(data = verts,aes(x,y),col="#9d0191",size=0.25)+
  geom_point(data = verts,aes(x,y),col="white",size=0.25,alpha=0.5)+
  geom_point(data=verts[verts$name!="",],aes(x,y), col="white", size=3,alpha=1)+
  labs(title="Force Directed Edge Bundling")+
  ggraph::theme_graph(background = "black")+
  theme(plot.title = element_text(color="white"))

p2 <- ggplot()+
  geom_polygon(data=states,aes(long,lat,group=group),col="white",size=0.1,fill=NA)+
  geom_path(data = hbundle,aes(x,y,group=group),col="#9d0191",size=0.05)+
  geom_path(data = hbundle,aes(x,y,group=group),col="white",size=0.005)+
  geom_point(data = verts,aes(x,y),col="#9d0191",size=0.25)+
  geom_point(data = verts,aes(x,y),col="white",size=0.25,alpha=0.5)+
  geom_point(data=verts[verts$name!="",],aes(x,y), col="white", size=3,alpha=1)+
  labs(title="Hammer Edge Bundling")+
  ggraph::theme_graph(background = "black")+
  theme(plot.title = element_text(color="white"))

alpha_fct <- function(x,b=0.01,p=5,n=20){
  (1-b)*(2/(n-1))^p * abs(x-(n-1)/2)^p+b
}

p3 <- ggplot()+
  geom_polygon(data=states,aes(long,lat,group=group),col="white",size=0.1,fill=NA)+
  ggforce::geom_bezier(data = sbundle,aes(x,y,group=group,alpha=alpha_fct(..index..*20)),n=20,col="#9d0191",size=0.1,show.legend = FALSE)+
  ggforce::geom_bezier(data = sbundle,aes(x,y,group=group,alpha=alpha_fct(..index..*20)),n=20,col="white",size=0.01,show.legend = FALSE)+
  geom_point(data = verts,aes(x,y),col="#9d0191",size=0.25)+
  geom_point(data = verts,aes(x,y),col="white",size=0.25,alpha=0.5)+
  geom_point(data=verts[verts$name!="",],aes(x,y), col="white", size=3,alpha=1)+
  labs(title="Stub Edge Bundling")+
  ggraph::theme_graph(background = "black")+
  theme(plot.title = element_text(color="white"))
```

## Disclaimer

Edge bundling is able to produce neat looking network visualizations.
However, they do not necessarily enhance readability. After
experimenting with several algorithms, it became also quite evident that
the algorithms are very sensitive to the parameter settings. Consult the
original literature (if they even provide any guidelines) or experiment
yourself and **do not expect any miracles**.
