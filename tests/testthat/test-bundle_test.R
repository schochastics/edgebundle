test_that("force bundle works", {
  g  <- igraph::graph_from_edgelist(matrix(c(1,12,2,11,3,10,4,9,5,8,6,7),ncol = 2,byrow = TRUE),FALSE)
  xy <- cbind(c(rep(0,6),rep(1,6)),c(1:6,1:6))
  expect_silent(dat <- edge_bundle_force(g,xy))
  expect_type(dat, "list")
})

test_that("stub bundle works", {
   g <- igraph::graph.star(10,"undirected")

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

  expect_silent(dat <- edge_bundle_stub(g,xy))
  expect_type(dat, "list")
})
