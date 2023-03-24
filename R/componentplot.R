#' @name componentplot
#' @title Draws a Plot Showing the Number of Components as a Function of the Cutoff Value.
#' @description The function takes a square sequence similarity matrix and calculates
#' the number of independent components in which the sequence similarity between the
#' members is greater than a specified value (between 0 and 1). It creates a plots showing
#' the number of components corresponding to a given cutoff value. The function also
#' depicts the inflection point with a vertical line. Upper and lower bounds can be
#' provided between which the inflection point will be found (if it exists).
#'
#' Version 0.0.1.
#' Author: Dr. Matthew Cserhati
#' Email: csmatyi@protonmail.com
#' March 19, 2023
#'
#' @importFrom grDevices colorRampPalette dev.off jpeg
#' @importFrom graphics box legend
#' @importFrom inflection stats igraph
#' @importFrom graphics abline lines
#' @importFrom loess predict
#' @importFrom bese graph_from_adjacency_matrix
#'
#' @param mx a square sequence similarity matrix
#' @param lower_bound lower bounds for calculating the inflection point in, default = 0
#' @param upper_bound upper bounds for calculating the inflection point in, default = 1
#' @return A plot showing the number of components as a function of the cutoff threshold.
#'
#' @references Mullner. <ArXiv:1109.2378>
#' @references Cserhati, Carter. (2020, Journal of Creation 34(3):41-50), <https://dl0.creation.com/articles/p137/c13759/j34-3_64-73.pdf>
#'
#' @examples
#' componentplot(xenarthra,0.75,0.9)
#' componentplot(xenarthra)
#'
#'@export
utils::globalVariables(c("exit","abline","loess","lines","predict","bese","graph_from_adjacency_matrix"))
componentplot <- function(mx,lower_bound=1,upper_bound=1) {
  if (!is.numeric(mx)) {
    print("Data frame contains non-numeric values!")
    exit()
  }
  if (lower_bound < -1 | lower_bound > 1 | upper_bound < -1 | upper_bound > 1) {
    print("One of the boundaries is out of range (0-1)!")
    exit()
  }

  mxv <- mx[upper.tri(mx)]
  minv <- floor(10*min(mxv))/10
  xs <- c()
  xc <- seq(minv,1,0.01)
  for (i in xc) {
    mx2 <- mx
    mx2[mx2 >= i] <- 1
    mx2[mx2 < i] <- 0
    mg <- igraph::graph_from_adjacency_matrix(mx2)
    nt <- igraph::count_components(mg)
    xs <- c(xs,nt)
  }

  xc2 <- xc[xc >= lower_bound & xc <= upper_bound]
  xs2 <- xs[xc >= lower_bound & xc <= upper_bound]
  b <- inflection::bese(xc2,xs2,0)

  lx <- loess(xs~xc, span=0.1)
  xs_pred <- predict(lx,xc)

  p <- plot(xc,xs,pch=20,xlab="cutoff",ylab="no. components",cex=0.5)
  lines(xc,xs_pred,col="blue")
  if (!is.nan(b$iplast)) {
    if (lower_bound > minv) {abline(v=lower_bound,col="green")}
    if (upper_bound < 1) {abline(v=upper_bound,col="green")}
    abline(v=b$iplast,col="red")
  }
}
