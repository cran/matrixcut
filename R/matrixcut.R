#' @name matrixcut
#' @title Calculates The Component Membership in a Sequence Similarity Matrix at a given Cutoff Value
#' @description The function takes a square sequence similarity matrix and calculates
#' those independent components in which the sequence similarity between the members
#' is greater than a specified value (between 0 and 1). The result provied by the function
#' is a list of species with their component membership.
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
#' @param cut a given cutoff value to calculate components for, default value = -1
#' @return The inflection point: cutoff value for the optimal number of clusters
#' @return A list of the species together with their cluster membership.
#'
#' @references Mullner. <ArXiv:1109.2378>
#' @references Cserhati, Carter. (2020, Journal of Creation 34(3):41-50), <https://dl0.creation.com/articles/p137/c13759/j34-3_64-73.pdf>
#'
#' @examples
#' matrixcut(primates,0.8)
#' matrixcut(primates)
#'
#'@export
utils::globalVariables(c("exit","abline","loess","lines","predict","bese","graph_from_adjacency_matrix"))
matrixcut <- function(mx, cut=-1) {
  if (!is.numeric(mx)) {
    print("Data frame contains non-numeric values!")
    exit()
  }
  if ((cut < 0 | cut > 1) && (cut != -1)) {
    print("Cutoff value is out of range!")
    exit()
  }

  inflection_point <- NULL
  if (cut == -1) {
    mxv <- mx[upper.tri(mx)]
    minv <- floor(10*min(mxv))/10
    xs <- c()
    xc <- seq(minv,1,0.001)
    for (i in xc) {
      nt <- sum(as.integer(mxv >= i))
      xs <- c(xs,nt)
    }
    b <- inflection::bese(xc,xs,0)
    inflection_point <- b$iplast
  } else if (cut >= 0 && cut <= 1) {
    inflection_point <- cut
  } else {
    print("Invalid cutoff value!")
    exit()
  }

  mx2 <- mx
  mx2[mx2 >= inflection_point] <- 1
  mx2[mx2 < inflection_point] <- 0
  mg <- igraph::graph_from_adjacency_matrix(mx2)
  mc <- igraph::components(mg)

  mc$membership
}
