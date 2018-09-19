#' Differential expression analysis of square groups.
#'
#' @docType methods
#' @name squareDiff
#' @rdname squareDiff
#'
#' @param expr Matrix like object, each row is gene and each column is sample.
#' @param x Specifying x-axis which should be valid row names.
#' @param y Specifying y-axis which should be valid row names.
#' @param x_cutoff Specifying cutoff to classify x into high and low groups.
#' @param y_cutoff Specifying cutoff to classify x into high and low groups.
#'
#' @return A list.
#'
#' @author Wubing Zhang
#' @import ggplot2
#' @export

squareDiff <- function(expr, x = "TRAF3", y = c("TAP1", "B2M", "HLA-A", "HLA-B", "HLA-C"),
                       x_cutoff = c(-0.5, 0.5), y_cutoff = c(-0.5, 0.5), ...){
  require(ggplot2)
  # View the distribution
  p = SquareView(t(expr), x, y, x_cutoff = x_cutoff, y_cutoff = y_cutoff, option = 2, ...)
  gg = data.frame(x = colMeans(expr[x,]), y = colMeans(expr[y,]))

  # Select four group of samples
  group1_idx = which(gg$x > x_cutoff[2] & gg$y > y_cutoff[2])
  group2_idx = which(gg$x > x_cutoff[2] & gg$y < y_cutoff[2])
  group3_idx = which(gg$x < x_cutoff[1] & gg$y > y_cutoff[1])
  group4_idx = which(gg$x < x_cutoff[1] & gg$y < y_cutoff[1])

  # Do the comparison
  DP_xPyN = runLimma(expr, group2_idx, group1_idx, pvalueCutoff = 0.05)
  DP_xNyP = runLimma(expr, group4_idx, group1_idx, pvalueCutoff = 0.05)
  xPyN_DN = runLimma(expr, group3_idx, group2_idx, pvalueCutoff = 0.05)
  xNyP_DN = runLimma(expr, group3_idx, group4_idx, pvalueCutoff = 0.05)
  return(list(ditr_p = p, DP_xPyN = DP_xPyN, DP_xNyP = DP_xNyP, xPyN_DN = xPyN_DN, xNyP_DN = xNyP_DN))
}
