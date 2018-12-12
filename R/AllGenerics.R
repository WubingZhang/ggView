#' Get summary of differential expression results
#'
#' @param x ExprDataSet object
#' @return A data frame which include differential expression analysis results
#' @export
getDEG <- function(x) {
  UseMethod("getDEG", x)
}
