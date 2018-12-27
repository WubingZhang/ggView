#' Get summary of differential expression results
#'
#' @param x ExprDataSet object
#' @return A data frame which include differential expression analysis results
#' @export
getDEG <- function(x) {
  UseMethod("getDEG", x)
}

#'  Coerce to a \code{DESeqDataSet} instance
#'
#' @name as.DESeqDataSet
#' @docType methods
#' @rdname as.DESeqDataSet
#'
#' @title Function to coerce a \code{ExprDataSet} object to a \code{DESeqDataSet} instance.
#' @param object A \code{ExprDataSet} instance.
#' @param ... additional parameter
#' @return A \code{DESeqDataSet} instance.
#' @export
#' @usage as.DESeqDataSet(object)
#' @author Wubing Zhang
as.DESeqDataSet <- function(x) {
  UseMethod("as.DESeqDataSet", x)
}
