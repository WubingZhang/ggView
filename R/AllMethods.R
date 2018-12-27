##' @method getDEG ExprDataSet
##' @export
getDEG.ExprDataSet <- function(x) {
  return(x@DEGRes)
}

##'  Coerce to a \code{DESeqDataSet} instance
##'
##' @name as.DESeqDataSet
##' @docType methods
##' @rdname as.DESeqDataSet
##'
##' @title Function to coerce a \code{ExprDataSet} object to a \code{DESeqDataSet} instance.
##' @param object A \code{ExprDataSet} instance.
##' @param ... additional parameter
##' @return A \code{DESeqDataSet} instance.
##' @export
##' @usage as.DESeqDataSet(object)
##' @author Wubing Zhang
as.DESeqDataSet.ExprDataSet <- function(object) {
  dds = DESeqDataSetFromMatrix(slot(object, "rawdata"),
        colData = slot(object, "SampleAnn"), design = ~Condition)
  return(dds)
}

##' @method dim ExprDataSet
##' @export
dim.ExprDataSet <- function(x) {
  dim(x@profile)
}

##' summary method for \code{ExprDataSet} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{ExprDataSet} instance.
##' @param ... additional parameter
##' @return A data frame
##' @exportMethod summary
##' @usage summary(object)
##' @author Wubing Zhang
setMethod("summary", signature(object="ExprDataSet"),
          function(object) {
            return(object@DEGRes)
          }
)
