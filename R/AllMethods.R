##' @method getDEG ExprDataSet
##' @export
getDEG.ExprDataSet <- function(x) {
  return(x@DEGRes)
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
