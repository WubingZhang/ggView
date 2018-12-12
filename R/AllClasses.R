##' Class "ExprDataSet"
##' This class represents the datasets of expression profiles.
##'
##'
##' @name ExprDataSet-class
##' @aliases ExprDataSet-class
##'   show,ExprDataSet-method summary,ExprDataSet-method
##'
##' @docType class
##' @slot profile Expression profile.
##' @slot normalizedDS Normalized dataset.
##' @slot sampleAnn For matrix input: a DataFrame with at least two columns.
##' `sampleName` and `condition`.
##' @slot geneTable A DataFrame with gene annotations.
##' @slot type "Microarray" or "RNASeq".
##' @slot GPL GPL ID.
##' @slot Process Description of normalization process.
##' @slot DEGRes Summary result of differential expressed genes.
##' @exportClass ExprDataSet
##' @author WubingZhang
##' @keywords classes
ExprDataSet <- setClass("ExprDataSet",
         representation=representation(
           profile        = "matrix",
           normlizedDS    = "matrix",
           sampleAnn    = "data.frame",
           geneTable      = "data.frame",
           type           = "character",
           GPL            = "character",
           Process        = "character",
           DEGRes         = "data.frame"
        ),
        prototype=prototype(type      = "RNASeq",
                            Process   = "Normalized"),
        validity=function(object)
        {
          if(!all(object@sampleAnn[,1] %in% colnames(object@profile))){
            return("Sample names in sampleAnn do not match colnames of expression profile.")
          }
          return(TRUE)
        }
)

