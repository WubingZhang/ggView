##' Calculate score across genes and samples
##'
##' This wrapper function combines filtering out genes with low reads in a number of samples (recommended for limma:voom) with normalization
##' @param gm normalized count matrix; rows are all genes in the signature that shall be summarized into one score; columns are samples
##' @param gset Gene sets.
##' @param summarizationFunction ("PC", default), Pearson, ssGSEA or mean (other value)
##' @return numeric vector or matrix.
##' @author Wubing Zhang
##' @import GSVA
##' @export
gsScore <- function(gm, gset=rownames(gm), summarizationFunction="PC") {
  if(tolower(summarizationFunction) == "ssgsea"){
    gss <- gsva(gm, gset, method="ssgsea")
    return(gss)
  }
  gm = gm[rownames(gm)%in%gset, , drop = FALSE]
  if(nrow(gm)<3) summarizationFunction="mean"
  if (tolower(summarizationFunction) == "pc") {
    pc <- prcomp(t(gm), retx=TRUE)
    gss <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(gm)))
  } else if (tolower(summarizationFunction) == "pearson"){
    corvalues = cor(t(gm))
    diag(corvalues) = 0
    corvalues[corvalues<0] = 0
    w = rowSums(corvalues) / sum(corvalues)
    gss <- (t(gm) %*% w)[,1]
  } else {
    gss <- colMeans(gm)
  }
  return(gss)
}
