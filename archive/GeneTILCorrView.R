#' Dot heatmap
#'
#' @docType methods
#' @name GeneTILCorrView
#' @rdname GeneTILCorrView
#'
#' @param gene A vector of gene names.
#' @param cibersort The path to cibersort Rdataset.
#' @param exprdir The directory of processed TCGA expression profiles.
#' @return An object created by \code{ggplot}, which can be assigned
#' and further customized.
#'
#' @author Wubing Zhang
#' @import reshape2
#' @import ggplot2
#' @export
#'
GeneTILCorrView <- function(gene, cibersort = "TCGA_CIBERSORT.rds", exprdir = "ProcessedDat/"){
  CancerType = list.files(exprdir)
  TCGA_CIBERSORT = readRDS(cibersort)
  TRAF3_expr = data.frame(Sample=c(), Expr=c(), Cancer=c(), stringsAsFactors = FALSE)
  for(CT in CancerType){
    path = paste0(exprdir, CT, "/", CT, ".RNAseq.TPM.rds")
    tmp = log2(readRDS(path)+1)
    tmpExpr = data.frame(Sample=colnames(tmp), Expr=colMeans(tmp[gene, , drop = FALSE]),
                         Cancer=CT, stringsAsFactors = FALSE)
    TRAF3_expr = rbind.data.frame(TRAF3_expr, tmpExpr)
  }
  TRAF3_expr = TRAF3_expr[TRAF3_expr$Sample%in%rownames(TCGA_CIBERSORT), ]
  TRAF3_expr = cbind.data.frame(TRAF3_expr, TCGA_CIBERSORT[TRAF3_expr$Sample, ])
  corvalues = sapply(unique(TRAF3_expr$Cancer), function(x){
    idx = TRAF3_expr$Cancer==x
    sapply(colnames(TCGA_CIBERSORT)[1:22], function(y){
      cor(TRAF3_expr$Expr[idx], TRAF3_expr[idx, y])
    })
  })
  corvalues[is.na(corvalues)] = 0
  corvalues = as.data.frame(corvalues)
  p = dotHeatView(corvalues)
  p = p + labs(color = "R")
  p
}
