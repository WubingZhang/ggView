#' Call differential expressed genes using limma
#'
#' @docType methods
#' @name runLimma
#' @rdname runLimma
#'
#' @param mat Matrix like object.
#' @param idx_ctrl Index of control samples.
#' @param idx_treat Index of Treatment samples.
#' @param pvalueCutoff Pvalue cutoff.
#'
#' @return A list contains limma result table and ggplot object (volcano figure).
#'
#' @author Wubing Zhang
#'
#' @import limma
#' @export

runLimma <- function(mat, idx_ctrl, idx_treat, pvalueCutoff = 0.25){
  tmp1 = mat[, idx_ctrl]
  tmp1[is.na(tmp1)] = rowMeans(tmp1, na.rm = TRUE)[which(is.na(tmp1), arr.ind = TRUE)[,1]]
  tmp2 = mat[, idx_treat]
  tmp2[is.na(tmp2)] = rowMeans(tmp2, na.rm = TRUE)[which(is.na(tmp2), arr.ind = TRUE)[,1]]
  mat = cbind(tmp1, tmp2)
  mat = mat[unique(which(!is.na(mat), arr.ind = TRUE)[,1]),]
  labels = factor(rep(c("Ctrl", "Treat"), c(length(idx_ctrl),length(idx_treat))))
  design = model.matrix(~-1+labels)
  colnames(design) <- levels(labels)
  rownames(design) = colnames(mat)
  contrast <- makeContrasts(TreatvsCtrl = Treat - Ctrl, levels = design)
  #"ls" for least squares or "robust" for robust regression
  fit = eBayes(contrasts.fit(lmFit(mat, design, na.rm=TRUE), contrast))
  res = topTable(fit, adjust.method="BH", coef="TreatvsCtrl", number = nrow(mat))
  res$nLogP = -log10(res$adj.P.Val)
  p = VolcanoView(res, top = 10, width = 5, height = 4, y_cutoff = pvalueCutoff)
  res = list(res = res, p=p)
  return(res)
}

