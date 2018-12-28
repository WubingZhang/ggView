#' Differential expression analysis
#'
#' @docType methods
#' @name DEAnalyze
#' @rdname DEAnalyze
#'
#' @param obj Matrix like object or an ExprDataSet instance.
#' @param SampleAnn Matrix like object (only when obj is a matrix),
#' the rownames should match colnames in obj, and the first column should be Condition.
#' @param type "Array" or "RNASeq", only needed when obj is matrix like object.
#' @param minS Same as that in `TransformCount`.
#' @param trans.method Same as `method` in `TransformCount`.
#'
#' @return An ExprDataSet instance.
#' @seealso \code{\link{ExprDataSet-class}}
#'
#' @author Wubing Zhang
#'
#' @importFrom limma eBayes topTable lmFit
#' @importFrom DESeq2 DESeqDataSetFromMatrix results DESeq
#' @export

DEAnalyze <- function(obj, SampleAnn = NULL, type = "Array",
                        minS = 2, trans.method = "vst"){
  if(is.matrix(obj) | is.data.frame(obj)){
    colname(SampleAnn)[1] = "Condition"
    expr <- as.matrix(obj[, rownames(SampleAnn)])
    obj = new("ExprDataSet", rawdata = expr, SampleAnn = SampleAnn, type = type)
  }

  # Fill missing values in expression profile
  if(sum(is.na(slot(obj, "rawdata")))>0){
    design = as.vector(slot(obj, "SampleAnn")[, 1])
    for(cond in unique(design)){
      idx = design==cond
      tmp = slot(obj, "rawdata")[, idx, drop = FALSE]
      tmp[is.na(tmp)] = rowMeans(tmp, na.rm = TRUE)[which(is.na(tmp), arr.ind = TRUE)[,1]]
      slot(obj, "rawdata")[, idx] = tmp
    }
  }

  # Normalizing dataset
  if(tolower(type) == "array"){
    slot(obj, "normlized") = limma::normalizeQuantiles(slot(obj, "rawdata"))
  }else {
    # Remove dropouts
    sel <- rowSums(slot(obj, "rawdata") >= 0) >= minS
    slot(obj, "rawdata") <- slot(obj, "rawdata")[sel, ]
    slot(obj, "normlized") = TransformCount(slot(obj, "rawdata"), method = trans.method)
  }

  # Differential expression analysis
  if(tolower(type) == "array"){
    requireNamespace("limma")
    design = model.matrix(~-1+Condition, slot(obj, "SampleAnn"))
    rownames(design) = colnames(slot(obj, "normlized"))
    design[,1] = 1
    #"ls" for least squares or "robust" for robust regression
    fit = eBayes(lmFit(slot(obj, "normlized"), design, na.rm=TRUE))
    res = topTable(fit, adjust.method="BH", coef=2, number = nrow(slot(obj, "normlized")))
    res = res[, c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
    colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
  }else if(tolower(type) == "rnaseq"){
    requireNamespace("DESeq2")
    # DESeq2
    dds = as.DESeqDataSet(obj)
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::lfcShrink(dds, coef = 2, quiet = TRUE)
    res$padj[is.na(res$padj)] = 1
    res = res[, c("baseMean", "log2FoldChange", "stat", "pvalue", "padj")]
    colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
  }else{
    stop("Data type error! ")
  }
  slot(obj, "DEGRes") = as.data.frame(res)
  return(obj)
}

