#' Differential expression analysis
#'
#' @docType methods
#' @name DEAnalyze
#' @rdname DEAnalyze
#'
#' @param obj Matrix like object or an ExprDataSet instance.
#' @param option "limma", "DESeq2", or others.
#' @param SampleAnn Matrix like object (only when obj is a matrix),
#' the first column should match colnames in obj, and the second column should be condition
#' @param type "Microarray" or "RNASeq", only needed when obj is not ExprDataSet instance.
#'
#' @return An ExprDataSet instance.
#' @seealso \code{\link{ExprDataSet-class}}
#'
#' @author Wubing Zhang
#'
#' @importFrom limma makeContrasts eBayes topTable contrasts.fit lmFit
#' @importFrom DESeq2 DESeqDataSetFromMatrix results DESeq
#' @export

DEAnalyze <- function(obj, option = "limma", SampleAnn = NULL, type = "Microarray"){
  if(is.matrix(obj) | is.data.frame(obj)){
    rownames(SampleAnn) = SampleAnn[,1]
    colnames(SampleAnn)[1:2] = c("Sample", "Condition")
    obj = as.matrix(obj[, SampleAnn$Sample])
    obj = ExprDataSet(profile = obj, sampleAnn = SampleAnn, type = type)
  }
  slot(obj, "Process") <- paste0(slot(obj, "Process"), "+", option)

  # Fill missing values in expression profile
  if(sum(is.na(slot(obj, "profile")))>0){
    expr = slot(obj, "profile")
    design = as.character(slot(obj, "sampleAnn")[, 2])
    for(cond in unique(design)){
      idx = design==cond
      tmp = expr[, idx, drop = FALSE]
      tmp[is.na(tmp)] = rowMeans(tmp, na.rm = TRUE)[which(is.na(tmp), arr.ind = TRUE)[,1]]
      expr[, idx] = tmp
    }
  }

  # Remove dropouts
  if(sum(rowSums(slot(obj, "profile"))<3)>10){
    expr = slot(obj, "profile")
    design = as.character(slot(obj, "sampleAnn")[, 2])
    tmp = sapply(unique(design), function(x){
      rowSums(expr[, design==x, drop = FALSE])>2
    })
    slot(obj, "profile") <- expr[rowSums(tmp)>1, ]
  }

  # Differential expression analysis
  if(tolower(option) == "limma"){
    requireNamespace("limma")
    design = model.matrix(~-1+Condition, slot(obj, "sampleAnn"))
    rownames(design) = colnames(slot(obj, "profile"))
    #"ls" for least squares or "robust" for robust regression
    fit = eBayes(lmFit(slot(obj, "profile"), design, na.rm=TRUE))
    res = topTable(fit, adjust.method="BH", coef=1, number = nrow(slot(obj, "profile")))
    res = res[, c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
    colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    slot(obj, "DEGRes") = res
  }

  if(tolower(option) == "deseq2"){
    requireNamespace("DESeq2")
    # DESeq2
    dds <- DESeq2::DESeqDataSetFromMatrix(slot(obj, "profile"), colData=slot(obj, "sampleAnn"), design=~Condition)
    dds <- DESeq2::DESeq(dds)
    res <- as.data.frame(DESeq2::results(dds, alpha = 0.1))
    res$padj[is.na(res$padj)] = 1
    res = res[, c("baseMean", "log2FoldChange", "stat", "pvalue", "padj")]
    colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    slot(obj, "DEGRes") = res
  }
  return(obj)
}

