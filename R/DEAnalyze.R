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
#' @param method Differential expression analysis method, e.g. limma, DESeq2, GFOLD.
#' @param app.dir The path to application (e.g. GFOLD).
#'
#' @return An ExprDataSet instance.
#' @seealso \code{\link{ExprDataSet-class}}
#'
#' @author Wubing Zhang
#'
#' @importFrom limma eBayes topTable lmFit
#' @importFrom DESeq2 DESeqDataSetFromMatrix results DESeq
#' @importFrom edgeR filterByExpr
#' @export

DEAnalyze <- function(obj, SampleAnn = NULL, type = "Array", method = "limma",
                      app.dir = "/Users/Wubing/Applications/gfold/gfold"){
  requireNamespace("edgeR")
  ## Create a new object ##
  if(is.matrix(obj) | is.data.frame(obj)){
    colnames(SampleAnn)[1] = "Condition"
    expr <- as.matrix(obj[, rownames(SampleAnn)])
    obj = new("ExprDataSet", rawdata = expr, SampleAnn = SampleAnn, type = type)
  }

  ## Fill missing values in expression profile ##
  if(sum(is.na(slot(obj, "rawdata")))>0){
    design = as.vector(slot(obj, "SampleAnn")[, 1])
    for(cond in unique(design)){
      idx = design==cond
      tmp = slot(obj, "rawdata")[, idx, drop = FALSE]
      tmp[is.na(tmp)] = rowMeans(tmp, na.rm = TRUE)[which(is.na(tmp), arr.ind = TRUE)[,1]]
      slot(obj, "rawdata")[, idx] = tmp
    }
  }

  ## Differential expression analysis ##
  design = model.matrix(~-1+Condition, slot(obj, "SampleAnn"))
  rownames(design) = colnames(slot(obj, "rawdata"))
  design[,1] = 1

  if(tolower(type) == "array"){
    requireNamespace("limma")
    slot(obj, "normlized") = normalizeQuantiles(slot(obj, "rawdata"))
    #"ls" for least squares or "robust" for robust regression
    fit = eBayes(lmFit(slot(obj, "normlized"), design, na.rm=TRUE))
    res = topTable(fit, adjust.method="BH", coef=2, number = nrow(slot(obj, "normlized")))
    res = res[, c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
    colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
  }else if(tolower(type) == "rnaseq"){
    idx_r <- filterByExpr(slot(obj, "rawdata"), design)
    slot(obj, "rawdata") = slot(obj, "rawdata")[idx_r, ]
    if(method == "DESeq2"){
      requireNamespace("DESeq2")
      slot(obj, "normlized") = TransformCount(slot(obj, "rawdata"), method = "vst")
      # DESeq2
      dds = as.DESeqDataSet(obj)
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::lfcShrink(dds, coef = 2, quiet = TRUE)
      res$padj[is.na(res$padj)] = 1
      res = res[, c("baseMean", "log2FoldChange", "stat", "pvalue", "padj")]
      colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    }
    if(method == "limma"){
      requireNamespace("limma")
      slot(obj, "normlized") = TransformCount(slot(obj, "rawdata"), method = "voom")
      # limma:voom
      dge <- DGEList(counts=slot(obj, "rawdata"))
      dge <- calcNormFactors(dge)
      dge <- voom(dge, design, plot=FALSE)
      fit <- eBayes(lmFit(dge, design))
      res = topTable(fit, adjust.method="BH", coef=2, number = nrow(slot(obj, "rawdata")))
      res = res[, c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
      colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    }
    if(method == "gfold"){
      slot(obj, "normlized") = TransformCount(slot(obj, "rawdata"), method = "voom")
      # GFOLD
      tmp = mapply(function(x){
        write.table(cbind(NA, slot(obj, "rawdata")[,x], NA, NA),
                    file=paste0(colnames(slot(obj, "rawdata"))[x], ".txt"),
                    sep="\t", col.names=FALSE)}, x=1:ncol(slot(obj, "rawdata")))
      lev = levels(slot(obj, "SampleAnn")$Condition)
      ctrlname = rownames(slot(obj, "SampleAnn"))[slot(obj, "SampleAnn")$Condition==lev[1]]
      treatname = rownames(slot(obj, "SampleAnn"))[slot(obj, "SampleAnn")$Condition==lev[2]]
      system(paste0(app.dir, " diff -s1 ", paste0(ctrlname, collapse=","),
                    " -s2 ", paste0(treatname, collapse=","), " -suf .txt -o gfold_tmp")
             )
      res = read.table("gfold_tmp", row.names=1, stringsAsFactors = FALSE)
      res = res[, c(5, 4, 2, 3, 3)]
      colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
      tmp = file.remove(paste0(ctrlname, ".txt"), paste0(treatname, ".txt"),
                        "gfold_tmp", "gfold_tmp.ext")
    }
  }else{
    stop("Data type error! ")
  }
  slot(obj, "DEGRes") = as.data.frame(res)
  return(obj)
}

