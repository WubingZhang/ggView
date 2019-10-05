#' Differential expression analysis
#'
#' @docType methods
#' @name DEAnalyze
#' @rdname DEAnalyze
#'
#' @param obj Matrix like object or an ExprDataSet instance.
#' @param SampleAnn Matrix like object (only when obj is a matrix),
#' the rownames should match colnames in obj, and the first column should be Condition.
#' @param type "Array", "RNASeq" or "msms", only needed when obj is matrix like object.
#' @param method Differential expression analysis method, e.g. limma, DESeq2, GFOLD,
#' glm.pois, glm.qlll, and glm.nb.
#' @param paired Boolean, specifying whether perform paired comparison.
#' @param app.dir The path to application (e.g. GFOLD).
#'
#' @return An ExprDataSet instance.
#' @seealso \code{\link{ExprDataSet-class}}
#'
#' @author Wubing Zhang
#'
#' @importFrom limma eBayes topTable lmFit
#' @importFrom DESeq2 DESeqDataSetFromMatrix results DESeq design
#' @import msmsTests
#' @export

DEAnalyze <- function(obj, SampleAnn = NULL,
                      type = "Array", method = "limma",
                      paired = FALSE,
                      app.dir = "/Users/Wubing/Applications/gfold/gfold"){
  requireNamespace("edgeR")
  #### Create a new object ####
  if(is.matrix(obj) | is.data.frame(obj)){
    colnames(SampleAnn)[1] = "Condition"
    if(paired) colnames(SampleAnn)[2] = "Sibs"
    expr <- as.matrix(obj[, rownames(SampleAnn)])
    obj = new("ExprDataSet", rawdata = expr, SampleAnn = SampleAnn, type = type)
  }

  #### Build design matrix ####
  if(paired){
    Sibs = factor(slot(obj, "SampleAnn")$Sibs)
    Condition = factor(slot(obj, "SampleAnn")$Condition)
    design = model.matrix(~Sibs+Condition)
  }else{
    design = model.matrix(~1+Condition, slot(obj, "SampleAnn"))
    rownames(design) = colnames(slot(obj, "rawdata"))
  }
  idx_r = rowSums(slot(obj, "rawdata"), na.rm = TRUE)!=0
  data = slot(obj, "rawdata")[idx_r, ]
  if(tolower(type) == "array"){
    requireNamespace("limma")
    slot(obj, "normlized") = normalizeQuantiles(data)
    #"ls" for least squares or "robust" for robust regression
    fit = eBayes(lmFit(slot(obj, "normlized"), design, na.rm=TRUE))
    res = topTable(fit, adjust.method="BH", coef=ncol(design), number = Inf)
    res = res[, c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
    colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
  }else if(tolower(type) == "rnaseq"){
    if(tolower(method) == "deseq2"){
      requireNamespace("DESeq2")
      slot(obj, "normlized") = TransformCount(data, method = "vst")
      # DESeq2
      dds = DESeqDataSetFromMatrix(data, colData = slot(obj, "SampleAnn"), design = design)
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::lfcShrink(dds, coef = ncol(design), quiet = TRUE)
      res$padj[is.na(res$padj)] = 1
      res = res[, c("baseMean", "log2FoldChange", "stat", "pvalue", "padj")]
      colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    }else if(tolower(method) == "limma"){
      requireNamespace("limma")
      slot(obj, "normlized") = TransformCount(slot(obj, "rawdata"), method = "voom")
      # limma:voom
      dge <- DGEList(counts=data)
      dge <- calcNormFactors(dge)
      dge <- voom(dge, design, plot=FALSE)
      fit <- eBayes(lmFit(dge, design))
      res = topTable(fit, adjust.method="BH", coef=ncol(design), number = nrow(slot(obj, "rawdata")))
      res = res[, c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
      colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    }else if(tolower(method) == "edger"){
      requireNamespace("edgeR")
      slot(obj, "normlized") = TransformCount(slot(obj, "rawdata"), method = "voom")
      dge <- DGEList(counts=data)
      dge <- calcNormFactors(dge)
      dge <- estimateDisp(dge, design, robust=TRUE)
      fit <- glmFit(dge, design)
      lrt <- glmLRT(fit)
      res <- topTags(lrt, n = nrow(slot(obj, "rawdata")))
      res = res$table[, c("logCPM", "logFC", "logFC", "PValue", "FDR")]
      colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    }else if(tolower(method) == "gfold"){
      slot(obj, "normlized") = TransformCount(slot(obj, "rawdata"), method = "voom")
      # GFOLD
      tmp = mapply(function(x){
        write.table(cbind(NA, data[,x], NA, NA),
                    file=paste0(colnames(data)[x], ".txt"),
                    sep="\t", col.names=FALSE)}, x=1:ncol(data))
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
    }else{
      stop("Method not available for RNA-seq data !!!")
    }
  }else if(tolower(type) == "msms"){
    if (tolower(method) == "limma"){
      requireNamespace("limma")
      slot(obj, "normlized") = data
      #"ls" for least squares or "robust" for robust regression
      fit = eBayes(lmFit(slot(obj, "normlized"), design))
      res = topTable(fit, adjust.method="BH", coef=ncol(design), number = Inf)
      res = res[, c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
      colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    }else if(grepl("^glm\\.", tolower(method))){
      fd <- data.frame(gene = rownames(exprSet),
                       row.names = rownames(exprSet),
                       stringsAsFactors = FALSE)
      MSnSet_obj <- MSnSet(exprs=slot(obj, "normlized"), fData=fd,
                           pData=slot(obj, "SampleAnn"))
      MSnSet_obj <- pp.msms.data(MSnSet_obj)  # pp.msms.data function used to deleted genes which all expression is 0.

      null.f <- "y~1"
      alt.f <- "y~Condition"
      div <- colSums(slot(obj, "normlized"), na.rm = TRUE)
      ### msmsTests method
      if(tolower(method)=="glm.pois"){
        res <- msms.glm.pois(MSnSet_obj, alt.f, null.f, div=div)
      }else if(tolower(method)=="glm.qlll"){
        res <- msms.glm.qlll(MSnSet_obj, alt.f, null.f, div=div)
      }else if(tolower(method)=="glm.nb"){
        res <- msms.edgeR(MSnSet_obj, alt.f, null.f, div=div)
      }
      res$baseMean = rowMeans(slot(obj, "normlized"))[rownames(res)]
      res$padj = p.adjust(res$p.value, method = "BH")
      res = res[, c(4,1:3,5)]
      colnames(res) = c("baseMean", "log2FC", "stat", "pvalue", "padj")
    }
  }else{
    stop("Data type error! ")
  }
  slot(obj, "DEGRes") = as.data.frame(res)
  return(obj)
}

