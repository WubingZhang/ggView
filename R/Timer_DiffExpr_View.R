#' Visualize the differential expression results in TCGA data (Timer)
#'
#' @docType methods
#' @name Timer_DiffExpr_View
#' @rdname Timer_DiffExpr_View
#'
#' @param gene A vector of gene names.
#' @param exprdir The directory of processed TCGA expression profiles.
#' @examples
#' args <- commandArgs(T)
#' if(length(args)==0 | args[1]=="-h"){
#'   message("Usage: Rscript Timer_DiffExpr_View.R gene outdir")
#'   q("no")
#' }
#' gene = args[1]
#' outdir = args[2]
#' @import ggplot2 dplyr ggpubr
#' @export
#'
Timer_DiffExpr_View <- function(gene, exprdir = "ProcessedDat/"){
  require(ggplot2)
  require(dplyr)
  require(ggpubr)
  CancerType = list.files(path = exprdir)
  message(Sys.time(), " # Read TPM ...")
  gg = data.frame(CancerType = c(), TPM = c(), stringsAsFactors = FALSE)
  for(ct in CancerType){
    tpm = readRDS(paste0(exprdir, "/", ct, "/", ct, ".RNAseq.TPM.rds"))
    gg = rbind.data.frame(gg, data.frame(CancerType = ct, TPM = tpm[gene, ],
                                         stringsAsFactors = FALSE))
  }
  gg = gg[!gg$CancerType %in% c("COADREAD", "GBMLGG"), ]
  gg$Cluster = "Normal"
  gg$Cluster[grepl("01$", rownames(gg))] = "Tumor"
  gg$Cluster[grepl("06$", rownames(gg))] = "Metastasis"

  ## Filter out subset of cancer types ##
  idx1 = (gg$Cluster=="Metastasis") & (gg$CancerType!="SKCM")
  idx2 = (gg$Cluster=="Normal") & (gg$CancerType=="SKCM")
  gg = gg[!(idx1|idx2), ]
  gg$Cluster = factor(gg$Cluster, levels = c("Normal", "Metastasis", "Tumor"))
  tmp = sapply(unique(gg$CancerType), function(x) mean(gg$TPM[gg$CancerType==x]))
  gg$CancerType = factor(gg$CancerType, levels = names(sort(tmp)))
  gg = gg[order(gg$CancerType, gg$Cluster, gg$TPM), ]
  gg$x = paste0(gg$CancerType, ".", gg$Cluster)
  gg$x = factor(gg$x, levels = unique(gg$x))

  ## Remove part of cancer types ##
  tmp = sapply(levels(gg$CancerType), function(x){
    tpm1 = median(gg$TPM[gg$CancerType==x & gg$Cluster=="Normal"])
    tpm2 = median(gg$TPM[gg$CancerType==x & gg$Cluster=="Tumor"])
    tpm3 = median(gg$TPM[gg$CancerType==x & gg$Cluster=="Metastasis"])
    ifelse(!is.na(tpm1), tpm1>tpm2, tpm3<tpm2)
  })
  gg = gg[!gg$CancerType%in%levels(gg$CancerType)[is.na(tmp)|tmp], ]


  ## Prepare dataframe for ggplot ##
  message(Sys.time(), " # Remove cancer types without normal control ...")
  tmp = table(gsub("\\..*", "", unique(gg$x)))
  gg2 = gg[gg$CancerType %in% names(tmp)[tmp==2], ]
  gg2$x = factor(as.character(gg2$x), levels = unique(as.character(gg2$x)))
  comparisons = lapply(unique(gg2$CancerType), function(x)
    as.character(unique(gg2$x[gg2$CancerType==x])))

  ## Prepare data for geom_rect ##
  rect <- data.frame(x=c(grep("\\.[NM]", levels(gg2$x)),
                         grep("HPVneg", levels(gg2$x)))) %>%
    mutate(xmin=x-0.45, xmax=x+1.45)

  ## Plot figure ##
  message(Sys.time(), " # Draw the figure ...")
  p = ggplot(gg2, aes(x = x, y = TPM))
  p = p + geom_rect(data=rect, aes(NULL,NULL, xmin=xmin, xmax=xmax),
                    ymin=-Inf, ymax=Inf, fill="Gray90")
  p = p + geom_boxplot(aes(color = Cluster), outlier.shape = NA)
  p = p + stat_compare_means(comparisons = comparisons, method="t.test",
                             label = "p.signif", hide.ns = TRUE,
                             label.y = quantile(gg2$TPM, 0.999),
                             tip.length = 0, bracket.size = 0)
  p = p + geom_jitter(aes(color = Cluster), size = 0.6, shape = 1, alpha=0.4)
  p = p + labs(x = NULL, y = "Log2TPM")
  p = p + scale_color_manual(values = c("#292CFB", "#9C31E9", "#FB1618"), guide = "none")
  p = p + ylim(0, quantile(gg2$TPM, 0.9995)+0.2)
  p = p + theme(legend.key = element_rect(fill = "transparent", color = NA))
  p = p + theme(text = element_text(colour="black",size = 10, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=12),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(),
                panel.background = element_blank())
  p = p + theme(legend.title=element_blank())
  p
}

