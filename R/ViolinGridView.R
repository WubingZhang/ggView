ViolinGridView <- function(dat, x, y, fill, rows, cols, label, label.y,
                           comparison = "Tumor",
                           alpha = 0.6, violin.size = 0.3, point.size = 0.1){
  requireNamespace("ggplot2")
  requireNamespace("Hmisc")
  colnames(dat)[colnames(dat)==rows] = "rows"
  colnames(dat)[colnames(dat)==cols] = "cols"

  #### Comparisons
  if(!is.null(comparison)){
    ## wilcox test of tumor vs other samples
    wilcox.p = apply(unique(dat[, c("rows", "cols")]), 1, function(x){
      wilcox.test(dat[dat$rows==x[1]&dat$cols==x[2]&dat[,fill]==comparison, y],
                  dat[dat$rows==x[1]&dat$cols==x[2]&dat[,fill]!=comparison, y])$p.value
    })
    wilcoxRes = unique(dat[, c("rows", "cols")])
    wilcoxRes$Pval = wilcox.p
    wilcoxRes$Sign = ""
    wilcoxRes$Sign[wilcoxRes$Pval<0.05] = "*"
    wilcoxRes$Sign[wilcoxRes$Pval<0.01] = "**"
    wilcoxRes$Sign[wilcoxRes$Pval<0.001] = "***"
    dat = merge(dat, wilcoxRes, by = c("rows", "cols"), all = TRUE)
    label.y = aggregate(dat[,y], list(rows = dat$rows), max, na.rm = TRUE)
    colnames(label.y) = c("rows", "label.y")
    dat = merge(dat, label.y, by = "rows", all = TRUE)
    label = "Sign"
    label.y = "label.y"
  }

  p = ggplot(dat, aes_string(x=x, y=y, fill=fill))
  p = p + geom_violin(trim = TRUE, alpha = alpha, size = violin.size)
  p = p + stat_summary(fun.data = mean_sdl, geom="pointrange", size = point.size)
  if(!is.null(label)){
    textDat = unique(dat[, c(x, label.y, label, "rows", "cols")])
    idx = duplicated(textDat[, c("rows", "cols")])
    textDat = textDat[!idx, ]
    p = p + geom_text(aes_string(x = x, y = label.y, label = label), data = textDat)
  }
  p = p + facet_grid(rows~cols, switch="both", scales = "free",
                     labeller = labeller(rows = row.label, cols = col.label))
  p = p + labs(x = NULL, y = NULL, color = NULL, fill = NULL)
  p = p + theme_bw()
  p = p + theme(axis.text = element_blank(), axis.ticks = element_blank())
  p = p + theme(panel.grid = element_blank())
  p = p + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5),
                strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                strip.background = element_blank(),
                panel.spacing = unit(-0.1*size, "lines"))
  p
}
