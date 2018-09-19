#' Compute correlation of between expression of two genesets
#'
#' @docType methods
#' @name corrView
#' @rdname corrView
#'
#' @param mat Matrix like object.
#' @param geneset1 Gene set 1 which should be row names of mat.
#' @param geneset2 Gene set 1 which should be row names of mat.
#' @param xlab Character specifying the label of x axis.
#' @param ylab Character specifying the label of y axis.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned
#' and further customized.
#'
#' @author Wubing Zhang
#' @import ggplot2
#' @export

corrView <- function(mat, geneset1, geneset2, xlab = geneset1[1],
                     ylab = geneset2[1]){
  geneset1 = intersect(geneset1, rownames(mat))
  geneset2 = intersect(geneset2, rownames(mat))
  gg = data.frame(x = colMeans(mat[geneset1, , drop = FALSE], na.rm = TRUE),
                  y = colMeans(mat[geneset2, , drop = FALSE], na.rm = TRUE))
  gg = gg[is.finite(gg$x) & is.finite(gg$y), ]
  tmp = cor.test(gg$x, gg$y)
  corr = paste('Correlation = ', round(tmp$estimate,3),
               "\nP.value = ", signif(tmp$p.value,3))
  p = ggplot(gg, aes(x, y))
  p = p + geom_point(size=0.5, color="blue")
  p = p + geom_smooth(method='lm',se=FALSE, size=0.5, color = "#c12603")
  # p = p + geom_abline(slope = 1, intercept = 0, color="gray50", linetype=2)
  p = p + labs(x=xlab, y=ylab)
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99))
  if(tmp$estimate>0){
    p = p + annotate("text", x = max(gg$x), y = min(gg$y),
                     label = corr, hjust=1, vjust=0)

  }else{
    p = p + annotate("text", x = max(gg$x), y = max(gg$y),
                     label = corr, hjust=1, vjust=1)
  }
  p
}
