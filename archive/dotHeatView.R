#' Dot heatmap
#'
#' @docType methods
#' @name dotHeatView
#' @rdname dotHeatView
#'
#' @param mat Matrix like object, each row is gene and each column is sample.
#' @param limit Two-length vector.
#' @return An object created by \code{ggplot}, which can be assigned
#' and further customized.
#'
#' @author Wubing Zhang
#' @import reshape2
#' @import ggplot2
#' @export

dotHeatView <- function(mat, limit = c(-2,2)){
  mat = as.data.frame(mat)
  mat$Gene = rownames(mat)
  gg = reshape2::melt(mat, id="Gene")
  gg$Gene = factor(gg$Gene, levels = rownames(mat))
  gg = gg[!is.na(gg$value), ]
  gg$value[gg$value>limit[2]] = limit[2]
  gg$value[gg$value<limit[1]] = limit[1]
  gg$size = abs(gg$value)
  # gg = gg[gg$size>0.1, ]
  p = ggplot(gg, aes(x=Gene, y=variable))
  p = p + geom_point(aes(color=value, size=size))
  p = p + labs(x=NULL, y=NULL, title=NULL, size=NULL, color="Z-score")
  p = p + scale_color_gradient2(low = "#033472", high = "#c12603")
  p = p + scale_size(guide = "none")
  p = p + theme(text = element_text(colour="black",size = 14),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"),
                axis.text.x=element_text(angle = 45, hjust=1, vjust = 1))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p
}
