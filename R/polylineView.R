#' polyline plot
#'
#' Polyline plot for each column
#'
#' @docType methods
#' @name polylineView
#' @rdname polylineView
#'
#' @param gg Data frame
#' @param x Can be integer (variable position) or string (variable name).
#' @param xlab Label of x-axis
#' @param ylab Label of y-axis
#' @param colab Label of color
#' @param main As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @export

polylineView <- function(gg, x = 0, xlab = NULL, ylab = "Count",
                         colab = NULL, main=NULL,
                         filename=NULL, width=5, height =4, ...){
  requireNamespace("ggplot2")
  # gg = gg[, sort(colnames(gg))]
  if(x==0){
    gg$x = 1:nrow(gg)
  }else{
    idx = ifelse(is.integer(x), x, which(colnames(gg)==x))
    colnames(gg)[idx] = "x"
  }
  gg = reshape2::melt(gg, id = "x")
  p = ggplot(gg, aes(x = x, y=value, color = variable))
  p = p + geom_point()
  p = p + geom_line()
  p = p + labs(x = xlab, y = ylab, color = colab, title = main)
  p = p + scale_x_continuous(breaks = unique(gg$x))
  p = p + scale_y_continuous(expand = c(0,0))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_blank())

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
