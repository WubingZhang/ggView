#' Use the customized theme
#' @param p A ggplot object
#' @param base The base size of the text in the figure
#' @return A ggplot object
#' @import ggplot
#' @export
theme_pubr <- function(p, base = 14){
  require(ggplot2)
  p = p + theme(text = element_text(colour="black",size = base, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size = base+4),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.key = element_blank(),
                legend.text = element_text(colour="black",size = base-6, family = "Helvetica"),
                legend.background = element_rect(color = NA, fill = NA, size = NA))
  p = p + theme(legend.title=element_blank())
  p
}
