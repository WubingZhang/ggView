#' Boxplot
#'
#' Boxplot for each column
#'
#' @docType methods
#' @name BoxView
#' @rdname BoxView
#'
#' @param gg Data frame
#' @param ctrl Index of control samples
#' @param treat Index of treat samples
#' @param label.c Label of control box
#' @param label.t Label of treatment box
#' @param ylab Label of y-axis
#' @param xlab Label of x-axis
#' @param colab Label of color
#' @param flip Whether flip coordinate
#' @param size Size of box and jitter
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
#' @importFrom reshape melt
#' @import ggpubr
#'
#' @export

BoxView <- function(gg, ctrl, treat, label.c = "Control", label.t = "Treatment", ylab="Expression", xlab=NULL,
                    colab=NULL, flip=FALSE, size = 1, main=NULL, filename=NULL, width=5, height =4, ...){
  gg = as.data.frame(gg)
  gg$Label = NA
  gg$Label[ctrl] = label.c
  gg$Label[treat] = label.t
  gg$Label = factor(gg$Label, levels = c(label.c, label.t))
  gg = reshape::melt(gg, id="Label")
  gg = gg[!is.na(gg$value), ]
  #=====Boxplot and compare paired samples====
  p = ggboxplot(gg, x="variable", y="value", color="Label", size=size,
                font.label=list(size=6), outlier.size=0)
  #Add p-values and significance levels to ggplots
  p = p + geom_jitter(aes(color=Label), size=size)
  p = p + stat_compare_means(aes(group = Label), label = "p.signif")
  p = p + labs(x=xlab, y=ylab, color=colab)
  if(flip) p = p + coord_flip()
  p = p + scale_color_manual(values = c("#365ebc", "#e41a1c"))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.position="right", legend.text = element_text(size=8))
  # p = p + ggrepel::geom_text_repel(aes(label=Gene, color=color))
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", dpi=600, width=width, height=height, ...)
  }
  return(p)
}
