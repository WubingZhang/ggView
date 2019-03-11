Scatter2DView <- function(data, x, y, color, label, label.list){
  requireNamespace("ggplot2")
  p = ggplot(data, aes(x, y, color = color))
  p = p + geom_point()
  p = p + scale_color_gradient2(low = "blue", mid = "yellow",
                                high = "red", midpoint = 2)
  p = p + labs(x = x, y = y, color = color)
  p = p + theme(text = element_text(colour="black", size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_blank())
  p
}
