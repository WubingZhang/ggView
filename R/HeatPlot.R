#' Typical heatmap from Jingxin
#'
#' @param dat Data frame.
#' @param x The variable of x.
#' @param y The variable of y.
#' @param rho The variable of color.
#' @param pvalue The pvalue column.
#' @param rho.limit A numberic vector with three values.
#' @param bar_anno Limit on lengend.
#' @param box_size A numberic value to specifying the box size.
#' @param barwidth A numberic value to specifying the width of colorbar.
#' @param barheight A numberic value to specifying the height of colorbar.
#'
#' @return ggplot object
#' @author Jingxin Fu
#'
HeatmapJView <- function(dat, x, y, rho, pvalue,
                         rho.limit=c(-1,0,1), bar_anno=c(-1,0,1),
                         box_size=9, barwidth = 0.8, barheight = 5){
  require(ggplot2)
  min_value = limit[1]
  mean_value = limit[2]
  max_value = limit[3]
  dat$sig = NA
  dat$sig[dat[,pvalue] > 0.1] ='FDR > 0.1'
  dat$sig[dat[,pvalue] < 0.1] = 'FDR < 0.1'
  dat$sig = factor(dat$sig, levels = c('FDR > 0.1','FDR < 0.1'))
  dat = na.omit(dat)

  ggplot(dat, aes_string(x, y)) +
    geom_point(aes_string(color = rho,shape = 'sig'),size=box_size) +
    scale_shape_manual(name=NULL, values = c('FDR > 0.1' = 7,'FDR < 0.1' = 15), drop = FALSE)+
    coord_equal()+labs(x=NULL,y=NULL) +
    scale_color_gradientn(colours = c("navy","white","firebrick3"),
                          values = rescale(c(min_value,mean_value,max_value)),
                          guide = "colorbar", limits=c(min_value,max_value),
                          breaks=bar_anno,na.value="white")+
    scale_fill_gradientn(colours = c("navy","white","firebrick3"),
                         values = rescale(c(min_value,mean_value,max_value)),
                         guide = "colorbar", limits=c(min_value,max_value),
                         breaks=bar_anno,na.value="white")+
    mytheme+guides(fill = guide_colorbar(barwidth = barwidth, barheight = barheight),
                   color = guide_colorbar(barwidth = barwidth, barheight = barheight))
}
