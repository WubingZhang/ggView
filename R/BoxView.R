#' Boxplot
#'
#' @docType methods
#' @name BoxView
#' @rdname BoxView
#'
#' @param gg A data frame.
#' @param x Character, specifying the column name for x plotting.
#' @param y Character, specifying the column name for y plotting.
#' @param color Character, specifying the column name/color.
#' @param fill Character, specifying the column name/fill color.
#' @param width Numeric, specifying the box width.
#' @param size Numeric, specifying size of the box.
#'
#' @param comparisons A list, customizing the comparisons.
#' @param test.method Character, specifying the statistical test method.
#' @param p.label Character, "p.signif" (shows the significance levels), or
#' "p.format" (shows the formatted p value).
#'
#' @param add.jitter Boolean, whether add jitter into the plot.
#' @param jitter.color Character, specifying the column name/color of the points.
#' @param jitter.size Numeric, specifying size of the jitter.
#' @param jitter.shape Numeric, specifying the shape of the dots.
#' @param jitter.width Numeric, specifying the width of the jitter region.
#'
#' @param label Character, specifying the column name for dot label.
#' @param toplabels Character, specifying the column name/color of the points.
#' @param label.size Numeric, specifying size of the labels.
#' @param label.force Numeric, force of repulsion between overlapping text labels. Defaults to 3.
#'
#' @param alpha A numeric, specifying the transparency of the boxes.
#' @param main Title of the figure.
#' @param xlab Title of x-axis
#' @param ylab Title of y-axis.
#' @param legend.position Position of legend, "none", "right", "top", "bottom", or
#' a two-length vector indicating the position.
#' @param ... Other parameters in geom_boxplot.
#'
#' @return An object created by \code{ggplot}, which can be assigned
#' and further customized.
#'
#' @author Wubing Zhang
#' @importFrom ggpubr compare_means
#' @export

BoxView <- function(gg, x, y,
                    color = "black",
                    fill = NA,
                    width = 0.5,
                    size = 1,
                    comparisons = NULL,
                    test.method = "wilcox.test",
                    p.label = c("p.signif", "p.format")[1],
                    add.jitter = FALSE,
                    jitter.color = color,
                    jitter.size = size,
                    jitter.shape = 16,
                    jitter.width = 0.2,
                    jitter.height = 0,
                    label = NULL,
                    toplabels = c(),
                    label.size = 5,
                    label.force = 3,
                    alpha = 1,
                    main = NULL,
                    xlab = x,
                    ylab = y,
                    legend.position = "none",
                    ...){
  gg = as.data.frame(gg)
  p = ggplot(gg, aes_string(x, y))

  ## Check if color is valid color
  boo <- try(col2rgb(color), silent=TRUE)
  boo1 = "try-error" %in% class(boo)
  boo <- try(col2rgb(fill), silent=TRUE)
  boo2 = "try-error" %in% class(boo)
  if(is.na(fill)) boo2 = TRUE

  if(color%in%colnames(gg)){ ## Customize colors in the data
    if(fill%in%colnames(gg))
      p = p + geom_boxplot(aes_string(color = color, fill = fill),
                           width = width, size = size, alpha = alpha,
                           outlier.shape = NA, ...)
    else if(!boo2)
      p = p + geom_boxplot(aes_string(color = color), fill = fill,
                           width = width, size = size, alpha = alpha,
                           outlier.shape = NA, ...)
    else
      p = p + geom_boxplot(aes_string(color = color),
                           width = width, size = size, alpha = alpha,
                           outlier.shape = NA, ...)
  }else if(!boo1){
    if(fill%in%colnames(gg))
      p = p + geom_boxplot(aes_string(fill = fill), color = color,
                           width = width, size = size, alpha = alpha,
                           outlier.shape = NA, ...)
    else if(!boo2)
      p = p + geom_boxplot(color = color, fill = fill,
                           width = width, size = size, alpha = alpha,
                           outlier.shape = NA, ...)
    else
      p = p + geom_boxplot(color = color, width = width, size = size,
                           alpha = alpha, outlier.shape = NA, ...)
  }else{
    if(fill%in%colnames(gg))
      p = p + geom_boxplot(aes_string(fill = fill), width = width, size = size,
                           alpha = alpha, outlier.shape = NA, ...)
    else if(!boo2)
      p = p + geom_boxplot(fill = fill, width = width, size = size, alpha = alpha,
                           outlier.shape = NA, ...)
    else
      p = p + geom_boxplot(width = width, size = size, alpha = alpha, outlier.shape = NA, ...)
  }
  ## Comparisons
  if(!is.null(comparisons)){
    p = p + ggpubr::stat_compare_means(comparisons = comparisons,
                                       method = test.method,
                                       label = p.label)
  }

  ## Add dots in the plot
  if(add.jitter){
    boo <- try(col2rgb(jitter.color), silent=TRUE)
    boo3 = "try-error" %in% class(boo)
    if(jitter.color %in% colnames(gg))
      p = p + geom_jitter(aes_string(color = jitter.color),
                          width = jitter.width, height = jitter.height,
                          shape=jitter.shape, size = jitter.size)
    else if(!boo3)
      p = p + geom_jitter(color = jitter.color,
                          width = jitter.width, height = jitter.height,
                          shape=jitter.shape, size = jitter.size)
    else
      p = p + geom_jitter(width = jitter.width, height = jitter.height,
                          shape=jitter.shape, size = jitter.size)
  }
  if(!is.null(label) && label%in%colnames(gg)){
    colnames(gg)[colnames(gg)==label] = "label"
    gg$label[!gg$label%in%toplabels] = ""
    p = p + ggrepel::geom_text_repel(data = gg,
                                     aes(label = label),
                                     size = label.size,
                                     force = label.force)
  }

  p = p + labs(x = xlab, y = ylab, title = main, color = NULL, fill = NULL)
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5),
                legend.position = legend.position)
  return(p)
}
