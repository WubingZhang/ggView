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
#' @param comparisons A list, customizing the comparisons.
#' @param test.method Character, specifying the statistical test method.
#' @param add.jitter Boolean, whether add jitter into the plot.
#' @param jitter.color Character, specifying the column name/color of the points.
#' @param jitter.size Numeric, specifying size of the jitter.
#' @param alpha A numeric, specifying the transparency of the dots.
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
                    color = NA,
                    fill = NA,
                    width = 0.6,
                    size = 1,
                    comparisons = NULL,
                    test.method = "t.test",
                    add.jitter = FALSE,
                    jitter.color = color,
                    jitter.size = size,
                    alpha = 0.6,
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
  if(color %in% colnames(gg)){ ## Customize colors in the data
    if(fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(color = color, fill = fill),
                           width = width, size = size, alpha = alpha, ...)
    else if(!boo2)
      p = p + geom_boxplot(aes_string(color = color), fill = fill,
                           width = width, size = size, alpha = alpha, ...)
    else
      p = p + geom_boxplot(aes_string(color = color), width = width,
                           size = size, alpha = alpha, ...)
  }else if(!boo1){ ## Customize single color
    if(fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(fill = fill), color = color,
                           width = width, size = size, alpha = alpha, ...)
    else if(!boo2)
      p = p + geom_boxplot(color = color, fill = fill, width = width,
                           size = size, alpha = alpha, ...)
    else
      p = p + geom_boxplot(color = color, width = width, size = size, alpha = alpha, ...)
  }else{ ## Use default color
    if(fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(fill = fill), width = width,
                           size = size, alpha = alpha, ...)
    else if(!boo2)
      p = p + geom_boxplot(fill = fill, width = width, size = size, alpha = alpha, ...)
    else
      p = p + geom_boxplot(width = width, size = size, alpha = alpha, ...)
  }
  ## Comparisons
  if(!is.null(comparisons)){
    p = p + ggpubr::compare_means(method = test.method, p.adjust.method = "fdr")
  }
  if(add.jitter){
    boo <- try(col2rgb(jitter.color), silent=TRUE)
    boo3 = "try-error" %in% class(boo)
    if(jitter.color %in% colnames(gg))
      p = p + geom_jitter(aes_string(color = jitter.color), size = jitter.size)
    else if(!boo3)
      p = p + geom_jitter(color = jitter.color, size = jitter.size)
    else
      p = p + geom_jitter(size = jitter.size)
  }
  p = p + labs(x = xlab, y = ylab, title = main)
  p = p + theme_bw(base_size = 12)
  p = p + theme(plot.title = element_text(hjust = 0.5),
                legend.position = legend.position)
  return(p)
}
