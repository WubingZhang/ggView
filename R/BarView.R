#' Bar plot
#'
#' Bar plot
#'
#' @docType methods
#' @name BarView
#' @rdname BarView
#'
#' @param df A data frame.
#' @param x A character, specifying the x-axis.
#' @param y A character, specifying the x-axis.
#' @param color A character, specifying the outline color.
#' @param fill A character, specifying the fill color.
#' @param bar.width A numeric, specifying the width of bar.
#' @param position "dodge" (default), "stack", "fill".
#' @param dodge.width A numeric, set the width in position_dodge.
#' @param main A charater, specifying the figure title.
#' @param xlab A character, specifying the title of x-axis.
#' @param ylab A character, specifying the title of y-axis.
#' @param ... Other parameters in geom_bar
#'
#' @author Wubing Zhang
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#'
#' @examples
#' mdata = data.frame(group=letters[1:5], count=sample(1:100,5))
#' BarView(mdata, x = "group", y = "count")
#' @import ggplot2 ggpubr
#' @export

BarView <- function(df, x = "x", y = "y", fill = "#FC6665", color = fill,
                    bar.width = 0.8, position = "dodge",
                    dodge.width = 0.8, main = NA,
                    xlab = NULL, ylab = NA, ...){
  requireNamespace("ggplot2")
  requireNamespace("ggpubr")
  ## Check if fill is valid color
  boo <- try(col2rgb(fill), silent=TRUE)
  boo = "try-error" %in% class(boo)
  boo2 <- try(col2rgb(color), silent=TRUE)
  boo2 = "try-error" %in% class(boo2)

  ## Use the order of x in the df
  df[,x] = factor(df[,x], levels = unique(df[,x]))

  if(boo){
    p <- ggplot(df, aes_string(x, y, fill=fill))
    if(boo2) p <- ggplot(df, aes_string(x, y, fill=fill, color = color))
    else p <- ggplot(df, aes_string(x, y, fill=fill), color = color)
    if(position=="dodge"){
      p <- p + geom_bar(width = bar.width, stat="identity",
                        position=position_dodge(width = dodge.width,
                                                preserve = "single"), ...)
    }else{
      p <- p + geom_bar(width = bar.width, stat="identity", position=position, ...)
    }
  }else{
    p <- ggplot(df, aes_string(x, y))
    if(boo2) p <- ggplot(df, aes_string(x, y, color = color), fill=fill)
    else p <- ggplot(df, aes_string(x, y), fill=fill, color = color)
    if(position=="dodge"){
      p <- p + geom_bar(width = bar.width, stat="identity",
                        position=position_dodge(width = dodge.width,
                                                preserve = "single"), ...)
    }else{
      p <- p + geom_bar(width = bar.width, stat="identity", position=position, ...)
    }
  }
  p = p + scale_y_continuous(expand = c(0,0))
  p = p + labs(fill = NULL)
  if(!(length(xlab)==1 && is.na(xlab))) p = p + labs(x=xlab)
  if(!(length(ylab)==1 && is.na(ylab))) p = p + labs(y=ylab)
  if(!(length(main)==1 && is.na(main))) p = p + labs(title=main)
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  return(p)
}
