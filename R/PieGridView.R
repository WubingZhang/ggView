#' A matrix of pie plots
#' @param dat A data frame for plotting.
#' @param y Character, specifying the percentage of each group.
#' @param fill Character, specifying the groups.
#' @param rows Character, specifying the row in the pie matrix.
#' @param cols Character, specifying the column in the pie matrix.
#' @param percentage Boolean specifying whether convert the frequency to percentage before plotting.
#' @param size Numeric, specifying the pie size.
#' @param row.label A named character vector, specifying the row labels.
#' @param col.label A named character vector, specifying the col labels.
#' @export
PieGridView <- function(dat, y = "Freq", fill = "Note",
                        rows = "Gene", cols = "Cancer",
                        percentage = TRUE, size = 2,
                        row.label = NULL, col.label = NULL){
  requireNamespace("ggplot2")
  colnames(dat)[colnames(dat)==rows] = "rows"
  colnames(dat)[colnames(dat)==cols] = "cols"

  if(percentage){
    tmp = aggregate(dat[,y], list(row=dat$rows, col=dat$cols), sum)
    rownames(tmp) = paste0(tmp$row, "_", tmp$col)
    dat[,y] = dat[,y] / tmp[paste0(dat$rows, "_", dat$cols), 3]
  }
  dat$x = ""
  if(is.null(row.label)){
    row.label = unique(dat$rows)
    names(row.label) = row.label
  }
  if(is.null(col.label)){
    col.label = unique(dat$cols)
    names(col.label) = col.label
  }
  p = ggplot(dat, aes_string(x="x", y=y, fill=fill))
  p = p + geom_bar(stat="identity", width=1)
  p = p + theme_void(base_size = 14)
  p = p + facet_grid(rows~cols, switch="both", labeller = labeller(rows = row.label, cols = col.label))
  p = p + coord_polar("y", start=0)
  p = p + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5),
                strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.spacing = unit(-0.1*size, "lines"))
  p
}
