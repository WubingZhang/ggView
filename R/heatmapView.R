#' Draw heatmap
#'
#' @docType methods
#' @name heatmapView
#' @rdname heatmapView
#'
#' @param mat Matrix like object, each row is gene and each column is sample.
#' @param limit Max value in heatmap
#' @param filename File path where to save the picture.
#' @param width Manual option for determining the output file width in inches.
#' @param height Manual option for determining the output file height in inches.
#'
#' @return Invisibly a pheatmap object that is a list with components.
#'
#' @author Wubing Zhang
#' @import pheatmap
#' @export

heatmapView <- function(mat, limit=2, cluster_rows=FALSE, cluster_cols=FALSE,
                        filename = NA, width = NA, height = NA, ...){
  mat[is.na(mat)] = 0
  mat[mat>limit] = limit
  mat[mat< -limit] = -limit
  breaks = seq(-limit, limit, length.out = 200)
  color = rev(colorRampPalette(c("#c12603", "white", "#033472"), space = "Lab")(199))
  pheatmap::pheatmap(mat, color=color, breaks=breaks, border_color=NA,
                     cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                     fontfamily = "Helvetica", filename=filename, width=width, height = height,...)
}
