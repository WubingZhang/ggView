#' K-means clustering and plot Heatmap
#'
#' @docType methods
#' @name KcHeatView
#' @rdname KcHeatView
#'
#' @param dd Matrix like object.
#' @param cluster_row Boolean, whether cluster rows.
#' @param cluster_col Boolean, whether cluster columns.
#' @param k Integer, how many clusters to be clustered.
#' @param filename Character, specifying output pdf name.
#' @param height Integer, sepcifying height of heatmap.
#' @param width Integer, sepcifying width of heatmap.
#' @param ... Other available parameters for pheatmap
#'
#' @return Invisibly a pheatmap object.
#'
#' @author Wubing Zhang
#'
#' @import pheatmap
#' @import RColorBrewer
#' @export

KcHeatView <- function(dd, limit=2, cluster_row=TRUE, cluster_col=TRUE, k=5,
                  filename=NA, height=0.12*nrow(dd)+3, width=0.12*ncol(dd)+7, ...){
  dd[is.na(dd)] = 0
  row_ann = NA; col_ann = NA
  cutrow = NA; cutcol = NA
  dd[dd>limit] = limit
  dd[dd< -limit] = -limit
  if(cluster_row){
    cluster_genes = kmeans(dd, centers = k, iter.max=100)
    clusters_G = sort(cluster_genes$cluster)
    row_ann = data.frame(Cluster=factor(paste0("Cluster", clusters_G),
                                        levels = unique(paste0("Cluster", clusters_G))))
    rownames(row_ann) = names(clusters_G)
    dd = dd[names(clusters_G), ]
    cutrow = k
  }

  if(cluster_col){
    cluster_samples = kmeans(t(dd), centers = k, iter.max=100)
    clusters_S = sort(cluster_samples$cluster)
    col_ann = data.frame(Cluster=factor(paste0("Cluster", clusters_S),
                                        levels=unique(paste0("Cluster", clusters_S))))
    rownames(col_ann) = names(clusters_S)
    dd = dd[, names(clusters_S)]
    cutcol = k
  }
  #--Heatmap---
  breaks = seq(-limit, limit, length.out = 200)
  color = rev(colorRampPalette(c("#c12603", "white", "#033472"), space = "Lab")(199))
  # breaks = seq(-max(abs(dd), na.rm = TRUE), max(abs(dd), na.rm = TRUE), length.out = 200)
  # color = rev(colorRampPalette(RColorBrewer::brewer.pal(5, "RdBu"))(199))
  if(cluster_row&cluster_col)
    pheatmap::pheatmap(dd, col=color, annotation_row = row_ann, annotation_col=col_ann,
             breaks = breaks, cutree_rows = cutrow, cutree_cols = cutcol,
             filename = filename, height = height, width = width,
             cluster_rows = FALSE, cluster_cols = FALSE, border_color=NA, scale="none", ...)
  else if(cluster_row)
    pheatmap::pheatmap(dd, col=color, annotation_row = row_ann,
                       breaks = breaks, cutree_rows = cutrow, cutree_cols = cutcol,
                       filename = filename, height = height, width = width,
                       cluster_rows = FALSE, cluster_cols = FALSE, border_color=NA, scale="none", ...)
  else if(cluster_col)
    pheatmap::pheatmap(dd, col=color, annotation_col = col_ann,
                       breaks = breaks, cutree_rows = cutrow, cutree_cols = cutcol,
                       filename = filename, height = height, width = width,
                       cluster_rows = FALSE, cluster_cols = FALSE, border_color=NA, scale="none", ...)
  else
    pheatmap::pheatmap(dd, col=color, breaks = breaks, cutree_rows = cutrow, cutree_cols = cutcol,
                       filename = filename, height = height, width = width,
                       cluster_rows = FALSE, cluster_cols = FALSE, border_color=NA, scale="none", ...)
}
