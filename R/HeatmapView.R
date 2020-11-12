#' Draw heatmap
#'
#' @docType methods
#' @name HeatmapView
#' @rdname HeatmapView
#'
#' @param mat Matrix like object
#' @param breaks A vector indicating numeric breaks
#' @param colors A vector of colors which correspond to values in breaks
#' @param na_col Color for NA values.
#'
#' @param cluster_rows Same as that in ComplexHeatmap::Heatmap.
#' @param row_dend_side Same as that in ComplexHeatmap::Heatmap.
#' @param cluster_cols Same as that in ComplexHeatmap::Heatmap.
#' @param column_dend_side Same as that in ComplexHeatmap::Heatmap.
#'
#' @param show_row_names Same as that in ComplexHeatmap::Heatmap.
#' @param row_names_side Same as that in ComplexHeatmap::Heatmap.
#' @param row_names_gp Same as that in ComplexHeatmap::Heatmap.
#' @param row_names_rot Same as that in ComplexHeatmap::Heatmap.
#' @param show_column_names Same as that in ComplexHeatmap::Heatmap.
#' @param column_names_side Same as that in ComplexHeatmap::Heatmap.
#' @param column_names_gp Same as that in ComplexHeatmap::Heatmap.
#' @param column_names_rot Same as that in ComplexHeatmap::Heatmap.
#'
#' @param top_ann A data frame. Each column will be treated as a simple annotation.
#' The data frame must have column names.
#' @param top_ann_col A list of colors which contain color mapping to df.
#' @param bott_ann Same as top_ann.
#' @param bott_ann_col A list of colors which contain color mapping to df.
#' @param left_ann Same as top_ann.
#' @param left_ann_col A list of colors which contain color mapping to df.
#' @param right_ann Same as top_ann.
#' @param right_ann_col A list of colors which contain color mapping to df.
#' @param show_ann_name Whether show annotation names.
#' @param show_legend Whether show annotation legends.
#'
#' @param row_split A vector or a data frame by which the rows are split.
#' But if cluster_rows is a clustering object, split can be a single number
#' indicating to split the dendrogram by cutree.
#' @param column_split Same as row_split.
#'
#' @param show_heatmap_legend Whether show legends.
#' @param legend_title Character specifyin the legend title.
#' @param legend_title_position Position of title relative to the legend.
#' topleft, topcenter, leftcenter-rot and lefttop-rot are only for vertical legend and
#' leftcenter, lefttop are only for horizontal legend.
#' @param legend_direction Vertical or horizontal?
#' @param legend_title_gp Graphic parameters of the title.
#' @param legend_labels_gp Graphic parameters for labels.
#' @param legend_height Height of the whole legend body. It is only used for vertical continous legend.
#' @param legend_width Width of the whole legend body. It is only used for horizontal continous legend.
#' @param legend_side Side to put heatmap legend
#' @param ... Other parameters in draw.
#'
#' @return No return.
#' @author Wubing Zhang
#'
#' @examples
#' dat = matrix(rnorm(100), 10)
#' rownames(dat) = letters[1:10]
#' colnames(dat) = letters[11:20]
#' rowann = data.frame(Group = rep(letters[1:2], each=5), index = 1:10)
#' colann = data.frame(Group = rep(letters[1:2], each=5), index = 11:20)
#' HeatmapView(dat, left_ann = rowann, top_ann = colann)
#'
#' @importFrom ComplexHeatmap columnAnnotation rowAnnotation Heatmap draw
#' @importFrom circlize colorRamp2
#' @export

HeatmapView <- function(mat,
                        breaks=c(-2, 0, 2),
                        colors = c("blue", "#EEEEEE", "red"),
                        na_col = "grey",
                        cluster_rows = FALSE,
                        row_dend_side = c("left", "right"),
                        cluster_cols = FALSE,
                        column_dend_side = c("top", "bottom"),
                        show_row_names = TRUE,
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 12),
                        row_names_rot = 0,
                        show_column_names = TRUE,
                        column_names_side = "bottom",
                        column_names_gp = gpar(fontsize = 12),
                        column_names_rot = 90,

                        top_ann = NULL,
                        top_ann_col = NULL,
                        bott_ann = NULL,
                        bott_ann_col = NULL,
                        left_ann = NULL,
                        left_ann_col = NULL,
                        right_ann = NULL,
                        right_ann_col = NULL,
                        show_ann_name = FALSE,
                        show_legend = TRUE,

                        row_split = NULL,
                        column_split = NULL,

                        show_heatmap_legend = TRUE,
                        legend_title = NULL,
                        legend_title_position = "lefttop",
                        legend_direction = "vertical",
                        legend_title_gp = gpar(fontsize = 12),
                        legend_labels_gp = gpar(fontsize = 12),
                        legend_height = 2,
                        legend_width = 0.3,
                        legend_side = "right",
                        ...){
  mat[mat>max(breaks)] = max(breaks)
  mat[mat<min(breaks)] = min(breaks)
  colPal = circlize::colorRamp2(breaks, colors)

  if(!is.null(top_ann)){
    if(is.null(top_ann_col))
      top_ann = ComplexHeatmap::columnAnnotation(df = top_ann,
                                 show_legend = show_legend,
                                 show_annotation_name = show_ann_name)
    else
      top_ann = ComplexHeatmap::columnAnnotation(df = top_ann, col = top_ann_col,
                                                 show_legend = show_legend,
                                                 show_annotation_name = show_ann_name)
  }
  if(!is.null(bott_ann)){
    if(is.null(bott_ann_col))
      bott_ann = ComplexHeatmap::columnAnnotation(df = bott_ann,
                                                  show_legend = show_legend,
                                                  show_annotation_name = show_ann_name)
    else
      bott_ann = ComplexHeatmap::columnAnnotation(df = bott_ann, col = bott_ann_col,
                                                  show_legend = show_legend,
                                                  show_annotation_name = show_ann_name)
  }
  if(!is.null(left_ann)){
    if(is.null(left_ann_col))
      left_ann = ComplexHeatmap::rowAnnotation(df = left_ann,
                                               show_legend = show_legend,
                                               show_annotation_name = show_ann_name)
    else
      left_ann = ComplexHeatmap::rowAnnotation(df = left_ann, col = left_ann_col,
                                               show_legend = show_legend,
                                               show_annotation_name = show_ann_name)
  }
  if(!is.null(right_ann)){
    if(is.null(right_ann_col))
      right_ann = ComplexHeatmap::rowAnnotation(df = right_ann,
                                                show_legend = show_legend,
                                                show_annotation_name = show_ann_name)
    else
      right_ann = ComplexHeatmap::rowAnnotation(df = right_ann, col = right_ann_col,
                                                show_legend = show_legend,
                                                show_annotation_name = show_ann_name)
  }

  p = ComplexHeatmap::Heatmap(mat,
                              col = colPal,
                              na_col = na_col,
                              cluster_rows = cluster_rows,
                              row_dend_side = row_dend_side,
                              cluster_columns = cluster_cols,
                              column_dend_side = column_dend_side,
                              show_row_names = show_row_names,
                              row_names_side = row_names_side,
                              row_names_gp = row_names_gp,
                              row_names_rot = row_names_rot,
                              show_column_names = show_column_names,
                              column_names_side = column_names_side,
                              column_names_rot = column_names_rot,
                              column_names_gp = column_names_gp,

                              left_annotation = left_ann,
                              right_annotation = right_ann,
                              top_annotation = top_ann,
                              bottom_annotation = bott_ann,

                              row_split = row_split,
                              column_split = column_split,
                              show_heatmap_legend = show_heatmap_legend,
                              heatmap_legend_param = list(
                                title = legend_title,
                                title_gp = legend_title_gp,
                                labels_gp = legend_labels_gp,
                                legend_height = unit(legend_height, "in"),
                                legend_width = unit(legend_width, "in"),
                                title_position = legend_title_position,
                                legend_direction = legend_direction))
  ComplexHeatmap::draw(p, padding = unit(c(5, 5, 5, 5), "mm"),
                       heatmap_legend_side = legend_side,
                       annotation_legend_side = legend_side, ...)
}
