#' Scatter plot
#'
#' Scatter plot supporting groups.
#'
#' @docType methods
#' @name ScatterView
#' @rdname ScatterView
#' @aliases ScatterView
#'
#' @param data Data frame.
#' @param x A character, specifying the x-axis.
#' @param y A character, specifying the y-axis.
#'
#' @param label An integer or a character specifying the column used as the label, default value is 0 (row names).
#' @param label.top Boolean, indicates whether label the name of top hits in the groups.
#' @param top Integer, specifying the number of top terms in the groups to be labeled.
#' @param toplabels Character vector, specifying terms to be labeled.
#'
#' @param model One of "none" (default), "ninesquare", "volcano", and "rank".
#' @param groups Specify the colored groups. Optional groups include "topleft", "topcenter",
#' "topright", "midleft", "midcenter", "midright", "bottomleft", "bottomcenter", "bottomright".
#' @param group_col A vector of colors for specified groups.
#' @param groupnames A vector of group names to show on the legend.
#'
#' @param auto_cut Boolean, take 1.5 fold standard deviation as cutoff.
#' @param auto_cut_x Boolean, take 1.5 fold standard deviation as cutoff on x-axis.
#' @param auto_cut_y Boolean, take 1.5 fold standard deviation as cutoff on y-axis.
#' @param auto_cut_diag Boolean, take 1.5 fold standard deviation as cutoff on diagonal.
#' @param x_cut An one or two-length numeric vector, specifying the cutoff used for x-axis.
#' @param y_cut An one or two-length numeric vector, specifying the cutoff used for y-axis.
#' @param slope A numberic value indicating slope of the diagonal cutoff.
#' @param intercept A numberic value indicating intercept of the diagonal cutoff.
#'
#' @param display_cut Boolean, indicating whether display the dashed line of cutoffs.
#' @param legend Whether show the color legend.
#'
#' @param main Title of the figure.
#' @param xlab Title of x-axis
#' @param ylab Title of y-axis.
#' @param ... Other available parameters in function 'geom_text_repel'.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{ScatterView}}
#'
#' @examples
#'
#'
#' @importFrom ggpubr theme_pubr
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#'

ScatterView<-function(data, x = "x", y = "y", label = 0,
                      label.top = TRUE, top = 0, toplabels = NULL,
                      model = c("none", "ninesquare", "volcano", "rank")[1],
                      groups = NULL, group_col = NULL, groupnames = NULL,
                      auto_cut = TRUE, auto_cut_x = auto_cut,
                      auto_cut_y = auto_cut, auto_cut_diag = auto_cut,
                      x_cut = NULL, y_cut = NULL, slope = 1, intercept = NULL,
                      display_cut = TRUE, legend = FALSE,
                      main = NULL, xlab = x, ylab = y, ...){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  requireNamespace("ggpubr", quietly=TRUE) || stop("need ggpubr package")
  df = as.data.frame(data, stringsAsFactors = FALSE)

  ## Add label column in the data frame.
  if(label==0) data$Label = rownames(data)
  else data$Label = as.character(data[, label])

  ## Compute the cutoff used for each dimension.
  model = tolower(model)
  if(length(model)>0){
    if(model == "ninesquare"){
      if(length(x_cut)==0 & auto_cut_x)
        x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
      if(length(y_cut)==0 & auto_cut_y)
        y_cut = c(-CutoffCalling(data[,y], 1.5), CutoffCalling(data[,y], 1.5))
      if(length(intercept)==0 & auto_cut_diag)
        intercept = c(-CutoffCalling(data[,y]-data[,x], 1.5), CutoffCalling(data[,y]-data[,x], 1.5))
    }
    if(model == "volcano"){
      if(length(x_cut)==0 & auto_cut_x)
        x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
      if(length(y_cut)==0 & auto_cut_y) y_cut = -log10(0.05)
    }
    if(model == "rank"){
      if(length(x_cut)==0 & auto_cut_x)
        x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
    }
  }

  ## Decide the colored groups
  if(length(groups)==0){
    if(length(model)>0){
      if(model == "ninesquare") groups = c("midleft", "topcenter", "midright", "bottomcenter")
      if(model == "volcano") groups = c("topleft", "topright")
      if(model == "rank") groups = c("topleft", "topright")
    }
  }
  avail_groups = c("topleft", "topcenter", "topright", "midleft", "midcenter",
                   "midright", "bottomleft", "bottomcenter", "bottomright")
  groups = intersect(groups, avail_groups)

  ## Annotate the groups in the data frame
  idx1 = idx3 = idx5 = FALSE
  idx2 = idx4 = idx6
  if(length(x_cut)>0){
    idx1 = data[,x] < min(x_cut)
    idx2 = data[,x] > max(x_cut)
  }
  if(length(y_cut)>0){
    idx3 = data[,y] < min(y_cut)
    idx4 = data[,y] > max(y_cut)
  }
  if(length(intercept)>0){
    idx5 = data[,y]<slope*data[,x]+min(intercept)
    idx6 = data[,y]>slope*data[,x]+max(intercept)
  }

  data$group="none"
  data$group[idx1 & idx4&idx6] = "topleft"
  data$group[(!idx1)&(!idx2) & idx4&idx6] = "topcenter"
  data$group[idx2 & idx4&idx6] = "topright"
  data$group[idx1&idx6 & (!idx3)&(!idx4)] = "midleft"
  data$group[(!idx1)&(!idx2) & (!idx3)&(!idx4) & (!idx5)&(!idx6)] = "midcenter"
  data$group[idx2 & (!idx3)&(!idx4) & idx5] = "midright"
  data$group[idx1 & idx3 & idx5] = "bottomleft"
  data$group[(!idx1)&(!idx2) & idx3 & idx5] = "bottomcenter"
  data$group[idx2 & idx3 & idx5] = "bottomright"
  data$group[!data$group%in%groups] = "none"
  ## Select the colors
  mycolour=c("#377eb8", "#ff7f00", "#a65628", "#4daf4a", "#005CB7",
             "#e41a1c", "#984ea3", "#f781bf", "#BABABA")
  names(mycolour) = c("topleft", "topcenter", "topright", "midleft", "midright",
                      "bottomleft", "bottomcenter", "bottomright", "none")
  ## Group names
  if(length(groupnames)!=length(groups)) groupnames = groups
  if(length(groups)>0) names(groupnames) = groups
  if(length(group_col)==length(groups)) mycolour[groups] = group_col
  ## Label top gene names ##
  data$rank = top + 1
  for(g in groups){
    idx1 = data$group==g
    data$rank[idx1] = rank(-(abs(data[,x])+abs(data[,y]))[idx1])
  }
  data$Label[data$rank>top & !(data$Label %in% toplabels)] = ""
  data$group=factor(data$group, levels = c(groups, "none"))

  gg = data
  ## Plot the scatter figure ##
  p = ggplot(gg, aes_string(x, y, label="Label", color = "group", fill = "group"))
  p = p + geom_point(shape = 21, alpha=0.8)
  p = p + scale_color_manual(values = mycolour, labels = groupnames)
  p = p + scale_fill_manual(values = mycolour, labels = groupnames)
  if(display_cut){
    if(length(x_cut)>0)
      p = p + geom_vline(xintercept = x_cut,linetype = "dotted")
    if(length(y_cut)>0)
      p = p + geom_hline(yintercept = y_cut,linetype = "dotted")
    if(length(intercept)>0)
      p = p + geom_abline(slope=slope, intercept=intercept, linetype = "dotted")
  }
  p = p + labs(x=xlab, y = ylab, title = main)
  if(label.top)
    p = p + ggrepel::geom_text_repel(...)
  p = p + ggpubr::theme_pubr()
  if(!legend) p = p + theme(legend.position = "none")
  return(p)
}


#' Quantile of normal distribution.
#'
#' Compute cutoff from a normal-distributed vector.
#'
#' @docType methods
#' @name CutoffCalling
#' @rdname CutoffCalling
#'
#' @param d A numeric vector.
#' @param scale Boolean or numeric, specifying how many standard deviation will be used as cutoff.
#'
#' @return A numeric value.
#' @examples
#' CutoffCalling(rnorm(10000))

CutoffCalling=function(d, scale=1){
  param=1
  if(is.logical(scale) & scale){
    param = round(length(d) / 20000, digits = 1)
  }else if(is.numeric(scale)){param = scale}

  Control_mean=0
  sorted_beta=sort(abs(d))
  temp=quantile(sorted_beta,0.68)
  temp_2=qnorm(0.84)
  cutoff=round(temp/temp_2,digits = 3)
  names(cutoff)=NULL
  cutoff=cutoff*param
  return(cutoff)
}
