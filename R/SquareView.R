#' Scatter View of MHC-I and PD-L1 level, and separate all samples by 9 square.
#'
#' @docType methods
#' @name SquareView
#' @rdname SquareView
#'
#' @param mat Data frame, each row is gene and each column is sample.
#' @param x Specifying x-axis which should be valid col names.
#' @param y Specifying y-axis which should be valid col names.
#' @param option 1 or 2.
#' @param label Specifying label text which should be valid col name.
#' @param top The number of top terms to be labeled.
#' @param genelist Interested term sets to be labeled.
#' @param x_cutoff Specifying cutoff to classify x into high and low groups.
#' @param y_cutoff Specifying cutoff to classify x into high and low groups.
#' @param groups Specifying name of each group.
#'
#' @return ggplot object.
#'
#' @author Wubing Zhang
#' @import ggplot2
#' @export

SquareView<-function(mat, x="Ctrl.Mean", y="IFNr.Mean", option = 1, label = 0,
                     top=10, genelist=c(), x_cutoff=c(-1,1), y_cutoff=c(-1,1),
                     groups=c("DN", "MP_PN", "DP", "MN_PP")){
  requireNamespace("ggExtra", quietly=TRUE) || stop("need ggExtra package")
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  mat = as.data.frame(mat, stringsAsFactors = FALSE)
  if(label==0)
    mat$Gene = rownames(mat)
  else
    mat$Gene = mat[, label]
  mat$x = rowMeans(mat[, x, drop = FALSE], na.rm = TRUE)
  mat$y = rowMeans(mat[, y, drop = FALSE], na.rm = TRUE)
  mat$group="Others"
  mat$text = mat$Gene

  idx1 = mat$x < x_cutoff[1]
  idx2 = mat$x > x_cutoff[2]
  idx3 = mat$y < y_cutoff[1]
  idx4 = mat$y > y_cutoff[2]
  if(option == 1){
    # PN
    mat$rank = top + 1
    idx = (!idx1)&(!idx2)&idx4
    mat$group[idx] = groups[1]
    mat$rank[idx] = rank((-mat$y)[idx])
    # MP
    idx = (!idx3)&(!idx4)&idx1
    mat$group[idx] = groups[2]
    mat$rank[idx] = rank(mat$x[idx])
    # PP
    idx = (!idx1)&(!idx2)&idx3
    mat$group[idx] = groups[3]
    mat$rank[idx] = rank(mat$y[idx])
    # MN
    idx = (!idx3)&(!idx4)&idx2
    mat$group[idx] = groups[4]
    mat$rank[idx] = rank(-mat$x[idx])
    mat$text[mat$rank>top] = NA
  }
  if(option == 2){
    # Double-positive regulators
    mat$rank = top + 1
    idx = idx1&idx3
    mat$group[idx] = groups[3]
    mat$rank[idx] = rank((mat$x+mat$y)[idx])
    # Double-negative regulators
    idx = idx2&idx4
    mat$group[idx] = groups[1]
    mat$rank[idx] = rank(-1*(mat$x+mat$y)[idx])
    # MHCI-positive PD-L1 negative regulators
    idx = idx1&idx4
    mat$group[idx] = groups[2]
    mat$rank[idx] = rank((mat$x-mat$y)[idx])
    # MHCI-negative PD-L1 positive regulators
    idx = idx2&idx3
    mat$group[idx] = groups[4]
    mat$rank[idx] = rank((mat$y-mat$x)[idx])
    mat$text[mat$rank>top] = NA
  }

  gg = mat
  gg$group=factor(gg$group, levels = c(groups,"Others"))
  #===============
  mycolour=c("#4daf4a", "#984ea3", "#ff7f00", "#005CB7", "aliceblue")
  idx1 = gg$Gene %in% genelist
  idx2 = !(gg$group=="Others" | is.na(gg$text))
  label_gg = gg[idx1|idx2,]
  col_label = rep("#004b84", nrow(label_gg))
  col_label[label_gg$group=="Others"]="gray60"
  p = ggplot(gg,aes(x=x, y=y, colour=group,fill=group))
  p = p + geom_point(shape=".",alpha=1/1,size = 1) + scale_color_manual(values=mycolour)
  p = p + geom_jitter(size = 1)
  p = p + geom_vline(xintercept = x_cutoff,linetype = "dotted")
  p = p + geom_hline(yintercept = y_cutoff,linetype = "dotted")
  p = p + guides(col = guide_legend(ncol = 3, byrow = TRUE))
  p = p + labs(x=x[1], y=y[2])
  if(top>0)
    p = p + ggrepel::geom_text_repel(aes(x=x,y=y,label=Gene), color=col_label, data=label_gg)

  p = p + annotate("text", color="red", x = x_cutoff[1], y=min(gg$y), vjust=-2, hjust = 1,
                   label=paste0(groups[3], ": ", length(which(mat$group==groups[3]))))
  p = p + annotate("text", color="red", x = x_cutoff[1], y=max(gg$y), vjust=2, hjust = 1,
                   label=paste0(groups[2], ": ", length(which(mat$group==groups[2]))))
  p = p + annotate("text",color="red", x = x_cutoff[2], y = min(gg$y), vjust=-2, hjust = -0.5,
                   label=paste0(groups[4], ": ", length(which(mat$group==groups[4]))))
  p = p + annotate("text",color="red", x = x_cutoff[2], y = max(gg$y), vjust=2, hjust = -0.5,
                   label=paste0(groups[1], ": ", length(which(mat$group==groups[1]))))
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.position="none")+theme(legend.title=element_blank())
  return(p)
}
