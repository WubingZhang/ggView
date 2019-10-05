#' @export
SurvView <- function(SampleAnn, bio = "sgTRAF3",
                     os = "os", event = "event", group = 4,
                     glabels = c("1st Quart", "2nd Quart", "3rd Quart", "4th Quart"),
                     main="CD8 T-effector, survival",
                     plotMedian=FALSE,
                     plot.nrisk=FALSE){
  require(survival)
  library(IMvigor210CoreBiologies)
  require(ggView)
  data(color_palettes)

  labCex <- 0.9
  namesCex <- 0.9
  legendCex <- 0.9
  titleCex <- 1
  axisCex <- 0.9
  titleF <- 1
  tmpDat <- SampleAnn[!is.na(SampleAnn[, bio]), c(bio, os, event)]
  colnames(tmpDat) = c(bio, "os", "event")
  tmpDat$group <- cut(tmpDat[, bio],
                      quantile(tmpDat[, bio], seq(0, 1, by = 1/group)),
                      include.lowest = TRUE)
  tmpDat$group <- factor(tmpDat$group, labels=glabels)
  fittedSurv <- survfit(Surv(time=os, event=event) ~ group, data=tmpDat)
  diffSurv <- survdiff(Surv(time=os, event=event) ~ group, data=tmpDat)
  ggView::plotSurvival2(survFit=fittedSurv,
                survDiff=diffSurv,
                diff.factor=as.factor(tmpDat$group),
                main=main,
                cols=rev(color_palettes[["irf_palette"]][round(seq(1,4,length.out = group))]),
                ltypes=rep(1, nlevels(tmpDat$group)),
                pval="logrank",
                #coxphData=coxphDat,
                plotMedian=plotMedian,
                plot.nrisk=plot.nrisk,
                xLab="OS (months)",
                lwd=3,
                cexMedian=axisCex - 0.1,
                cexLegend=legendCex,
                cex.lab=labCex)

  print(coxph(Surv(os, event) ~ as.integer(tmpDat$group), data=tmpDat))
}
