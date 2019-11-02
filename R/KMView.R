#' Survival analysis and plot the kaplan-meier curves
#' @param survdata A data frame including the clinical information.
#' @param bio A character or an integer, specifying which variable/column
#' in the survdata should be used to predict the risk of patients.
#' @param os A character or an integer, specifying the variable/column of overall survival.
#' @param event A character or an integer, specifying the variable/column of outcome.
#' @param event A character or an integer, specifying the variable/column of outcome.
#' @param cut A numeric vector for separating the patients into different groups,
#' only available when `bio` is a continuous variable.
#' @param optimalCut Boolean, indicating whether determin the ootimal cutpoint.
#' @param pval.method Character, indicating the method for the pvalue calculation,
#' only logrank and cox are available now.
#' @param pval.pos Character, indicating the method for the pvalue calculation.
#' @param labels A character vector, specifying the labels for each patient group.
#' @param ... Other available parameters in ggsurvplot.
#'
#' @return A ggplot instance that allows further custimization.
#'
#' @author Wubing Zhang
#'
#' @import survival survminer
#' @export
KMView <- function(survdata, bio, os = "os", event = "event",
                     cut = c(0, 0.5, 1), optimalCut = TRUE,
                     pval.method = c("cox", "logrank"),
                     pval.pos = c(0.9, 0.95),
                     labels = 1:(length(cut)-1), ...){
  requireNamespace("survival")
  requireNamespace("survminer")
  labels = factor(labels, levels = labels)
  tmpDat <- survdata[!is.na(survdata[, bio]), c(bio, os, event)]
  colnames(tmpDat) = c("bio", "os", "event")
  if(optimalCut){
    tmp <- surv_cutpoint(tmpDat, time = "os", variables = "bio")
    tmpDat<- surv_categorize(tmp, labels = labels[1:2])
    fittedSurv <- survfit(Surv(time=os, event=event) ~ bio, data=tmpDat)
    if(pval.method=="cox"){
      cox <- coxph(Surv(time=os, event=event) ~ bio, data=tmpDat)
      pval <- signif(summary(cox)$coefficients[1,5], digits=3)
    }else{
      diffSurv <- survdiff(Surv(time=os, event=event) ~ bio, data=tmpDat)
      pval <- signif(1 - pchisq(diffSurv$chisq, length(diffSurv$n) - 1), digits=3)
    }
  }else{
    if(is.numeric(tmpDat[1, "bio"]))
      tmpDat$group <- cut(tmpDat[, "bio"],
                          quantile(tmpDat[, "bio"], cut),
                          labels = labels,
                          include.lowest = TRUE)
    tmpDat$group <- factor(tmpDat$group)
    fittedSurv <- survfit(Surv(time=os, event=event) ~ group, data=tmpDat)
    diffSurv <- survdiff(Surv(time=os, event=event) ~ group, data=tmpDat)
    pval <- signif(1 - pchisq(diffSurv$chisq, length(diffSurv$n) - 1), digits=3)
  }
  p = ggsurvplot(fittedSurv, data = tmpDat, surv.median.line = "hv",
                 legend.labs = labels, ...)
  p = p$plot + annotate("text", x = quantile(tmpDat$os, pval.pos[1], na.rm = TRUE),
                        y = pval.pos[2], label = paste("pval:", pval), hjust=0, vjust=1)
  p = p + labs(color = NULL)
  p
}
