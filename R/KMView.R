#' Survival analysis and plot the kaplan-meier curves
#' @param survdata A data frame including the clinical information.
#' @param bio A character or an integer, specifying which variable/column
#' in the survdata should be used to predict the risk of patients.
#' @param os A character or an integer, specifying the variable/column of overall survival.
#' @param event A character or an integer, specifying the variable/column of outcome.
#' @param adjust A vector of factors to include in the cox model, default c("age", "gender", "race").
#' @param cut A numeric vector for separating the patients into different groups,
#' only available when `bio` is a continuous variable.
#' @param labels A character vector, specifying the labels for each patient group.
#' @param pval.method  logrank (default) or cox, indicating the method for the pvalue calculation.
#' @param pval.pos Character, indicating the method for the pvalue calculation.
#' @param ... Other available parameters in ggsurvplot or ggadjustedcurves.
#'
#' @return A ggplot instance that allows further custimization.
#'
#' @author Wubing Zhang
#'
#' @import survival survminer
#' @export
KMView <- function(survdata, bio, os = "os", event = "event",
                   adjust = c("age", "gender", "race"),
                   cut = c(0, 0.5, 1), labels = NULL,
                   pval.method = c("cox", "logrank"),
                   pval.pos = c(0, 0.2),
                   ...){
  requireNamespace("survival")
  requireNamespace("survminer")
  # Remove variables with too many NAs
  survdata = survdata[, colSums(is.na(survdata))<nrow(survdata)/2]
  # Remove samples with NAs
  idx = !(is.na(survdata[,os])|is.na(survdata[,event])|is.na(survdata[,bio]))
  survdata = survdata[idx, ]
  # Remove variables with unique value
  idx = apply(survdata, 2, function(x) sum(!is.na(unique(x))))
  survdata = survdata[, idx>1]
  tmpDat = survdata[, colnames(survdata)%in%c(os, event, bio, adjust)]
  tmpDat = tmpDat[, c(os, event, bio, setdiff(colnames(tmpDat), c(os, event, bio)))]
  colnames(tmpDat)[1:2] = c("os", "event")

  if(is.numeric(tmpDat[1, bio])){
    interestTerm = bio
    if(length(unique(tmpDat[, bio]))>3) tmpDat[, bio] = scale(tmpDat[, bio])
    cox <- coxph(Surv(time=os, event=event) ~ ., data=tmpDat)
    breaks = quantile(tmpDat[,interestTerm], cut)
    ## Generate labels automatically
    if(is.null(labels)) labels = paste0("level", 1:(length(breaks)-1))
    labels = factor(labels, levels = labels)

    ## Group all samples
    tmpDat[,interestTerm] = cut(tmpDat[,interestTerm], breaks,
                                labels = labels, include.lowest = TRUE)

    ## Calculate the p-value for plotting
    tmpres <- summary(cox)
    pval <- signif(tmpres$coefficients[interestTerm,5], 3)
    hr <- round(tmpres$conf.int[interestTerm,1], 3)
    hr.low <- round(tmpres$conf.int[interestTerm,3], 3)
    hr.high <- round(tmpres$conf.int[interestTerm,4], 3)
  }else{
    if(is.null(labels)) labels = unique(tmpDat[, bio])
    names(labels) = unique(tmpDat[, bio])
    tmpDat[, bio] = as.factor(tmpDat[, bio])
    interestTerm = bio

    ## Calculate the p-value for plotting
    cox <- coxph(Surv(time=os, event=event) ~ ., data=tmpDat)
    idx = grepl(interestTerm, rownames(summary(cox)$coefficients))
    tmpres <- summary(cox)
    pval <- signif(tmpres$coefficients[idx,5], 3)
    hr <- round(tmpres$conf.int[idx,1], 3)
    hr.low <- round(tmpres$conf.int[idx,3], 3)
    hr.high <- round(tmpres$conf.int[idx,4], 3)
  }
              
  ## Calculate the p-value for plotting
  if(pval.method[1]=="logrank"){
    if(any(adjust%in%colnames(tmpDat))){
      warning("Coxph p-value is used because of the adjustment factors in the analysis!!!")
    }else{
      errflag = FALSE
      diffSurv <- tryCatch(survdiff(Surv(time=os, event=event) ~ ., data=tmpDat),
                           error = function(e) errflag <<- TRUE)
      if(errflag) return(1)
      pval <- signif(1 - pchisq(diffSurv$chisq, length(diffSurv$n) - 1), digits=3)
      cox <- coxph(Surv(time=os, event=event) ~ ., data=tmpDat)
      tmpres <- summary(cox)
      hr <- round(tmpres$conf.int[1,1], 3)
      hr.low <- round(tmpres$conf.int[1,3], 3)
      hr.high <- round(tmpres$conf.int[1,4], 3)
    }
  }

  ## Plotting the KM-plot
  if(any(adjust%in%colnames(tmpDat))){
    newcox = coxph(Surv(time=os, event=event) ~ ., data=tmpDat)
    p = ggadjustedcurves(newcox, variable = interestTerm, data = tmpDat,
                         legend.labs = labels, ...)
    p = p + annotate("text", x = quantile(tmpDat$os, pval.pos[1], na.rm = TRUE), y = pval.pos[2], 
                     label = paste0("P = ", pval, "\nHR = ", 
                                   hr, " [", hr.low, ", ", hr.high, "]"), hjust=0, vjust=1)
    p = p + labs(color = NULL)
    p
  }else{
    colnames(tmpDat)[colnames(tmpDat)==interestTerm] = "interestTerm"
    fittedSurv <- survfit(Surv(time=os, event=event) ~ interestTerm,
                          data=tmpDat[, c("os", "event", "interestTerm")])
    p = ggsurvplot(fittedSurv, data = tmpDat, surv.median.line = "hv",
                   legend.labs = labels, ...)
    p = p$plot + annotate("text", x = quantile(tmpDat$os, pval.pos[1], na.rm = TRUE), y = pval.pos[2], 
                          label = paste0("P = ", pval, "\nHR = ", 
                                        hr, " [", hr.low, ", ", hr.high, "]"), hjust=0, vjust=1)
    p = p + labs(color = NULL)
    p
  }
  return(p)
}
