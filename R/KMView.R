#' Survival analysis and plot the kaplan-meier curves
#' @param survdata A data frame including the clinical information.
#' @param bio A character or an integer, specifying which variable/column
#' in the survdata should be used to predict the risk of patients.
#' @param os A character or an integer, specifying the variable/column of overall survival.
#' @param event A character or an integer, specifying the variable/column of outcome.
#' @param adjust A vector of factors to include in the cox model, default c("age", "gender", "race").
#' @param interact A factor that interact with `bio`.
#' @param cut A numeric vector for separating the patients into different groups,
#' only available when `bio` is a continuous variable.
#' @param labels A character vector, specifying the labels for each patient group.
#' @param optimalCut Boolean, indicating whether determin the ootimal cutpoint.
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
                   adjust = c("age", "gender", "race"), interact = NULL,
                   cut = c(0, 0.5, 1), labels = NULL,
                   optimalCut = TRUE,
                   pval.method = c("logrank", "cox"),
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
    ## Add interaction term into the tmpDat
    if(length(interact)!=0 && interact[1]%in%colnames(tmpDat)){
      tmpDat$Interaction = tmpDat[,bio] * tmpDat[,interact[1]]
      tmpDat$Partner = tmpDat[,interact[1]]
    }
    interestTerm = bio
    if("Interaction" %in% colnames(tmpDat)) interestTerm = "Interaction"
    cox <- coxph(Surv(time=os, event=event) ~ ., data=tmpDat)
    ## Select the cut points
    if(optimalCut | is.null(cut)){
      if(any(adjust%in%colnames(tmpDat))|"Interaction"%in%colnames(tmpDat)){
        ## Select the best cut point based on the log-rank test.
        perm_pvals = lapply(seq(0.1, 0.9, length.out = 100), function(x){
          tmpCut = quantile(tmpDat[,interestTerm], x)
          tmpDat[, interestTerm][tmpDat[, interestTerm]<tmpCut] = "low"
          tmpDat[, interestTerm][tmpDat[, interestTerm]!="low"] = "high"
          tmpDat[, interestTerm] = factor(tmpDat[, interestTerm],
                                          levels = c("low", "high"))
          errflag = FALSE
          diffSurv <- tryCatch(survdiff(Surv(time=os, event=event) ~ ., data=tmpDat),
                                 error = function(e) errflag <<- TRUE)
          if(errflag) return(1)
          pval <- signif(1 - pchisq(diffSurv$chisq, length(diffSurv$n) - 1), digits=3)
        })
        cut = c(0, seq(0.1, 0.9, length.out = 100)[which.min(unlist(perm_pvals))], 1)
        breaks = quantile(tmpDat[,interestTerm], cut)
      }else{
        tmp <- surv_cutpoint(tmpDat, time = "os", event="event", variables = interestTerm)
        breaks = c(min(tmpDat[,interestTerm]), tmp$cutpoint$cutpoint, max(tmpDat[,interestTerm]))
      }
    }
    ## Generate labels automatically
    if(is.null(labels)) labels = paste0("level", 1:(length(breaks)-1))
    labels = factor(labels, levels = labels)

    ## Group all samples
    tmpDat[,interestTerm] = cut(tmpDat[,interestTerm], breaks,
                                labels = labels, include.lowest = TRUE)

    ## Calculate the p-value for plotting
    pval <- signif(summary(cox)$coefficients[interestTerm,5], 3)
  }else{
    if(is.null(labels)) labels = unique(tmpDat[, bio])
    names(labels) = unique(tmpDat[, bio])
    tmpDat[, bio] = as.factor(tmpDat[, bio])
    ## Add interaction term into the tmpDat
    if(length(interact)!=0 && interact[1]%in%colnames(tmpDat)){
      tmpDat$Partner = as.factor(tmpDat[,interact[1]])
      tmpDat$Interaction = as.numeric(tmpDat[,bio]) * as.numeric(tmpDat$Partner)
    }
    interestTerm = bio
    if("Interaction" %in% colnames(tmpDat)) interestTerm = "Interaction"

    ## Calculate the p-value for plotting
    cox <- coxph(Surv(time=os, event=event) ~ ., data=tmpDat)
    idx = grepl(interestTerm, rownames(summary(cox)$coefficients))
    pval <- signif(min(summary(cox)$coefficients[idx,5]), digits=3)
  }
  ## Calculate the p-value for plotting
  if(pval.method[1]=="logrank"){
    if(any(adjust%in%colnames(tmpDat))){
      warning("Coxph p-value is used because of the adjustment factors in the analysis!!!")
      break
    }
    if("Interaction"%in%colnames(tmpDat)){
      warning("Coxph p-value is used because of the interaction term in the analysis!!!")
      break
    }
    errflag = FALSE
    diffSurv <- tryCatch(survdiff(Surv(time=os, event=event) ~ ., data=tmpDat),
                         error = function(e) errflag <<- TRUE)
    if(errflag) return(1)
    pval <- signif(1 - pchisq(diffSurv$chisq, length(diffSurv$n) - 1), digits=3)
  }

  ## Plotting the KM-plot
  if(any(adjust%in%colnames(tmpDat)) | "Interaction"%in%colnames(tmpDat)){
    newcox = coxph(Surv(time=os, event=event) ~ ., data=tmpDat)
    p = ggadjustedcurves(newcox, variable = interestTerm, data = tmpDat,
                         legend.labs = labels, ...)
    p = p + annotate("text", x = quantile(tmpDat$os, pval.pos[1], na.rm = TRUE),
                     y = pval.pos[2], label = paste("pval:", pval), hjust=0, vjust=1)
    p = p + labs(color = NULL)
    p
  }else{
    fittedSurv <- survfit(Surv(time=os, event=event) ~ tmpDat[,interestTerm], data=tmpDat)
    p = ggsurvplot(fittedSurv, data = tmpDat, surv.median.line = "hv",
                   legend.labs = labels, ...)
    p = p$plot + annotate("text", x = quantile(tmpDat$os, pval.pos[1], na.rm = TRUE),
                          y = pval.pos[2], label = paste("pval:", pval), hjust=0, vjust=1)
    p = p + labs(color = NULL)
    p
  }
  return(p)
}
