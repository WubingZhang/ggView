##' Custom version of Kaplan Meier plot
##'
##' Custom display of survival data, including plotting lines for median survival and adding number at risk
##' @param survFit survfit object; survfit(Surv(TTE, Cens) ~ group, data=df)
##' @param survDiff surfdiff object; survdiff(Surv(TTE, Cens) ~ group, data=df)
##' @param diff.factor factor, the grouping factor used for survfit call
##' @param main character; main title for plot
##' @param cols character vector; colors, default to darkgreen, darkmagenta, cyan4, darkorange, darkred
##' @param xLab character; x-axis label for plot
##' @param yLab character; y-axis label for plot
##' @param ltypes numeric vector; line types (lty) for plot
##' @param plotMedian logical; shall the median survival for each group be plotted
##' @param pval type of p-value computed, either "logrank" or "coxph" (the latter
##' can only be computed on two groups)
##' @param coxphData coxph object; coxphData <- coxph(Surv(TTE, Cens) ~ group, data=df)
##' @param legPos legend position for legend()
##' @param mar plotting margins
##' @param plot.nrisk logical; shall number of samples for each group and time point
##' be plotted underneath graph
##' @param nrisk.interval numeric; spacing of time intervals for which number of samples
##'                is given; defaults to 2
##' @param timemark logical; if set to TRUE (default), censoring marks are added to survival curves
##' @param lwd line width for survival curves
##' @param cexMedian numeric; relative size of font for median survival times
##' @param cexLegend numeric; relative size of legend font
##' @param ... other arguments passed to plot.survfit
##' @author Yuanyuan Xiao, Dorothee Nickles
##' @export
##' @import survival
plotSurvival2 <- function(survFit,
                          survDiff,
                          diff.factor,
                          main,
                          cols,
                          xLab="PFS (months)",
                          yLab="Probability of survival",
                          ltypes=1,
                          plotMedian=FALSE,
                          pval=c("logrank", "coxph", "none"),
                          coxphData,
                          legPos="topright",
                          mar=c(12,9,3,2),
                          plot.nrisk=FALSE,
                          nrisk.interval=2,
                          timemark=TRUE,
                          lwd=2,
                          cexMedian=0.75,
                          cexLegend=0.8,
                          ...) {


  stopifnot(is(survFit, "survfit"))
  stopifnot(is(survDiff, "survdiff"))
  stopifnot(is(diff.factor, "factor"))

  ## simple checks to decrease chance that fitted data used indeed
  ## same data and same calls
  fitCall <- as.character(survFit$call)
  diffCall <- as.character(survDiff$call)
  stopifnot(fitCall[2] == diffCall[2])
  stopifnot(fitCall[3] == diffCall[3])

  group.labels <- levels(diff.factor)

  ## setting colors
  if (missing(cols)) {
    cols <- c("darkgreen", "darkmagenta", "cyan4", "darkorange", "darkred")[1:nlevels(diff.factor)]
  }
  stopifnot(length(cols) == length(group.labels))

  ## computing number of samples at each time point
  if (plot.nrisk) {
    time.pt <- seq(0, max(survFit$time), nrisk.interval)
    ix = 0
    n.risk = c()
    for (kk in 1:(length(survFit$strata))) {
      fit.n.risk = survFit$n.risk[(ix+1) : (ix + survFit$strata[kk])]
      fit.time = survFit$time[(ix+1) : (ix + survFit$strata[kk])]
      tmp = findInterval(time.pt, fit.time)
      n.risk=rbind(n.risk,
                   ifelse(tmp < length(fit.time), fit.n.risk[tmp+1], 0)
      )
      ix = ix + survFit$strata[kk]
    }
    dimnames(n.risk)[[2]] = time.pt

    if (mar[1]<4+length(group.labels)) {
      mar[1] <- 4+length(group.labels)
    }
    org.mar <- par()$mar
    par(mar=mar)
  }

  ## plot curves and axes
  plot(survFit,
       xlab=xLab,
       lwd=lwd,
       lty=ltypes,
       col=cols,
       ylab=yLab,
       main=main,
       axes=FALSE,
       mark.time=timemark,
       ...)
  box()
  axis(1,
       at=seq(0, max(survFit$time), nrisk.interval),
       seq(0, max(survFit$time), nrisk.interval)
  )
  axis(2,
       at=seq(0, max(survFit$time), 0.1),
       seq(0, max(survFit$time), 0.1),
       las=2)
  abline(h=0,
         col="darkgrey")

  ## add the survival medians
  if (plotMedian) {
    median.surv <- summary(survFit)$table[,"median"]
    wMed <- which(!is.na(median.surv))
    median.surv <- median.surv[!is.na(median.surv)]
    median.surv <- round(median.surv, digits=2)
    for (i in 1:length(median.surv)) {
      lines(x=rep(median.surv[i], 6),
            y=seq(0,0.5,0.1),
            lty=2,
            col="darkgrey")
      lines(x=seq(0, median.surv[i], length.out=6),
            y=rep(0.5,6),
            lty=2,
            col="darkgrey")
      text(x=median.surv[i],
           y=-0.02,
           labels=median.surv[i],
           cex=cexMedian,
           col=cols[wMed][i])
    }
  }

  ## add the number of samples underneath plot
  if (plot.nrisk) {
    for (i in 1:length(group.labels)) {
      mtext(side=1,
            at=-1.7,
            line=i+3.5,
            text=group.labels[i],
            col=cols[i],
            adj=1,
            cex=0.6)
      mtext(side=1,
            at=time.pt,
            line=i+3.5,
            text=n.risk[i,],
            col=cols[i],
            cex=0.8)
    }
  }

  ## compute statistics
  if (pval == "logrank") {
    pval <- signif(1 - pchisq(survDiff$chisq, length(survDiff$n) - 1), digits=2)
    group.labels <- c(group.labels, paste("log rank pval:", pval))
  }
  if (pval == "coxph") {
    stopifnot(!missing(coxphData))
    stopifnot(is(coxphData, "coxph"))
    coxCall <- as.vector(coxphData$formula)
    stopifnot(fitCall[2] == coxCall)
    stopifnot(sum(survFit$n) == coxphData$n)
    stopifnot(length(group.labels) == 2)

    tmp <- round(getHRandCIfromCoxph(coxphData),
                 digits=4)
    group.labels <- c(group.labels,
                      paste("pval:", tmp[1,"P"]),
                      paste("HR:",
                            paste0(tmp[1,"HR"], " (", tmp[1,"CIl0.95"], ";",
                                   tmp[1,"CIu0.95"], ")")))
  }

  ## plot legend
  legend(legPos,
         group.labels,
         lwd=2,
         col=c(cols, "white", "white"),
         lty=ltypes,
         legend=group.labels,
         bty="n",
         cex=cexLegend)

  ## reset mar that was changed to allow adding numbers underneath plot
  if (plot.nrisk) {
    par(mar=org.mar)
  }

}
