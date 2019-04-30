#' @export
Figure1H <- function(m2, SampleAnn, sig = "sgTRAF3",
                     irf = "Best Confirmed Overall Response",
                     ml = "FMOne mutation burden per MB"){
  require(IMvigor210CoreBiologies)
  data(color_palettes)

  plotDots <- FALSE
  binary <- FALSE
  # font sizes
  labCex <- 0.9
  namesCex <- 0.9
  legendCex <- 0.9
  titleCex <- 1
  axisCex <- 0.9
  titleF <- 1

  orgMar <- par()$mar


  tmpDat <- SampleAnn
  if(! sig %in% colnames(tmpDat))
    tmpDat[, sig] <- scale(m2[sig, ], center=TRUE, scale=TRUE)

  tmpDat <- tmpDat[tmpDat[, irf] != "NE",]
  tmpDat[, irf] <- droplevels(tmpDat[, irf])

  par(mar=c(5.5, 4.1, 2, 0.5))

  pval <- signif(wilcox.test(tmpDat[, sig] ~ tmpDat[, "binaryResponse"])$p.value, 2)
  print(paste("Wilcoxon test P", sig, "score by binary Response:", pval))
  pval <- signif(t.test(tmpDat[, sig] ~ tmpDat[, "binaryResponse"])$p.value, 2)
  print(paste("T test P", sig, "score by binary Response:", pval))

  yLab <- paste(sig, "expression")

  a <- boxplot(tmpDat[, sig] ~ tmpDat[, irf],
               ylab=yLab,
               col=color_palettes[["irf_palette"]][levels(tmpDat[, irf])],
               cex.axis=axisCex,
               cex.names=namesCex,
               cex.lab=labCex,
               ylim=c(min(tmpDat[, sig]), max(tmpDat[, sig])),
               whisklty = 1)
  if (plotDots) {
    points(y=tmpDat[, sig],
           x=jitter(as.numeric(tmpDat[, irf]), factor=0.7),
           col=alpha("darkgrey", 0.7), pch=16)
  }
  mtext(a$n, side=3, at=1:nlevels(tmpDat[, irf]), line=0, cex=0.75)
  mtext(paste0(sig, ", response"), side=3, at=2.5, line=1, font=titleF, cex=titleCex)

  yrange <- par("usr")[4] - par("usr")[3]
  yunit <- yrange/60
  segments(x0=1, x1=2,
           y0=par("usr")[4] - yunit * 4,
           y1=par("usr")[4] - yunit * 4,
           col="black", xpd=TRUE, lwd=1)
  segments(x0=3, x1=4,
           y0=par("usr")[4] - yunit * 4,
           y1=par("usr")[4] - yunit * 4,
           col="black", xpd=TRUE, lwd=1)
  segments(x0=1.5, x1=1.5,
           y0=par("usr")[4] - yunit * 4,
           y1=par("usr")[4] - yunit * 2,
           col="black", xpd=TRUE, lwd=1)
  segments(x0=3.5, x1=3.5,
           y0=par("usr")[4] - yunit * 4,
           y1=par("usr")[4] - yunit * 2,
           col="black", xpd=TRUE, lwd=1)
  segments(x0=1.5, x1=3.5,
           y0=par("usr")[4] - yunit * 2,
           y1=par("usr")[4] - yunit * 2,
           col="black", xpd=TRUE, lwd=1)
  text(y=par("usr")[4] - yunit ,
       x=2.5, labels=pLevel(pval),
       cex=ifelse(pLevel(pval) == "n.s.", axisCex - 0.1, axisCex))
}
