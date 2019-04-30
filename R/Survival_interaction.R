#' Survival interaction in TIDE
#' @import survival
#' @export
#'
Survival_interaction <- function(mat, # expression, mutation, or RPPA matrix
                                 pivot,  # CD8A, CD8B, GZMA, GZMB, PRF1 expression
                                 survival, # clinical information
                                 CTL_genes = c('CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1')){
  suppressPackageStartupMessages(library(survival))

  # remove negative survival values
  survival = survival[survival[,1] > 0,]

  common = Reduce(intersect, list(rownames(mat),rownames(survival),rownames(pivot)))
  print(paste(length(common), "samples"))

  # stop at low sample numbers
  if(length(common) < 20) stop("two few samples")

  not_included = setdiff(CTL_genes, colnames(pivot))

  if(length(not_included) > 0) stop(paste(c("pivot genes", not_included, "are not included"), ' '))

  pivot = rowMeans(pivot[common,CTL_genes])
  mat = mat[common,,drop=F]
  survival = survival[common,,drop=F]

  # stop at low death rate
  death_rate = sum(survival[,2])/dim(survival)[1]
  if(death_rate < 0.1) stop(paste("death rate", death_rate, "is too low"))

  # split up survival and background
  surv = Surv(survival[,1], survival[,2])

  if(dim(survival)[2] > 2){
    B = survival[,3:dim(survival)[2], drop=F]
  }else{
    B = survival[,c(), drop=F]
  }

  # build up regression data space
  B = cbind(B, pivot, pivot, pivot)
  B = as.data.frame(B)
  N_B = dim(B)[2]

  # test pivot effect first
  coxph.pivot = summary(coxph(surv~., data=B[,c(-N_B, -(N_B-1)), drop=F]))$coef

  colnames(B)[N_B-1] = "partner"
  colnames(B)[N_B] = "Interaction"

  # iterate over features
  features = colnames(mat)
  N = length(features)

  result_interaction = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
  result_main = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
  result_partial = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
  result_base = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))

  for (i in 1:N)
  {
    title = features[i]
    partner = mat[,i]

    B[,N_B-1] = partner
    B[,N_B] = partner * pivot

    # region 1: model with interaction
    errflag = F
    coxph.fit = tryCatch(coxph(surv~., data=B),
                         error = function(e) errflag <<- T,
                         warning = function(w) errflag <<- T)

    if(!errflag)
    {
      reg.summary = summary(coxph.fit)$coef
      result_interaction[i,] = reg.summary["Interaction", c("z", "Pr(>|z|)")]
      result_main[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
    }

    # region 2: model without interaction
    errflag = F
    coxph.fit = tryCatch(coxph(surv~., data=B[,-N_B]),
                         error = function(e) errflag <<- T,
                         warning = function(w) errflag <<- T)

    if(!errflag)
    {
      reg.summary = summary(coxph.fit)$coef
      result_partial[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
    }

    # region 3: model with only base effect
    errflag = F
    coxph.fit = tryCatch(coxph(surv~., data=B[,c(-N_B, -(N_B-2)), drop=F]),
                         error = function(e) errflag <<- T,
                         warning = function(w) errflag <<- T)

    if(!errflag)
    {
      reg.summary = summary(coxph.fit)$coef
      result_base[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
    }
  }

  # Adjust pvalue.
  result_interaction = result_interaction[!is.na(result_interaction[,"p"]),,drop=F]
  result_interaction = cbind(result_interaction, p.adjust(result_interaction[,"p"], method="fdr"))
  result_main = result_main[!is.na(result_main[,"p"]),,drop=F]
  result_main = cbind(result_main, p.adjust(result_main[,"p"], method="fdr"))
  result_partial = result_partial[!is.na(result_partial[,"p"]),,drop=F]
  result_partial = cbind(result_partial, p.adjust(result_partial[,"p"], method="fdr"))
  result_base = result_base[!is.na(result_base[,"p"]),,drop=F]
  result_base = cbind(result_base, p.adjust(result_base[,"p"], method="fdr"))

  return(list(coxph.pivot = coxph.pivot, result_interaction = result_interaction,
              result_main = result_main, result_partial = result_partial,
              result_base = result_base))
}
