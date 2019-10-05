#' Synergy analysis
#'
getSynergy <- function(data, base, partner, Time = "Time", Status = "Status"){
  suppressPackageStartupMessages(library(survival))
  data = as.data.frame(data, stringsAsFactors = FALSE)
  if(!all(c(base, partner, Time, Status) %in% colnames(data))){
    stop("Features aren't included in the data:",
         paste(setdiff(c(base, partner, Time, Status), colnames(data))))
  }
  data[, Time] = as.numeric(data[, Time])
  data[, Status] = as.numeric(data[, Status])
  data = data[!(is.na(data[,Time])|is.na(data[,Status])), ]
  surv = Surv(data[,Time], data[,Status])
  # fit_base = summary(coxph(surv~., data=data[, c(base), drop = FALSE]))$coef
  # fit_partner = summary(coxph(surv~., data=data[, c(base, partner)]))$coef
  data[[paste0(base, "_", partner)]] = data[, base] * data[, partner]
  idx = c(base, partner, paste0(base, "_", partner))
  errflag = FALSE
  fit_interaction = tryCatch(summary(coxph(surv~., data=data[, idx]))$coef,
                             error = function(e) errflag <<- TRUE,
                             warning = function(w) errflag <<- TRUE)
  if(!errflag){
    rownames(fit_interaction) = gsub("`", "", rownames(fit_interaction))
    res = as.vector(t(fit_interaction[idx, c("z", "Pr(>|z|)")]))
  }else{
    res = rep(NA, 6)
  }
  names(res) = paste0(rep(c("base", "partner", "interaction"), each=2), rep(c(".z", ".p"), 2))
  return(res)
}
