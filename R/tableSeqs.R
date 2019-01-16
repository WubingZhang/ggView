#' Table creation for sequences
#'
#' Build a contingency table of the counts at each position of sequences.
#'
#' @docType methods
#' @name tableSeqs
#' @rdname tableSeqs
#'
#' @param sequences A character vector with the same length of sequences.
#' @param len Integer, specifying the read length.
#' @return A data frame.
#' @author Wubing Zhang
#' @export
#'
tableSeqs <- function(sequences, len = nchar(sequences[1])){
  sequences = sequences[nchar(sequences)==len]
  tmp = as.data.frame(strsplit(sequences, ""))
  tmp = as.matrix(tmp)
  gg = c()
  for(item in unique(c(tmp))){
    gg = cbind(gg, rowSums(tmp==item))
    gg = as.data.frame(gg)
  }
  colnames(gg) = unique(c(tmp))
  return(gg)
}
