#' Extract reads from fastq file
#'
#' @docType methods
#' @name extractReads
#' @rdname extractReads
#'
#' @param fastq_file Path to fastq file.
#'
#' @return A character vector, in which each item is a read sequence.
#' @author Wubing Zhang
#' @export
#'
extractReads <- function(fastq_file){
  gzcon <- gzfile(fastq_file, "rt")
  fq <- readLines(gzcon)
  close(gzcon)
  fa <- fq[seq(2, length(fq), 4)]
  return(fa)
}
