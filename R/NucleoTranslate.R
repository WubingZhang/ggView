#' Nucleotide Sequence Translation
#' Translate nucleic acid sequence to corresponding peptide sequences.
#'
#' @docType methods
#' @name NucleoTranslate
#' @rdname NucleoTranslate
#'
#' @param sequence A character, specifying the nucleic acid sequence.
#'
#' @return A character of peptide sequence.
#'
#' @examples
#' sequence = "AATGAGGTTGGTGCGGAGTTTCATCCTACTATGTTG"
#' NucleoTranslate(sequence)
#' @author Wubing Zhang
#' @export

NucleoTranslate <- function(sequence){
  # Read DNA codon table
  CodonPath = file.path(system.file("extdata", package = "ggView"),
                        "DNACodonTable.gz")
  CodonFile = gzfile(CodonPath)
  CodonTable = read.table(CodonFile, sep = "\t", header = TRUE,
                          row.names = 1, stringsAsFactors = FALSE)

  # Separate DNA codons from sequence and translate
  tmp = unlist(strsplit(sequence, ""))
  codons <- sapply(seq(3, length(tmp), 3), function(x){
    paste0(tmp[(x-2):x], collapse = "")
  })
  Peptide = CodonTable[codons, "SLC"]
  Peptide = Peptide[!is.na(Peptide)]
  idx = length(Peptide)
  if(idx>0){
    if(sum(Peptide=="Stop")>0) idx = which(Peptide=="Stop")[1]-1
    Peptide = paste0(Peptide[1:idx], collapse = "")
  }
  return(Peptide)
}
