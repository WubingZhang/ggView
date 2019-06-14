#' Basic function to convert mouse to human gene names
#' @param x A vector of genes.
#' @param idType Entrez Or Symbol
#' @param ref Reference table
#' @export
convertMouseGeneList <- function(x, idType = "Symbol",
                                 ref = "~/Jobs/Archive/GeneLists/GeneAnnotation/HOM_MouseHumanSequence.txt"){
  if(!file.exists(ref))
    data.table::fread("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
                      header = TRUE, stringsAsFactors = FALSE) -> reftable
  else
    read.table(ref, sep = "\t", header = TRUE, stringsAsFactors = FALSE) -> reftable
  colnames(reftable)[5] = "Entrez"
  idx1 = reftable[, idType] %in% x
  idx2 = grepl("human", reftable$Common.Organism.Name)
  idx3 = reftable$HomoloGene.ID %in% reftable$HomoloGene.ID[idx1]
  tmp1 = reftable[idx1, c("HomoloGene.ID", idType)];
  tmp = reftable[(idx2&idx3), c("HomoloGene.ID", idType)]
  colnames(tmp1)[2] = "Mouse"; colnames(tmp)[2] = "Human"
  res = merge(tmp1, tmp, by = "HomoloGene.ID")[,-1]

  # t(sapply(unique(reftable$HomoloGene.ID), function(i){
  #   idx1 = reftable$HomoloGene.ID==i
  #   idx2 = grepl("mouse", reftable$Common.Organism.Name)
  #   idx3 = grepl("human", reftable$Common.Organism.Name)
  #   c(ID=i, mSymbol=paste0(reftable$Symbol[idx1&idx2], collapse = ","),
  #     hSymbol = paste0(reftable$Symbol[idx1&idx3], collapse = ","),
  #     mEntrez=paste0(reftable$Entrez[idx1&idx2], collapse = ","),
  #     hEntrez = paste0(reftable$Entrez[idx1&idx3], collapse = ","))
  # })) -> newReftable
  # newReftable = as.data.frame(newReftable, stringsAsFactors = FALSE)
  # saveRDS(newReftable, "~/Jobs/Archive/GeneLists/GeneAnnotation/HOM_MouseHuman.rds")
  return(res)
}

#' Basic function to convert human to mouse gene names
#' @param x A vector of genes.
#' @param idType Entrez Or Symbol
#' @param ref Reference table
#' @export
convertHumanGeneList <- function(x, idType = "Symbol",
                                 ref = "~/Jobs/Archive/GeneLists/GeneAnnotation/HOM_MouseHumanSequence.txt"){
  if(!file.exists(ref))
    data.table::fread("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
                      header = TRUE, stringsAsFactors = FALSE) -> reftable
  else
    read.table(ref, sep = "\t", header = TRUE, stringsAsFactors = FALSE) -> reftable
  colnames(reftable)[5] = "Entrez"

  idx1 = reftable[, idType] %in% x
  idx2 = grepl("mouse", reftable$Common.Organism.Name)
  idx3 = reftable$HomoloGene.ID %in% reftable$HomoloGene.ID[idx1]
  tmp1 = reftable[idx1, c("HomoloGene.ID", idType)]
  tmp = reftable[(idx2&idx3), c("HomoloGene.ID", idType)]
  colnames(tmp1)[2] = "Human"; colnames(tmp)[2] = "Mouse"
  res = merge(tmp1, tmp, by = "HomoloGene.ID")[,-1]
  return(res)
}
