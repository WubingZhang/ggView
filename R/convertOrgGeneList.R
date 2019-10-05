#' Basic function to convert mouse to human gene names
#' @param x A vector of genes.
#' @param idType Entrez Or Symbol
#' @param ref Reference table
#' @import biomaRt
#' @export
convertMouseGeneList <- function(x, idType = "Symbol",
                                 ref = "~/Jobs/Archive/GeneLists/GeneAnnotation/HOM_MouseHumanSequence.txt"){
  if(!file.exists(ref))
    data.table::fread("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
                      header = TRUE, stringsAsFactors = FALSE) -> reftable
  else
    read.table(ref, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) -> reftable
  reftable = as.data.frame(reftable, stringsAsFactors = FALSE)
  colnames(reftable)[5] = "Entrez"
  idx1 = reftable[, idType] %in% x
  idx2 = grepl("human", reftable[,2])
  idx3 = reftable[,1] %in% reftable[idx1,1]
  tmp1 = reftable[idx1, c("HomoloGene ID", idType)];
  tmp = reftable[(idx2&idx3), c("HomoloGene ID", "Symbol", "Entrez")]
  colnames(tmp1)[2] = "Mouse"
  res = merge(tmp1, tmp, by = 1)[,-1]

  require(biomaRt)
  attributes <- filters <- "mgi_symbol"
  if(tolower(idType) == "entrez") attributes <- filters <- "entrezgene_id"
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = attributes, mart = mouse,
                   filters = filters, values = x,
                   attributesL = c("hgnc_symbol", "entrezgene_id"), martL = human)
  colnames(genesV2) = c("Mouse", "Symbol", "Entrez")
  genesV2$Symbol[genesV2$Symbol==""] = TransGeneID(genesV2$Entrez[genesV2$Symbol==""], "Entrez", "Symbol")
  genesV2 = genesV2[!is.na(genesV2$Symbol), ]
  res = rbind.data.frame(res, genesV2)

  y = setdiff(x, res[,1])
  y_entrez = TransGeneID(toupper(y), "Symbol", "Entrez")
  y_symbol = TransGeneID(y_entrez, "Entrez", "Symbol")
  names(y_entrez) = names(y_symbol) = NULL
  tmp = data.frame(Mouse = y, Symbol = y_symbol, Entrez = y_entrez, stringsAsFactors = FALSE)
  res = rbind.data.frame(res, tmp)

  idx = duplicated(paste0(res[,1], "_", res[,2], "_", res[,3])) | is.na(res[,2])
  res = res[!idx, ]
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
