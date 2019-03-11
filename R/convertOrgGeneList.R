# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                   values = x , mart = mouse, attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  genesV2 = genesV2[!duplicated(genesV2$MGI.symbol), ]
  rownames(genesV2) = genesV2$MGI.symbol
  
  return(genesV2[x,"HGNC.symbol"])
}

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                   values = x , mart = human, attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  
  genesV2 = genesV2[!duplicated(genesV2$HGNC.symbol), ]
  rownames(genesV2) = genesV2$HGNC.symbol
  
  return(genesV2[x,"MGI.symbol"])
}
