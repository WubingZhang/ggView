##' Filtering a count matrix to remove low read genes
##'
##' This function takes a count matrix and removes genes that don't reach a minimum number of reads in a certain number of samples.
##' @param m raw count matrix (genes are rows, columns are samples)
##' @param minSamples numeric; minimum number of samples in which CPM threshold has to be reached; defaults to a tenth of available samples
##' @param minCpm numeric; the minimum counts per million (CPM) that has to be reached
##' @return filtered count matrix
##' @author Dorothee Nickles
##' @import edgeR
##' @export
filterNvoom <- function(m, minSamples=ncol(m)/10, minCpm=0.25) {
  tmp <- edgeR::cpm(m)
  sel <- apply(tmp, 1, function(x) {
    sum(x >= minCpm) >= minSamples
  })
  return(m[sel,])
}

#' Filter and normalize a count matrix
#' This wrapper function combines filtering out genes with low reads
#' in a number of samples with normalization.
#' @docType methods
#' @name TransformCount
#' @rdname TransformCount
#'
#' @param m raw count matrix (genes are rows, columns are samples)
#' @param minS numeric; minimum number of samples in which count threshold has to be
#' reached; defaults to a tenth of available samples
#' @param method the normalization method, one of "none", "TPM", "voom", "vst", "rlog".
#' @param idType same as in 'Count2TPM', only required when method = "TPM".
#' @param org same as in 'Count2TPM', only required when method = "TPM".
#'
#' @return filtered and normalized count matrix
#' @author Wubing Zhang
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom
#' @importFrom DESeq2 varianceStabilizingTransformation rlog
#' @export

TransformCount <- function(m, minS=ncol(m)/10, method=c("TPM", "voom", "vst", "rlog")[1],
                           idType = "Ensembl", org="hsa") {
  # pre-filtering
  sel <- rowSums(m >= 0) >= minS
  m_sel <- m[sel, ]

  # Count data transformations
  if (method == "none"){
    v = m_sel
  }else if (method == "TPM") { v <- Count2TPM(m_sel, idType=idType, org=org)
  } else if(method == "voom") {
    v <- filterNvoom(m, minSamples=minS, minCpm=0.25)
    # dge <- edgeR::DGEList(m_sel)
    # dge <- edgeR::calcNormFactors(dge)
    # v <- limma::voom(dge)$E
  } else if(method == "vst"){
    v <- DESeq2::varianceStabilizingTransformation(m_sel, blind = FALSE)
  } else if(method == "rlog"){
    v <- DESeq2::rlog(m_sel, blind = FALSE)
  }else{
    stop("Undefined normalization method !")
  }
  return(v)
}


#' Convert read counts to transcripts per million (TPM)
#'
#' @docType methods
#' @name Count2TPM
#' @rdname Count2TPM
#'
#' @param countMat A read count matrix, with geneid as rownames and sample as columns.
#' @param idType Type of gene id.
#' @param org Organism, hsa or mmu.
#'
#' @return A tpm expression profile.
#'
#' @author Wubing Zhang
#' @import biomaRt
#' @export
#'
Count2TPM <- function(countMat, idType = "Ensembl", org="hsa")
{
  requireNamespace("biomaRt")
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  type = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "start_position", "end_position")
  if(org=="mmu") type[3] = "mgi_symbol"
  # listEnsemblArchives()
  # listMarts()
  # listAttributes()
  ds = datasets[grepl(org, datasets)]
  mart <- useMart(host = "www.ensembl.org", biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
  ensembl = getBM(attributes=type, mart = mart)
  ensembl$Length <- abs(ensembl$end_position - ensembl$start_position)
  if(toupper(idType) == "ENSEMBL"){
    len <- ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), "Length"]
    rownames(countMat) = ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), 3]
  }
  else if(toupper(idType) == "SYMBOL")
    len <- ensembl[match(rownames(countMat), ensembl[,3]), "Length"]
  else if(toupper(idType) == "ENTREZ")
    len <- ensembl[match(rownames(countMat), ensembl[,2]), "Length"]
  else
    stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")

  na_idx = which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
    countMat = countMat[!is.na(len),]
    len = len[!is.na(len)]
  }
  tmp <- countMat / len
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
  TPM = TPM[!duplicated(rownames(TPM)),]
  return(TPM)
}
