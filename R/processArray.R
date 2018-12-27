#' Method to preprocess Microarray dataset from GEO
#'
#' @docType methods
#' @name processArray
#' @rdname processArray
#'
#' @param expr A file path to GSE_series_matrix.
#' @param GPL A file path to GPL soft file.
#' @param symbol The colname of gene symbol in `GPL`.
#' @param org 'hsa' or 'mmu'.
#' @return A normalized expression profile.
#'
#' @author Wubing Zhang
#' @importFrom GEOquery getGEO Table
#' @importFrom data.table fread
#' @importFrom limma normalizeQuantiles
#' @export

processArray <- function(expr="GSE5821_series_matrix.txt",
                         GPL="GPL96.soft", symbol=NA, org = "hsa"){
  #===Determine the parameters using shell command
  Sample_title = system(paste0("grep -n '!Sample_title' ", expr), intern = TRUE)
  Sample_title = unlist(strsplit(Sample_title, "\t\""))
  Sample_title = gsub("\"", "", Sample_title)
  Sample_title[1] = "Gene"
  skip = system(paste0("grep -n '!series_matrix_table_begin' ", expr), intern = TRUE)
  skip = as.integer(gsub(":!series_matrix_table_begin", "", skip))
  data_process = system(paste0("grep '!Sample_data_processing' ", expr), intern = TRUE)
  if(length(data_process)==0) data_process = ""
  data_process = paste0(data_process, collapse = " ")
  quantile = TRUE
  if(grepl("quantile|RMA", data_process))  quantile = FALSE

  #===Read expression data and map probes to genes based on GPL annotation
  expr = data.table::fread(expr, sep = "\t", header = TRUE,
               stringsAsFactors = FALSE, skip = skip, fill = TRUE)
  expr = as.data.frame(expr, stringAsFactors=FALSE)
  colnames(expr) = Sample_title
  if(grepl("soft|annot", GPL)){
    gpl <- GEOquery::Table(GEOquery::getGEO(filename = GPL))
  }else{
    gpl <- data.table::fread(GPL, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  gpl = as.data.frame(gpl, stringAsFactors=FALSE)
  idx = is.na(gpl$ID) | gpl$ID==""
  gpl = gpl[!idx, ]
  rownames(gpl) = gpl$ID
  if(is.na(symbol)) symbol = which(grepl("symbol", colnames(gpl), ignore.case = TRUE))[1]
  expr[,1] = gpl[as.character(expr[,1]), symbol]
  if(symbol=="gene_assignment"){
    expr[,1] = sapply(expr[,1], function(x){
      tmp = unlist(strsplit(x, " // "))
      if(length(tmp)>1) return(tmp[2])
      else return(NA)
    })
  }
  #=====Convert gene symbols to official symbols
  expr[,1] = gsub("///.*", "", expr[,1])
  if(!is.na(TransGeneID(expr[1,1], "Symbol", "Symbol", org))){
    require(MAGeCKFlute)
    expr[,1] = TransGeneID(expr[,1], "Symbol", "Symbol", org)
  }
  idx = duplicated(expr[,1]) | is.na(expr[,1])
  expr = expr[!idx, ]
  rownames(expr) = expr[,1]
  expr = expr[,-1]
  expr = matrix(as.numeric(as.matrix(expr)), nrow = nrow(expr),
               dimnames = list(rownames(expr), colnames(expr)))
  #======Normalize the expression data
  if(quantile) expr = limma::normalizeQuantiles(expr)
  expr = as.matrix(expr)
  if(max(expr, na.rm = TRUE)>50) expr = log2(expr+0.0000001)
  expr = as.data.frame(expr, stringAsFactors=FALSE, check.names=FALSE)
  return(expr)
}

