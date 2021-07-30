#' Normalize counts to gene length
#'
#' @param counts Counts data
#' @param lengths Genes' lengths. If null use a precalculated gene lengths
#'
#' @return TPM values
CountToTPM <- function(counts, gene.length = NULL) {
  id.type <- autoDetectGeneIdType(rownames(counts)[1])

  if (is.null(gene_lengths)) {
    # data("gene_lengths")
    gene_lengths <- MEIGASE::gene_lengths
    gene_lengths <- gene_lengths[which(!is.na(gene_lengths[, id.type])), ]
    gene.length <- gene_lengths[, "lengths"]
    names(gene.length) <- tolower(gene_lengths[, id.type])
  }

  names(gene.length) <- tolower(names(gene.length))
  A <- intersect(tolower(rownames(counts)), names(gene.length))
  if (length(A) == 0) {
    stop("When gene_length is NULL, only HGNC symbol, ENTREZ or ENSEMBL gene IDs are supported for counts data.")
  } else {
    counts <- counts[match(A, tolower(rownames(counts))), ]
    gene.length <- gene.length[A]

    rate <- counts / lengths
    TPM <- apply(rate, 2, function(x) 1e6 * x / sum(x))
  }

  return(TPM)
}

#' Detect gene id type
#'
#' @param id Gene id in counts data
#'
#' @return  Gene id Type
autoDetectGeneIdType <- function(id) {
  type <- "HGNC_symbol"
  if (grepl("^[Ee][Nn][Ss][A-Za-z]{0,3}[Gg][0-9]+", id)) {
    type <- "Ensembl_ID"
  } else if (grepl("^[0-9]+$", id)) type <- "Entrez_ID"
  return(type)
}
