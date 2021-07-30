#' Calculate Shrunken log2 foldchanges (LFC) based on DESeq2
#'
#' @param  raw.count.profile RNA-seq data in raw count with row as genes, column as samples
#' @param  sample.group Expression profile sample column label, used to calculate the
#' differentially expressed log2FC between samples
#' @param  case/control In default, case = "1", control = "0", or the specific character
#' representing mutant (case) and wild type (control) according to 'sample.group'.
#'
#' @return

get.shrunken.log2FC <- function(raw.count.profile, sample.group, case, control) {

  ###
  #' 1. Determine whether the expression profile is raw count (non-negative integer)
  #'
  ################################
  if (any(!is.integer(as.matrix(raw.count.profile)))) {
    stop("** This is not a integer profile ! ")
  }
  if (any(raw.count.profile < 0)) {
    stop("** This is not a non-negative profile ! ")
  }

  ###
  #' 2. Calculate log2FC  of gene expression differences between MUT and WT samples
  #'
  ################################


  # object construction
  sample.group <- factor(sample.group, levels = c(control, case))
  dds <- DESeq2::DESeqDataSetFromMatrix(raw.count.profile, DataFrame(sample.group), ~sample.group)

  # standard analysis
  dds <- DESeq2::DESeq(dds, betaPrior = TRUE, parallel = TRUE)
  res <- DESeq2::results(dds)

  # moderated log2 fold changes
  shrink.log2FC <- res$log2FoldChange
  names(shrink.log2FC) <- rownames(res)

  # drop NA
  shrink.log2FC <- na.omit(shrink.log2FC)
  return(shrink.log2FC)
}
