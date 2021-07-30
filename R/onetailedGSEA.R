#'  One tailed GSEA
#' @param  log2FC The pre-calculated log2FC is used for gene rank
#' @param  sig.df Multiple signature gene sets to be enriched.
#' @param  min.sz At least how many immunophenotype genes are matched in the expression
#' profile, min.sz = 15 by default.
#' @param  pvalueCutoff The threshold of GSEA p value, if set as 1, all enrichment
#' results will be output
#' @param  pAdjustMethod Method used to adjust p values, "BH" by default.
#' @param  is.fgsea Whether to call fgsea function directly, default is FALSE.
#'
#' @return GSEA result

onetailedGSEA <- function(log2FC, sig.df, min.sz = 15, pvalueCutoff = 1, pAdjustMethod = "BH", is.fgsea = FALSE) {

  ###
  # 1.Rank genes based on abs(log2FC)  --------------------
  #'
  #' ###############################
  log2FC <- abs(log2FC)
  log2FC <- log2FC[order(log2FC, decreasing = TRUE)]


  ###
  #' 2. Perform GSEA using clusterProfiler
  #'
  ################################


  sig.df <- unique(sig.df)

  ### Perform GSEA
  if (is.fgsea) {
    sig.list <- split(sig.df[, "variable"], sig.df[, "setAnno"])
    GSEAres <- fgsea::fgsea(
      pathways = sig.list,
      stats = log2FC,
      minSize = min.sz
    )
    colnames(GSEAres) <- gsub("pathway", "ID", colnames(GSEAres))
    colnames(GSEAres) <- gsub("pval", "pvalue", colnames(GSEAres))
    colnames(GSEAres) <- gsub("padj", "p.adjust", colnames(GSEAres))
    colnames(GSEAres) <- gsub("^ES$", "enrichmentScore", colnames(GSEAres))
    colnames(GSEAres) <- gsub("size", "setSize", colnames(GSEAres))
    colnames(GSEAres) <- gsub("leadingEdge", "leading_edge", colnames(GSEAres))
  } else {
    GSEAres <- clusterProfiler::GSEA(log2FC,
      minGSSize = min.sz,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      TERM2GENE = sig.df
    )
  }

  return(GSEAres)
}
