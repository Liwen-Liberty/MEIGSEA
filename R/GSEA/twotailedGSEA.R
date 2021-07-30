#' two-tailed GSEA strategy
#'
#' @param  log2FC The pre-calculated log2FC is used for gene rank
#' @param  sig.df Multiple signature gene sets to be enriched.
#' @param  min.sz At least how many immunophenotype genes are matched in the expression
#' profile, min.sz = 15 by default.
#' @param  pvalueCutoff The threshold of GSEA p value, if set as 1, all enrichment
#' results will be output
#' @param  pAdjustMethod Method used to adjust p values, "BH" by default.
#' @param  is.fgsea Whether to call fgsea function directly, default is FALSE.
#' @param  is.fisherPcombine Whether to integrate p-values from two-tail GSEA,default is TRUE
#'
#' @return GSEA result

twotailedGSEA <- function(log2FC, sig.df, min.sz = 15, pvalueCutoff = 1, pAdjustMethod = "BH", is.fgsea = FALSE, is.fisherPcombine = TRUE) {

  ###
  # 1. Rank genes based on log2FC  --------------------
  #'
  #' ###############################
  log2FC <- log2FC[order(log2FC, decreasing = TRUE)]


  ###
  # 2. Perform GSEA using clusterProfiler --------------------
  #'
  #' ###############################


  sig.df <- unique(sig.df)

  ### Perform GSEA
  if (is.fgsea) {
    sig.list <- split(sig.df[, "gene"], sig.df[, "label"])
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
  cat("*** GSEA calculation finished ! \n")

  ###
  # 3. Intergrate p values enriched by positive and negative gene sets  --------------------
  #'
  #' ###############################

  #' ######
  # (1) Identify gene sets that have both positively and negatively subsets  --------------------
  #' ######
  sig <- table(gsub("_.{4}tive", "", GSEAres$ID))
  bidirection.geneset <- names(sig)[sig > 1]

  #' ######
  # (2) Determine if p values need to be integrated --------------------
  #' ######

  if (is.fisherPcombine & length(bidirection.geneset) > 0) {
    twoTailRes <- twotailedGSEA.res.comb(GSEAres)
  } else {
    twoTailRes <- as.data.frame(GSEAres)[c("ID", "NES", "pvalue")]
    twoTailRes$ID <- gsub("_.{4}tive", "", GSEAres$ID)
  }

  GSEAresult <- list(GSEAres = GSEAres, twoTailRes = twoTailRes)
  return(GSEAresult)
}
