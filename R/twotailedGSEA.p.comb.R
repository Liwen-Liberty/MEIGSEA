#'  Intergrate p values enriched by positive and negative gene sets using Fisher’s method
#' @param  GSEA.res GSEA result from two-tail GSEA strategy
#'
#' @return Intergrated p value and GSEA results of positive and negative gene sets


twotailedGSEA.res.comb <- function(GSEA.res) {


  ###
  # 1. Identify gene sets that have both positively and negatively subsets --------------------
  #'
  #' ###############################

  sig <- table(gsub("_.{4}tive", "", GSEA.res$ID))
  bidirection.geneset <- names(sig)[sig > 1]



  ###
  # 2. Intergrate p values using Fisher’s method --------------------
  #'
  #' ###############################

  if (length(bidirection.geneset) > 0) {
    #######
    # (1) Construct a matrix including p values enriched by positive and negative gene sets
    #######
    p.matrix <- sapply(bidirection.geneset, function(sig) {
      tmp.2p <- GSEA.res$pvalue[which(GSEA.res$ID %in% paste(sig, c("positive", "negative"), sep = "_"))]
      return(tmp.2p)
    })
    p.matrix <- t(p.matrix)

    #  Intergrate p values
    fisher.p <- fisher.method(p.matrix,
      p.corr = "BH",
      zero.sub = 1e-05,
      na.rm = TRUE
    )
    rownames(fisher.p) <- bidirection.geneset



    ###
    # 3. Refresh the enrichment  results--------------------
    #'
    #' ###############################

    # Extract the enrichment result of other signatures that are either positive or negative
    twoTailRes <- as.data.frame(GSEA.res)[which(gsub("_.{4}tive", "", GSEA.res$ID) %in% names(sig)[sig == 1]), c("ID", "NES", "pvalue")]

    # For nagative subsets, change the direction for NES
    twoTailRes$NES[grep("negative", twoTailRes$ID)] <- -1 * twoTailRes$NES[grep("negative", twoTailRes$ID)]

    # Remove the direction of ID
    twoTailRes$ID <- gsub("_.{4}tive", "", twoTailRes$ID)

    #######
    # Set the direction of intergrated enrichment score NES
    #'
    #######

    comb.NES.matrix <- sapply(bidirection.geneset, function(sig) {
      tmp.2res <- GSEA.res[which(GSEA.res$ID %in% paste(sig, c("positive", "negative"), sep = "_")), c("NES", "pvalue")]
      if (all(tmp.2res$pvalue > 0.05)) {
        # Neither positive nor negative subsets were significant, set NES=NA
        comb.NES <- NA
      } else if (any(tmp.2res$pvalue > 0.05)) {
        # Any one of the subsets is significant, set NES=significant NES
        comb.NES <- tmp.2res$NES[which(tmp.2res$pvalue < 0.05)]
      } else {
        # Both the subsets were significant

        if (purrr::reduce(tmp.2res$NES, `*`) > 0) {
          # NESs with different sign, set NES = positive subsets' NES
          comb.NES <- Inf * sign(tmp.2res$NES[1])
        } else {
          # NESs with same sign, set NES = Inf
          comb.NES <- GSEA.res$NES[which(GSEA.res$ID == paste0(sig, "_positive"))]
        }
      }

      return(comb.NES)
    })

    twoTailRes <- rbind.data.frame(twoTailRes, data.frame(
      ID = rownames(fisher.p),
      NES = comb.NES.matrix,
      pvalue = fisher.p$"p.value"
    ))
  } else {
    twoTailRes <- NULL
    cat("No bi-directions signature ! \n")
  }
  return(twoTailRes)
}
