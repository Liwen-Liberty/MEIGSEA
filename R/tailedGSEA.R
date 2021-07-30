#' one/two-tail GSEA strategy
#'
#' @param  exp.profile RNA-seq data in raw count or microarray expression profile
#' with row as genes, column as samples.
#' @param  is.rawcount Whether or not the 'raw.count.profile' is RNA-seq raw count.
#' For microarray expression data as input, is.rawcount = FALSE.
#' @param  sample.group Expression profile sample column label, used to calculate the
#' differentially expressed log2FC between samples
#' @param  case/control In default, case = "1", control = "0", or the specific character
#' representing mutant (case) and wild type (control) according to 'sample.group'.
#' @param  sig.df Multiple signature gene sets to be enriched.
#' @param  min.sz At least how many immunophenotype genes are matched in the expression
#' profile, min.sz = 15 by default.
#' @param  pvalueCutoff  The threshold of GSEA p value, if set as 1, all enrichment
#' results will be output
#' @param  pAdjustMethod Method used to adjust p values, "BH" by default.
#' @param  is.fgsea Whether to call fgsea function directly, default is FALSE.
#' @param  is.fisherPcombine Whether to integrate p-values from two-tail GSEA,default is FALSE
#'
#' @return Enrichment results of two GSEA strategies, and log2FC


tailedGSEA <- function(exp.profile, is.rawcount = TRUE, sample.group, case, control, sig.df, min.sz = 15, pvalueCutoff = 1, pAdjustMethod = "BH",
                       is.fgsea = FALSE, is.fisherPcombine = FALSE) {

  ###
  #' 1. Calculate log2FC  of gene expression differences between MUT and WT samples
  #'
  ################################
  if (is.rawcount) {
    #  Check whether the expression profile is raw count (non-negative integer)
    if (any(!is.integer(as.matrix(exp.profile)))) {
      stop("** This is not a integer profile ! ")
    }
    if (any(exp.profile < 0)) {
      stop("** This is not a non-negative profile ! ")
    }

    ## ! Calculate Shrunken log2 foldchanges (LFC) based on DESeq2
    log2FC <- get.shrunken.log2FC(raw.count.profile = exp.profile, sample.group = sample.group, case = case, control = control)
  } else {
    # Calculate the mean of sample expression
    group.mean <- tapply(1:ncol(exp.profile), sample.group, function(sub.exp) {
      tmp.mean <- rowMeans(exp.profile[, sub.exp])
      return(tmp.mean)
    })
    # Avoid log2(0) or 0 in the denominator
    group.mean[[case]][which(group.mean[[case]] == 0)] <- 1e-05
    group.mean[[control]][which(group.mean[[control]] == 0)] <- 1e-05

    log2FC <- log2(group.mean[[case]] / group.mean[[control]])
  }



  ###
  #' 2. One-tail GSEA
  #'
  ################################
  # Delete gene sets with size<min.sz
  tmp.sig.df <- sig.df[, c("setAnno", "gene")]
  geneset.size <- table(tmp.sig.df[, "setAnno"])
  if (any(geneset.size < min.sz)) {
    tmp.sig.df <- tmp.sig.df[-which(tmp.sig.df[, "setAnno"] %in% names(geneset.size)[which(geneset.size < min.sz)]), ]
  }

  if (nrow(tmp.sig.df) > 0) {
    # Note: the gene sets larger than 500 are filtered out
    onetailedGSEA.res <- onetailedGSEA(log2FC = log2FC, sig.df = tmp.sig.df, min.sz = min.sz, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, is.fgsea = is.fgsea)
    cat("*** onetail GSEAres finished ! \n")
  } else {
    onetail.GSEAres <- NULL
    cat("*** onetail GSEAres filed due to null signature set ! \n")
  }

  ###
  #' 3. Two-tail GSEA
  #'
  ################################
  # Delete gene sets with size<min.sz
  tmp.sig.df <- sig.df[, c("label", "gene")]
  geneset.size <- table(tmp.sig.df[, "label"])
  if (any(geneset.size < min.sz)) {
    tmp.sig.df <- tmp.sig.df[-which(tmp.sig.df[, "label"] %in% names(geneset.size)[which(geneset.size < min.sz)]), ]
  }

  if (nrow(tmp.sig.df) > 0) {
    twotailedGSEA.res <- twotailedGSEA(
      log2FC = log2FC, sig.df = tmp.sig.df, min.sz = min.sz, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, is.fgsea = is.fgsea,
      is.fisherPcombine = TRUE
    )
    cat("*** twotail GSEAres finished ! \n")
  } else {
    twotailedGSEA.res <- NULL
    cat("*** twotail GSEAres filed due to null signature set ! \n")
  }


  # return results
  GSEAresult <- list(log2FC = log2FC, oneTailRes = onetailedGSEA.res, twoTailRes = twotailedGSEA.res)
  return(GSEAresult)
}
