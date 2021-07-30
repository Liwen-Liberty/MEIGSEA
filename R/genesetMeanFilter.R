#' Refine immunophenotype signature
#' @param  exp.profile Expression profile with row as genes, column as samples.
#' @param  signature.df A data frame containing immunophenotype signature
#' @param  cor.cutoff Setting the cutoff of correlation coefficient when refining the immunophenotype
#'
#' @return Refined signature infomation

genesetMeanFilter <- function(exp.profile, signature.df, cor.method = "spearman", cor.cutoff = 0.5) {
  #######
  # (1) Delete genes that are not in the expression profile --------------------
  #'
  #######
  idx <- match(signature.df$Gene, rownames(exp.profile))
  rm.genes <- signature.df$Gene[which(is.na(idx))]
  if (length(rm.genes) != 0) signature.df <- signature.df[-which(is.na(idx)), ]
  if (dim(signature.df)[1] == 0) stop("All signature genes are not in the expression profile")

  #######
  # (2) Calculate the correlation between expression of each gene and the mean expression value of the gene set
  #'
  #######
  cor.res <- tapply(signature.df$Gene, signature.df$setAnno, function(geneset) {
    signature.exp <- exp.profile[geneset, ]
    geneset.MeanCor <- genesetMeanCor(geneset.exp = signature.exp, cor.method = cor.method, cor.cutoff = cor.cutoff)
    return(geneset.MeanCor)
  })

  cor.res.df <- frameInlist2addOneFrame(cor.res, "setAnno")
  cor.res.df <- cor.res.df[which(cor.res.df$signif == "sig"), ]

  #######
  # (3) Store Result --------------------
  #'
  #######
  # 1）Ignore the direction, and only significant genes are retained
  signatureGene.list <- split(cor.res.df$gene, cor.res.df$setAnno)

  # 2）Distinguish the direction, each signature has significantly positive and negative genes divided into two list
  cor.res.df$label <- paste(cor.res.df$setAnno, cor.res.df$direction, sep = "_")

  signatureGene.list <- list(df = cor.res.df, list = signatureGene.list)
  return(signatureGene.list)
}
