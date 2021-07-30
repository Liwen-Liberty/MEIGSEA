
#' Calculate the correlation between expression of each gene and the mean expression
#' level of the gene set
#'
#' @param  geneset.exp Expression profile with row as genes, column as samples.
#' @param  cor.method Method used to compute the correlation, either "spearman",
#' "pearson" or "kendall", "spearman" by default
#' @param  padjust.method Method used to adjust p values, "BH" by default.
#' @param  cor.cutoff Setting the cutoff of correlation coefficient when
#' refining the immunophenotype
#'
#' @return A data frame including gene, cor and p.adjust

genesetMeanCor <- function(geneset.exp, cor.method = "spearman", padjust.method = "BH", cor.cutoff = NULL) {
  # Calculate the mean expression level of the gene set
  geneset.mean <- colMeans(geneset.exp, na.rm = TRUE)

  # Calculate the correlation between expression of each gene and the mean
  # expression level of the gene set
  cor.withMean <- apply(as.matrix(geneset.exp), 1, function(x) {
    cor.res <- cor.test(x, geneset.mean, method = cor.method, adjust = padjust.method)
    return(data.frame(cor = signif(cor.res$estimate, 3), p.ad = signif(cor.res$p.value, 3)))
  })

  cor.withMean <- do.call(rbind.data.frame, cor.withMean)
  cor.withMean$gene <- rownames(geneset.exp)
  # significant level
  cor.withMean$signif <- "non-sig"
  cor.withMean$signif[which(cor.withMean$p.ad < 0.05)] <- "sig"
  # correlation direction
  cor.withMean$direction <- "positive"
  cor.withMean$direction[which(cor.withMean$cor < 0)] <- "negative"

  # return filtered result
  if (!is.null(cor.cutoff)) cor.withMean <- cor.withMean[which(abs(cor.withMean$cor) >= cor.cutoff), ]

  return(cor.withMean)
}
