#' Obtain co-mutated gene pairs for conditional association analysis
#'
#' @param mutationCalling.map A matrix of gene mutation status (0 or 1, in character),
#' row as genes, column as samples.
#' @param sig.corscore.df Significant associations between mutations and immunophenotypes
#'
#' @return co-mutated gene pairs

getRemainGenepairs <- function(mutationCalling.map, sig.corscore.df) {

  ## 1）Obtain co-mutated gene pairs
  fisher.res <- get_mut_inter_fisher(mutationMat = mutationCalling.map, outfile = NULL, min_count = 9)
  sig.genepairs <- strsplit(rownames(fisher.res)[which(fisher.res$pvalue < 0.05)], split = "_")

  ## 2）Screen gene pairs satisfying the conditional association analysis
  remain.genepairs.list <- lapply(sig.genepairs, function(pair) {
    remain.pair <- data.frame(condition = NA, test = NA)

    condition.mutationCalling.map <- mutationCalling.map[pair[2], which(mutationCalling.map[pair[1], ] == 0)]
    if (length(which(condition.mutationCalling.map == 1)) >= 10) remain.pair <- rbind.data.frame(remain.pair, c(pair[1], pair[2]))

    condition.mutationCalling.map <- mutationCalling.map[pair[1], which(mutationCalling.map[pair[2], ] == 0)]
    if (length(which(condition.mutationCalling.map == 1)) >= 10) remain.pair <- rbind.data.frame(remain.pair, c(pair[2], pair[1]))

    remain.pair <- na.omit(remain.pair)
    return(remain.pair)
  })
  remain.genepairs.list <- remain.genepairs.list[which(sapply(remain.genepairs.list, nrow) == 2)]
  remain.genepairs <- do.call(rbind.data.frame, remain.genepairs.list)

  ## 3）Screen gene pairs whose gene mutations both associated with a specific  immunophenotype
  immunophenotypes <- unique(sig.corscore.df$signature)
  remain.genepairs <- do.call(rbind, lapply(immunophenotypes, function(x) {
    immunr.related.genes <- unique(sig.corscore.df[which(sig.corscore.df$signature == x), "driver"])
    remain.genepairs[(remain.genepairs$condition %in% immunr.related.genes) & (remain.genepairs$test %in% immunr.related.genes), ]
  }))
  dup.rows <- which(duplicated(remain.genepairs))
  if (length(dup.rows) > 0) remain.genepairs <- remain.genepairs[-dup.rows, ]
  if (nrow(remain.genepairs) == 0) remain.genepairs <- NULL

  return(remain.genepairs)
}
