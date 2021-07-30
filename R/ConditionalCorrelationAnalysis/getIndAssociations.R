#' Obtain Indepedent associations
#'
#' @param remain.genepairs Co-mutated gene pairs for conditional association analysis
#' @param sig.corscore.df Significant associations between mutations and immunophenotypes
#' @param raw.count.profile Expression profile with row as genes, column as samples.
#' @param mutationCalling.map A matrix of gene mutation status (0 or 1, in character), row as genes, column as samples.
#' @param TPM.profile Given tpm expression profile.
#' @param sig.cutoff Cutoff for significant associations, sig.cutoff = 0.01 by default.
#' See the MEIGSEA function for details on the remaining parameters
#'
#' @return Indepedent associations


getIndAssociations <- function(remain.genepairs, sig.corscore.df, raw.count.profile, mutationCalling.map, TPM.profile, sig.cutoff = 0.01, signature.df,
                               case = "1", control = "0", cor.cutoff = NULL, cor.method = "spearman", min.sz = 15, perm.times = 100, pAdjustMethod = "BH") {

  # Conditional correlation analysis
  ConditionalRelation <- list()
  for (i in 1:nrow(remain.genepairs)) {
    genepair <- remain.genepairs[i, ]
    IntegratedGSEA.res <- getConditionalRelation(genepair,
      mutationCalling.map = mutationCalling.map, raw.count.profile = raw.count.profile, tpm.exp = TPM.profile, signature.df = signature.df,
      case = case, control = control, cor.cutoff = cor.cutoff, cor.method = cor.method, min.sz = min.sz, perm.times = perm.times, pAdjustMethod = pAdjustMethod
    )
    ConditionalRelation <- c(ConditionalRelation, list(IntegratedGSEA.res))
  }
  ConditionalRelation <- do.call(rbind, ConditionalRelation)

  # Add two columns indicating the condition and test genes
  ConditionalRelation$condition <- unlist(lapply(strsplit(ConditionalRelation$driver, "_"), function(x) x[1]))
  ConditionalRelation$test <- unlist(lapply(strsplit(ConditionalRelation$driver, "_"), function(x) x[2]))
  # Recalculated the adjusted p value
  ConditionalRelation$p.adj <- p.adjust(ConditionalRelation$p.value, method = pAdjustMethod)

  # Convert to list format
  immune.types <- unique(ConditionalRelation$signature)
  ConditionalRelation <- lapply(immune.types, function(x) {
    sig.drivers <- sig.corscore.df[which(sig.corscore.df$signature == x), "driver"]
    per.ConditionalRelation <- ConditionalRelation[which(ConditionalRelation$signature == x), ]
    subset(per.ConditionalRelation, condition %in% sig.drivers & test %in% sig.drivers)
  })
  names(ConditionalRelation) <- immune.types

  # Obtain Indepedent associations
  sig.corscore.df <- do.call(rbind, lapply(immune.types, function(i) {
    per.corscore.df <- sig.corscore.df[which(sig.corscore.df$signature == i), ]
    per.ConditionalRelation <- ConditionalRelation[[i]]
    # If p value is not significant, the association between test gene mutation and the immunophenotype is affected by the condition gene mutation
    affected.test.genes <- unique(per.ConditionalRelation[which(per.ConditionalRelation$p.adj > sig.cutoff | per.ConditionalRelation$weight == 0), "test"])
    if (length(affected.test.genes) > 0) {
      # condition genes affecting test gene
      condition.genes <- lapply(affected.test.genes, function(test) {
        paste(per.ConditionalRelation[which((per.ConditionalRelation$p.adj > sig.cutoff | per.ConditionalRelation$weight == 0) & per.ConditionalRelation$test == test), "condition"], collapse = ",")
      })
      names(condition.genes) <- affected.test.genes
      per.corscore.df$FILTER[which(per.corscore.df$driver %in% affected.test.genes)] <- "Dependent"
      for (test.genes in affected.test.genes) per.corscore.df$CoMutGenes[which(per.corscore.df$driver == test.genes)] <- condition.genes[[test.genes]]
    }
    per.corscore.df
  }))

  return(sig.corscore.df)
}


#' Conditional correlation analysis
#'
#' @param genepair a specific co-mutated gene pair
#' @param raw.count.profile expression profile with row as genes, column as samples.
#' @param mutationCalling.map A matrix of gene mutation status (0 or 1, in character), row as genes, column as samples.
#' @param TPM.profile Given tpm expression profile.
#' See the MEIGSEA function for details on the remaining parameters
#'
#' @return Indepedent associations

getConditionalRelation <- function(genepair, raw.count.profile, mutationCalling.map, tpm.exp, signature.df, case = "1", control = "0",
                                   cor.cutoff = NULL, cor.method = "spearman", min.sz = 15, perm.times = 100, pAdjustMethod = "BH") {
  condition.gene <- genepair[[1]]
  test.gene <- genepair[[2]]
  ## 1)  Extracte samples whose condition gene mutation state is 0
  sample.idx <- which(mutationCalling.map[condition.gene, ] == "0") # 样本标签
  # mutation data
  sample.group <- matrix(mutationCalling.map[test.gene, sample.idx], nrow = 1)
  dimnames(sample.group) <- list(paste0(genepair, collapse = "_"), colnames(mutationCalling.map)[sample.idx])
  # expression data
  raw.count.profile <- raw.count.profile[, sample.idx]
  tpm.exp <- tpm.exp[, sample.idx]

  ## 2) Perform Integrated GSEA strategy to get associations between test gene mutation and all immunophenotypes
  IntegratedGSEA.res <- IntegratedGSEA(
    raw.count.profile = raw.count.profile, tpm.exp = tpm.exp, mutationCalling.map = sample.group, signature.df = signature.df,
    case = case, control = control, cor.cutoff = cor.cutoff, cor.method = cor.method, min.sz = min.sz, perm.times = perm.times, pAdjustMethod = pAdjustMethod
  )

  corscore.df <- IntegratedGSEA.res$corscore.df

  return(corscore.df)
}
