#' Conditional Correlation Analysis
#'
#' @param corscore.df Association between mutations and immunophenotypes calculated by
#' integrated GSEA algorithms
#' @param raw.count.profile Expression profile with row as genes, column as samples.
#' @param mutationCalling.map A matrix of gene mutation status (0 or 1, in character),
#' row as genes, column as samples.
#' @param tpm.exp Given tpm expression profile.
#' @param sig.cutoff Cutoff for significant associations,sig.cutoff = 0.01 by default.
#' See the MEIGSEA function for details on the remaining parameters
#'
#' @return Significant correlation pairs


CondCorrAnal <- function(corscore.df, raw.count.profile, is.rawcount = TRUE, gene.length = NULL, tpm.exp = NULL, mutationCalling.map,
                         sig.cutoff = 0.01, signature.df, case = "1", control = "0", cor.cutoff = NULL, cor.method = "spearman", min.sz = 15, perm.times = 100, pAdjustMethod = "BH") {

  # Convert mutation data into numeric form
  mutationCalling.map <- as.matrix(mutationCalling.map)
  mode(mutationCalling.map) <- "numeric"

  # Get significantly associations
  sig.corscore.df <- corscore.df[which(corscore.df$p.adj < sig.cutoff & corscore.df$weight > 0), ]
  sig.corscore.df$FILTER <- "PASS"
  sig.corscore.df$CoMutGenes <- NA

  # Obtain co-mutated gene pairs for conditional association analysis
  remain.genepairs <- getRemainGenepairs(mutationCalling.map, sig.corscore.df)
  if (!is.null(remain.genepairs)) {
    # raw count to TPM if needed
    if (is.null(tpm.exp) & is.rawcount) {
      TPM.profile <- CountToTPM(raw.count.profile, gene.length = NULL)
      cat("** TPM transformed finished ! \n")
    } else if (is.null(tpm.exp) & (!is.rawcount)) {
      TPM.profile <- raw.count.profile
    } else {
      TPM.profile <- tpm.exp
    }
    # Get independent associations
    sig.corscore.df <- getIndAssociations(remain.genepairs, sig.corscore.df, raw.count.profile, mutationCalling.map,
      TPM.profile = TPM.profile, sig.cutoff = sig.cutoff, signature.df = signature.df,
      case = case, control = control, cor.cutoff = cor.cutoff, cor.method = cor.method, min.sz = min.sz, perm.times = perm.times, pAdjustMethod = pAdjustMethod
    )
  }

  return(sig.corscore.df)
}
