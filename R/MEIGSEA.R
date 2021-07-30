#'  Main Function for identifying somatic gene mutations associated with
#' immunophenotype gene sets.
#'
#' @param  raw.count.profile RNA-seq data in raw count or microarray expression
#' profile with row as genes, column as samples.
#' @param  is.rawcount whether or not the 'raw.count.profile' is RNA-seq raw count.
#' For microarray expression data as input, is.rawcount = FALSE.
#' @param  gene.length Gene length annotation file, the first column is the gene name,the second
#' column is the gene length, if 'gene.length' is not provided, the 'tpm.exp' needs to be given.
#' @param  tpm.exp Given tpm expression profile (optional). If the 'gene.lengths' given,
#' tpm.exp is not needed.
#' @param  mutationCalling.map A matrix of gene mutation status (0 or 1, in character),
#' row as genes, column as samples.
#' @param  signature.df signature A data frame containing immunophenotype signature:
#' the first column (named 'Gene') is the gene ID, the second column (named 'setAnno')
#' is the immunophenotype name.
#' @param  case/control In default, case = "1", control = "0", or the specific character
#' representing mutant (case) and wild type (control) according to 'mutationCalling.map'.
#' @param  cor.cutoff Setting the cutoff of correlation coefficient when refining the
#' immunophenotype, cor.cutoff = NULL by default, means all genes with statistical
#' significance will be retained
#' @param  cor.method method used to compute the correlation, either "spearman", "pearson"
#' or "kendall", "spearman" by default
#' @param  min.sz At least how many immunophenotype genes are matched in the expression
#' profile, min.sz = 15 by default.
#' @param  perm.times Perumutation times when the normalization of ssGSEA score is needed,
#' 100 by default.
#' @param  pAdjustMethod Method used to adjust p values, "BH" by default.
#' @param cond.corr.anal Whether to perform conditional association analysis,
#' cond.corr.anal=TRUE by default.
#' @param sig.cutoff cutoff for significant associations, sig.cutoff = 0.01 by default.
#' @param  output Whether to store the results locally, output = TRUE by default.
#' @param  output.dir Character representing the storage path of results.
#' @param  interested.mut vector of interested gene mutations, for which the detailed
#' correlation results will be given as independent files
#'
#' @return  a list including detailed corrlation results, refined immunophenotype signature
#' and all correlated mutation-signature pairs

MEIGSEA <- function(raw.count.profile, is.rawcount = TRUE, gene.length = NULL, tpm.exp = NULL, mutationCalling.map, signature.df,
                    case = "1", control = "0",
                    cor.cutoff = NULL, cor.method = "spearman",
                    min.sz = 15, perm.times = 100, pAdjustMethod = "BH",
                    cond.corr.anal = TRUE, sig.cutoff = 0.01,
                    output = TRUE, output.dir = getwd(), interested.mut = NULL) {

  #  Identify somatic gene mutations associated with immunophenotype gene sets
  original.res <- IntegratedGSEA(
    raw.count.profile = raw.count.profile, is.rawcount = is.rawcount, gene.length = gene.length, tpm.exp = tpm.exp, mutationCalling.map = mutationCalling.map,
    signature.df = signature.df, case = case, control = control, cor.cutoff = cor.cutoff, cor.method = cor.method, min.sz = min.sz, perm.times = perm.times, pAdjustMethod = pAdjustMethod
  )

  corscore.df <- original.res$corscore.df
  rownames(corscore.df) <- NULL

  if (cond.corr.anal) {
    # Conditional Correlation Analysis
    sig.corscore.df <- CondCorrAnal(corscore.df,
      raw.count.profile = raw.count.profile, is.rawcount = is.rawcount, gene.length = gene.length, tpm.exp = tpm.exp, mutationCalling.map = mutationCalling.map,
      sig.cutoff = sig.cutoff, signature.df = signature.df, case = case, control = control, cor.cutoff = cor.cutoff, cor.method = cor.method, min.sz = min.sz, perm.times = perm.times, pAdjustMethod = pAdjustMethod
    )
    corscore.df$FILTER <- "Non-significant"
    corscore.df$CoMutGenes <- NA
    corscore.df <- t(apply(corscore.df, 1, function(x) {
      sig.index <- which(sig.corscore.df$signature == x[["signature"]] & sig.corscore.df$driver == x[["driver"]])
      if (length(sig.index) != 0) {
        x["FILTER"] <- sig.corscore.df[sig.index, "FILTER"]
        x["CoMutGenes"] <- sig.corscore.df[sig.index, "CoMutGenes"]
      }
      x
    }))
    corscore.df
  }

  # output files
  if (output) {
    # Correlated mutation-signature pairs
    write.table(corscore.df, sprintf("%s/%s", output.dir, "correlation_pairs.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    # Data frame of refined immunophenotype signature
    write.table(original.res$genesetMeancor.res$df, sprintf("%s/%s", output.dir, "refined_immunophenotype.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    # Significant correlation pairs
    if (cond.corr.anal) write.table(corscore.df[which(corscore.df[, "FILTER"] == "PASS"), 1:11], sprintf("%s/%s", output.dir, "sig_correlation_pairs.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  # Detailed corrlation results of interested gene mutations
  if (!is.null(interested.mut)) {
    sapply(interested.mut, function(g) {
      write.table(original.res$GSEAbased.cor.list[[g]]$comb.relations.df, sprintf("%s/%s", output.dir, paste0("detailed subMethods results_", g, ".txt")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })
  }

  # return results
  if (cond.corr.anal) original.res$corscore.df <- corscore.df
  return(original.res)
}
