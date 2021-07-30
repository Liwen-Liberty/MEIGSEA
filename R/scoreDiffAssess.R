#' ssGSEA strategy
#'
#' @param  exp.profile Expression profile with row as genes, column as samples.
#' @param  sig.list Multiple signature gene sets to be enriched.
#' @param  min.sz At least how many immunophenotype genes are matched in the
#' expression profile, min.sz = 15 by default.
#' @param  score.norm Whether ssGSEA scores need to be standardized based on random
#' background, default is FALSE.
#' @param  perm.times Permutation times when the normalization of ssGSEA score
#' is needed, 100 by default.
#' @param  is.weightInteg Whether or not the ssGSEA scores need to be weighted
#' @param  sig.weight When is.weightInteg = TRUE, set the weight vector
#' for the different signatures
#' @param  sig.group When is.weightInteg = TRUE, set the group to which different signatures belong
#' @param  sample.group Expression profile sample column label, used to calculate
#' ssGSEA score differences between samples.
#' @param  case/control In default, case = "1", control = "0", or the specific character
#' representing mutant (case) and wild type (control) according to 'sample.group'.
#' @param  pAdjustMethod Method used to adjust p values, "BH" by default.
#'
#' @return ssGSEA score for each sample and differences in ssGSEA scores between mutant and wild type


scoreDiffAssess <- function(exp.profile, sig.list, min.sz = 15,
                            score.norm = FALSE, perm.times = 1000,
                            is.weightInteg = FALSE, sig.weight = NULL, sig.group = NULL,
                            sample.group, case, control, pAdjustMethod = "BH") {
  ###
  #' 1. Calculate the ssGSEA score for each sample
  #'
  ################################

  if (all(table(sig.group) != 2)) {
    score.norm <- FALSE
  }

  #' ######
  # (1)parameter setting --------------------
  #' ######
  param <- list(
    expr_mat = as.matrix(exp.profile),
    gset.idx.list = sig.list,
    min.sz = min.sz,
    method = "ssgsea",
    ssgsea.norm = FALSE,
    score.norm = FALSE,
    perm.times = perm.times
  )

  #' ######
  # (2) Calculate the ssGSEA score--------------------
  #' ######

  if (param$score.norm) {
    # Standardized ssGSEA scores based on random background
    ssGSEAScore <- ssGSEAnormalization(param$expr_mat, param$"gset.idx.list", param$perm.times)
    rownames(ssGSEAScore) <- colnames(exp.profile)
  } else {
    # Calculate ssGSEA score
    ssGSEAScore <- GSVA::gsva(expr = param$expr_mat, gset.idx.list = param$"gset.idx.list", method = param$method, ssgsea.norm = param$"ssgsea.norm")
    # Scale ssGSEA score
    ssGSEAScore <- apply(ssGSEAScore, 1, scale)
    rownames(ssGSEAScore) <- colnames(exp.profile)
  }


  ###
  #' 2. If necessary, combine the different signature scores weighted
  #'
  ################################
  if (is.weightInteg) {
    ssGSEAScore <- t(ssGSEAScore) * sig.weight
    ssGSEAScore <- rowsum(ssGSEAScore, sig.group)
    ssGSEAScore <- t(ssGSEAScore)
  }


  ###
  #' 3. Compare the ssGSEA score between mutant and WT samples using Wilcoxon rank sum test
  #'
  ################################

  data.merge <- cbind.data.frame(ssGSEAScore, sample_group = sample.group)
  wilcox.res <- lapply(colnames(ssGSEAScore), function(s) {
    tmp.data <- data.merge[, c("sample_group", s)]
    tmp.data$sample_group <- factor(tmp.data$sample_group, levels = c(case, control))
    tmp.result <- coin::wilcox_test(as.formula(paste0(s, " ~ sample_group")), data = tmp.data)
    refine.res <- data.frame(signature = s, Zvalue = coin::statistic(tmp.result, type = "standardized"), p.value = coin::pvalue(tmp.result))
    return(refine.res)
  })
  wilcox.res <- do.call(rbind.data.frame, wilcox.res)
  colnames(wilcox.res)[2] <- "Zvalue"


  # return results
  comb.res <- list(wilcox.res = wilcox.res, ssGSEAScore = ssGSEAScore)
  return(comb.res)
}
