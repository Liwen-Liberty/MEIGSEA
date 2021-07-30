#'  Identifying somatic gene mutations associated with immunophenotype gene sets.
#'
#' @param  raw.count.profile RNA-seq data in raw count or microarray expression profile
#' with row as genes, column as samples.
#' @param  is.rawcount whether or not the 'raw.count.profile' is RNA-seq raw count.
#' For microarray expression data as input, is.rawcount = FALSE.
#' @param  gene.length Gene length annotation file, the first column is the gene name, the second
#' column is the gene length, if 'gene.length' is not provided, the 'tpm.exp' needs to be given.
#' @param  tpm.exp Given tpm expression profile (optional). If the 'gene.lengths' given,
#' tpm.exp is not needed.
#' @param  mutationCalling.map A matrix of gene mutation status (0 or 1, in character),
#' row as genes, column as samples.
#' @param  signature.df  A data frame containing immunophenotype signature: the first column
#' (named 'Gene') is the gene ID, the second column (named 'setAnno') is the immunophenotype name.
#' @param  case/control In default, case = "1", control = "0", or the specific character
#' representing mutant (case) and wild type (control) according to 'mutationCalling.map'.
#' @param  cor.cutoff Setting the cutoff of correlation coefficient when refining the
#' immunophenotype, cor.cutoff = NULL by default, means all genes with statistical significance
#' will be retained
#' @param  cor.method method used to compute the correlation, either "spearman", "pearson"
#' or "kendall", "spearman" by default
#' @param  min.sz At least how many immunophenotype genes are matched in the expression profile,
#' min.sz = 15 by default.
#' @param  perm.times Perumutation times when the normalization of ssGSEA score is needed,
#' 100 by default.
#' @param  pAdjustMethod Method used to adjust p values, "BH" by default.
#'
#' @return  A list including detailed corrlation results, refined immunophenotype signature
#' and all correlated mutation-signature pairs



IntegratedGSEA <- function(raw.count.profile, is.rawcount = TRUE, gene.length = NULL, tpm.exp = NULL, mutationCalling.map, signature.df,
                           case = "1", control = "0",
                           cor.cutoff = NULL, cor.method = "spearman",
                           min.sz = 15, perm.times = 100, pAdjustMethod = "BH") {



  ###
  # 1.Whether the samples of expression data and mutation data match --------------------
  #'
  #' ###############################
  if (!identical(colnames(raw.count.profile), colnames(mutationCalling.map))) {
    cat("Sample names do not match between expression profile and mutation profile!")
    return(NULL)
  }


  ###
  # 2.raw count to TPM if needed--------------------
  #'
  #' ###############################
  if (is.null(tpm.exp) & is.rawcount) {
    TPM.profile <- CountToTPM(raw.count.profile, gene.length = NULL)
    cat("** TPM transformed finished ! \n")
  } else if (is.null(tpm.exp) & (!is.rawcount)) {
    TPM.profile <- raw.count.profile
  } else {
    TPM.profile <- tpm.exp
  }



  ###
  # 3.refining immunophenotype signature --------------------
  #'
  ################################
  genesetMeancor.res <- genesetMeanFilter(exp.profile = TPM.profile, cor.method = cor.method, signature.df = signature.df, cor.cutoff = cor.cutoff)
  # If too few genes remained(dissatisfy the calculation) due to sample size and other reasons, all genes before refining will be used
  if (all(table(genesetMeancor.res$df$setAnno) < min.sz)) {
    genesetMeancor.res$df <- data.frame(gene = signature.df$Gene, setAnno = signature.df$setAnno)
    genesetMeancor.res$df$label <- paste0(genesetMeancor.res$df$setAnno, "_positive")
    genesetMeancor.res$list <- split(genesetMeancor.res$df$gene, genesetMeancor.res$df$setAnno)
  } else {
    # Check signature with both positive and negative subset < min.sz
    check.size <- tapply(genesetMeancor.res$df$label, genesetMeancor.res$df$setAnno, table)
    check.size <- sapply(check.size, function(x) {
      return(ifelse(all(x < min.sz), T, F))
    })
    if (length(which(check.size)) > 0) {
      restore.setAnno <- names(check.size)[which(check.size)]

      genesetMeancor.res$df <- genesetMeancor.res$df[-which(genesetMeancor.res$df$setAnno %in% restore.setAnno), ]
      restore.signature.df <- subset(signature.df, setAnno %in% restore.setAnno)
      refine.restore.signature.df <- data.frame(
        gene = restore.signature.df$Gene, setAnno = restore.signature.df$setAnno,
        label = paste0(restore.signature.df$setAnno, "_positive"),
        direction = "positive", signif = "sig", cor = 0, p.ad = 0
      )
      genesetMeancor.res$df <- rbind.data.frame(genesetMeancor.res$df, refine.restore.signature.df[, colnames(genesetMeancor.res$df)])
      genesetMeancor.res$list <- split(genesetMeancor.res$df$gene, genesetMeancor.res$df$setAnno)
    }
  }
  rownames(genesetMeancor.res$df) <- NULL
  cat("** gene set refined ! \n")



  ###
  # 4.identifing the gene mutation-signature correlations --------------------
  #' get the detailed results from each sub-method
  ################################
  GSEAbased.cor.list <- lapply(rownames(mutationCalling.map), function(g) {
    sample.group <- as.character(mutationCalling.map[g, ])

    #' ######
    # (1) one/two-tail GSEA strategy  --------------------
    #' ######
    one_twoGSEA.res <- tailedGSEA(
      exp.profile = raw.count.profile, is.rawcount = is.rawcount, sample.group = sample.group, case = case, control = control,
      sig.df = genesetMeancor.res$df, min.sz = min.sz, pAdjustMethod = pAdjustMethod,
      is.fgsea = FALSE, is.fisherPcombine = FALSE
    )

    #' ######
    # (2) ssGSEA strategy --------------------
    #' ######
    check.size <- sapply(genesetMeancor.res$list, length)
    if (any(check.size < min.sz)) {
      genesetMeancor.res$list <- genesetMeancor.res$list[-which(check.size < min.sz)]
    }
    print(paste0("*** length of signature list is ", length(genesetMeancor.res$list)))

    if (length(genesetMeancor.res$list) > 0) {
      ssGSEA.wilcox.res <- scoreDiffAssess(
        exp.profile = TPM.profile, sig.list = genesetMeancor.res$list, min.sz = min.sz,
        score.norm = FALSE, perm.times = perm.times,
        is.weightInteg = FALSE, sig.weight = NULL, sig.group = NULL,
        sample.group = sample.group, case = case, control = control, pAdjustMethod = pAdjustMethod
      )
    } else {
      ssGSEA.wilcox.res <- NULL
      cat("*** ssGSEA filed due to null signature set ! \n")
    }
    cat("*** ssGSEA wilcox finished ! \n")


    #' ######
    # (3) Format of results --------------------
    #' ######

    ### ssGSEA results
    wilcoxRes.df <- data.frame(ssGSEA.wilcox.res$wilcox.res, method = "ssGSEA")
    colnames(wilcoxRes.df) <- gsub("Zvalue", "statistic", colnames(wilcoxRes.df))

    ### oneTail.GSEA results
    oneTailRes <- as.data.frame(one_twoGSEA.res$oneTailRes)
    if (nrow(oneTailRes) > 0) {
      GSEAres.df <- data.frame(oneTailRes[, c("ID", "NES", "pvalue")], method = "oneTail.GSEA")
    } else {
      GSEAres.df <- data.frame(ID = NA, NES = NA, pvalue = NA, method = "oneTail.GSEA") # oneTailGSEA had no results
    }

    ### twoTail.GSEA results
    twoTailRes <- one_twoGSEA.res$twoTailRes$twoTailRes
    if (!is.null(twoTailRes)) {
      GSEAres.df <- rbind.data.frame(GSEAres.df, data.frame(twoTailRes[, c("ID", "NES", "pvalue")], method = "twoTail.GSEA"))
    } else {
      GSEAres.df <- rbind(GSEAres.df, data.frame(ID = NA, NES = NA, pvalue = NA, method = "twoTail.GSEA")) # twoTailGSEA had no results
    }
    colnames(GSEAres.df)[c(1, 2, 3)] <- c("signature", "statistic", "p.value")


    comb.relations.df <- rbind.data.frame(wilcoxRes.df, GSEAres.df)
    all.res <- list(one_twoGSEA.res = one_twoGSEA.res, ssGSEA.wilcox.res = ssGSEA.wilcox.res, comb.relations.df = comb.relations.df)

    return(all.res)
  })
  names(GSEAbased.cor.list) <- rownames(mutationCalling.map)


  ###
  # 5.p values integration and correlation score calculation --------------------
  #'
  ################################

  corscore.list <- lapply(GSEAbased.cor.list, function(corRes.oneGene) {
    corRes.oneSig.list <- lapply(unique(corRes.oneGene$comb.relations.df$signature), function(one.sig) {
      corRes.oneSig <- subset(corRes.oneGene$comb.relations.df, signature == one.sig)
      ### confidence coefficient (CC) of correlation
      corRes.oneSig$CC <- ifelse(corRes.oneSig$p.value < 0.05, 1, 0)
      ### P values integration by Fisher’s method
      fisher.p <- fisher.method(matrix(corRes.oneSig$p.value, nrow = 1),
        p.corr = pAdjustMethod,
        zero.sub = 1e-05,
        na.rm = TRUE
      )
      #######
      # Integrate results of ssGSEA and two-tail GSEA to determine the direction of association
      #######
      tmp.2res <- corRes.oneSig[which(corRes.oneSig$method %in% c("ssGSEA", "twoTail.GSEA")), c("statistic", "p.value")]
      if (all(tmp.2res$p.value > 0.05)) {
        direction <- NA
      } else if (any(tmp.2res$p.value > 0.05)) {
        direction <- ifelse(tmp.2res$statistic[which(tmp.2res$p.value < 0.05)] > 0, "positive", "negative")
      } else {
        if (purrr::reduce(tmp.2res$statistic, `*`) > 0) {
          direction <- ifelse(tmp.2res$statistic[1] > 0, "positive", "negative")
        } else {
          direction <- "unclear"
        }
      }

      tmp.res <- data.frame(
        signature = corRes.oneSig$signature[1],
        twoTailedGSEA.NES = ifelse(is.element("twoTail.GSEA", corRes.oneSig$method), corRes.oneSig$statistic[which(corRes.oneSig$method == "twoTail.GSEA")], NA),
        ssGSEA.Z.value = ifelse(is.element("ssGSEA", corRes.oneSig$method), corRes.oneSig$statistic[which(corRes.oneSig$method == "ssGSEA")], NA),
        fisher.p,
        direction = direction,
        weight = sum(corRes.oneSig$CC) / 3
      )

      tmp.res$weightedcorscore <- tmp.res$weight * tmp.res$S
      return(tmp.res)
    })
    corRes.oneSig.df <- do.call(rbind.data.frame, corRes.oneSig.list)
    corRes.oneSig.df$p.adj <- p.adjust(corRes.oneSig.df$p.value, method = pAdjustMethod)
    return(corRes.oneSig.df)
  })
  names(corscore.list) <- names(GSEAbased.cor.list)

  corscore.df <- frameInlist2addOneFrame(corscore.list, "driver")
  rownames(corscore.df) <- NULL

  return(list(GSEAbased.cor.list = GSEAbased.cor.list, corscore.df = corscore.df, genesetMeancor.res = genesetMeancor.res))
}

# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
fisher.method <- function(pvals, method = c("fisher"), p.corr = c(
                            "bonferroni", "BH",
                            "none"
                          ), zero.sub = 0.00001, na.rm = FALSE, mc.cores = NULL) {
  stopifnot(method %in% c("fisher"))
  stopifnot(p.corr %in% c("none", "bonferroni", "BH"))
  stopifnot(all(pvals >= 0, na.rm = TRUE) & all(pvals <= 1, na.rm = TRUE))
  stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 1)
  if (is.null(dim(pvals))) {
    stop("pvals must have a dim attribute")
  }
  p.corr <- ifelse(length(p.corr) != 1, "BH", p.corr)
  ## substitute p-values of 0
  pvals[pvals == 0] <- zero.sub
  if (is.null(mc.cores)) {
    fisher.sums <- data.frame(do.call(rbind, apply(pvals, 1, fisher.sum,
      zero.sub = zero.sub, na.rm = na.rm
    )))
  } else {
    fisher.sums <- parallel::mclapply(1:nrow(pvals), function(i) {
      fisher.sum(pvals[i, ], zero.sub = zero.sub, na.rm = na.rm)
    }, mc.cores = mc.cores)
    fisher.sums <- data.frame(do.call(rbind, fisher.sums))
  }

  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1 - pchisq(fisher.sums$S, df = 2 * fisher.sums$num.p)
  fisher.sums$p.adj <- switch(p.corr,
    bonferroni = p.adjust(fisher.sums$p.value, "bonferroni"),
    BH = p.adjust(fisher.sums$p.value, "BH"),
    none = fisher.sums$p.value
  )
  return(fisher.sums)
}
# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
fisher.sum <- function(p, zero.sub = 0.00001, na.rm = FALSE) {
  if (any(p > 1, na.rm = TRUE) || any(p < 0, na.rm = TRUE)) {
    stop("You provided bad p-values")
  }
  stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 1)
  p[p == 0] <- zero.sub
  if (na.rm) {
    p <- p[!is.na(p)]
  }
  S <- -2 * sum(log(p))
  res <- data.frame(S = S, num.p = length(p))
  return(res)
}
