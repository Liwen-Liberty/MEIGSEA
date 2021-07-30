#' Use contingency table test to identify  co-occurrence or mutual exclusion between mutations
#'
#' @param mutationMat Mutation matrix
#' @param outfile File path to the output upper triangle heatmap file
#' @param min_count Set the minimum number of mutations, and only test the genes whose mutations are larger than min_count
#'
#' @return Test results, including odds ratio and p value

get_mut_inter_fisher <- function(mutationMat = NULL, outfile = NULL, min_count = 2) {
  # Extracte genes with mutation number greater than min_count
  ex_idx <- which(rowSums(mutationMat) > min_count)
  if (length(ex_idx) > 0) {
    mutationMat_sub <- mutationMat[ex_idx, ]

    test_genes <- rownames(mutationMat_sub)
    test_gene_combn <- combn(test_genes, 2)
    res_list <- list()
    for (comb_set_idx in 1:dim(test_gene_combn)[2]) {
      # Build a fourfold table
      test_table <- matrix(NA, nrow = 2, ncol = 2)
      rownames(test_table) <- paste0(test_gene_combn[1, comb_set_idx], c("_MUT", "_WT"))
      colnames(test_table) <- paste0(test_gene_combn[2, comb_set_idx], c("_MUT", "_WT"))

      test_table_tmp <- table(rowSums(t(mutationMat_sub[test_gene_combn[, comb_set_idx], ])))
      test_table[1, 1] <- test_table_tmp["2"]
      test_table[2, 2] <- test_table_tmp["0"]

      test_table_tmp2 <- table(mutationMat_sub[test_gene_combn[1, comb_set_idx], ] - mutationMat_sub[test_gene_combn[2, comb_set_idx], ])
      test_table[1, 2] <- test_table_tmp2["1"]
      test_table[2, 1] <- test_table_tmp2["-1"]
      # convert NA to 0
      na_idx <- which(is.na(test_table))
      if (length(na_idx) > 0) test_table[na_idx] <- 0

      # statistic test, get odd ratio and p value
      res_tmp <- getChisqTest(testTable = test_table)
      res_list[[comb_set_idx]] <- as.data.frame(res_tmp[c("oddsratio", "pvalue")])
    }
    res <- Reduce("rbind", res_list)
    rownames(res) <- apply(test_gene_combn, 2, paste, collapse = "_")

    ## Whether to output a heat map
    if (!is.null(outfile)) {

      # Extract the gene names
      gene_names <- unique(unlist(strsplit(rownames(res), split = "_")))
      # Build the heatmap matrix and the matching P value matrix
      heatmap_arr <- matrix(NA, nrow = length(gene_names), ncol = length(gene_names))
      pvalue_arr <- matrix(NA, nrow = length(gene_names), ncol = length(gene_names))
      rownames(heatmap_arr) <- colnames(heatmap_arr) <- gene_names
      rownames(pvalue_arr) <- colnames(pvalue_arr) <- gene_names
      # Extract the index
      heatmap_arr_idx <- Reduce("rbind", strsplit(rownames(res), split = "_"))
      for (i in 1:dim(heatmap_arr_idx)[1]) {
        heatmap_arr[heatmap_arr_idx[i, 1], heatmap_arr_idx[i, 2]] <- res[i, "oddsratio"]
        pvalue_arr[heatmap_arr_idx[i, 1], heatmap_arr_idx[i, 2]] <- res[i, "pvalue"]
        heatmap_arr[heatmap_arr_idx[i, 1], heatmap_arr_idx[i, 1]] <- 1
        pvalue_arr[heatmap_arr_idx[i, 1], heatmap_arr_idx[i, 1]] <- 1
      }
      #  log2 conversion
      log2heatmap_arr <- log2(heatmap_arr)
      # replace infinity value
      if (length(which(log2heatmap_arr == Inf)) > 0) {
        log2heatmap_arr[which(log2heatmap_arr == Inf)] <- 1000
      }
      if (length(which(log2heatmap_arr == (-Inf))) > 0) {
        log2heatmap_arr[which(log2heatmap_arr == (-Inf))] <- -1000
      }
      # color setting
      heat_colors <- colorRampPalette(c("blue", "white", "white", "red"))(30)
      color_breaks <- seq(-4, 4, length = 31)
      # output
      pdf(outfile)
      print(pheatmap::pheatmap(mat = t(log2heatmap_arr), color = heat_colors, breaks = color_breaks, cluster_rows = FALSE, cluster_cols = FALSE))
      dev.off()
    }
    return(res)
  } else {
    return(NA)
  }
}
