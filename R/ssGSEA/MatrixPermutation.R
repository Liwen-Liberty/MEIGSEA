#' Construct random matrix
#'
#' @param  mat The original matrix
#' @param  is.row TRUE indicate perturb row label, otherwise perturb column label
#'
#' @return random matrix

MatrixPermutation <- function(mat, is.row) {
  if (is.row) {
    # Obtain the random row names
    permutation.rowTags <- sample(rownames(mat), size = nrow(mat), replace = FALSE)
    permutation.mat <- mat[permutation.rowTags, ]
    rownames(permutation.mat) <- rownames(mat)
  } else {
    # Obtain the random column names
    permutation.colTags <- sample(colnames(mat), size = ncol(mat), replace = FALSE)
    permutation.mat <- mat[, permutation.colTags]
    colnames(permutation.mat) <- colnames(mat)
  }

  return(permutation.mat)
}
