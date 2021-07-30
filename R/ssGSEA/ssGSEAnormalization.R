#' Normalization of ssGSEA score
#'
#' @param  expressionProfile Expression profile with row as genes, column as samples.
#' @param  signatureList Multiple signature gene sets to be enriched.
#' @param  perm.times Permutation times when the normalization of ssGSEA score is needed, 100 by default.
#'
#' @return A matrix of the normalized ssGASEA score

ssGSEAnormalization <- function(expressionProfile, signatureList, perm.times = 1000) {
  ###
  #' 1.Construst random expression profile 'permutation.exp'
  #'
  ################################

  sample.number <- ncol(expressionProfile)

  if (sample.number < 50) {
    # If the sample size < 50, randomize the whole expression profile
    permutation.exp <- matrix(rep(as.matrix(expressionProfile), times = (perm.times + 1)), nrow = nrow(expressionProfile))
    for (i in 1:perm.times) {
      tmp <- MatrixPermutation(as.matrix(expressionProfile), is.row = TRUE)

      permutation.exp[, (i * sample.number + (1:sample.number))] <- as.matrix(tmp)
    }
    rownames(permutation.exp) <- rownames(expressionProfile)
  }


  ###
  #' 2. Calculate ssGSEA score
  #'
  ################################
  score.list <- lapply(1:sample.number, function(i) {
    # random expression profile for one sample
    if (sample.number < 50) {
      exp.forOneSample <- permutation.exp[, ((0:perm.times) * sample.number + i)]
    } else {
      random.idx <- replicate(perm.times, sample(1:nrow(expressionProfile)))
      exp.forOneSample <- cbind(expressionProfile[, i], matrix(expressionProfile[random.idx, i], nrow = nrow(expressionProfile)))
    }

    # Calculate ssGSEA score
    score.forOneSample <- CustomizedssGSEA(expressionProfile = exp.forOneSample, signatureList = signatureList)

    return(score.forOneSample)
  })


  ###
  #' 3. Normalization of ssGSEA score
  #'
  ################################
  zscore.matrix <- sapply(score.list, function(score.forOneSample) {
    check.score <- score.forOneSample * score.forOneSample[, 1]
    score.forOneSample[which(check.score < 0)] <- NA

    # z-score normalization
    zScore.forOneSample <- apply(score.forOneSample, 1, function(x) {
      if (length(which(!is.na(x))) < 2) {
        zScore <- x[1]
      } else {
        zScore <- scale(x)[1]
      }
      return(zScore)
    })
    return(zScore.forOneSample)
  })
  return(t(zscore.matrix))
}
