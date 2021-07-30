#' Contingency table test
#'
#' @param testTable Fourfold table
#'
#' @return Test results, including odds ratio,  p value and test method

getChisqTest <- function(testTable = NULL) {
  testTable <- as.matrix(testTable)
  # Calculate theoretical frequency
  thFreq <- function(o1, o2, total) {
    (max(o1, o2) / total) * min(o1, o2)
  }
  n1 <- thFreq(sum(testTable[1, ]), sum(testTable[, 1]), sum(testTable))
  n2 <- thFreq(sum(testTable[1, ]), sum(testTable[, 2]), sum(testTable))
  n3 <- thFreq(sum(testTable[2, ]), sum(testTable[, 1]), sum(testTable))
  n4 <- thFreq(sum(testTable[2, ]), sum(testTable[, 2]), sum(testTable))
  theorFreq <- c(n1, n2, n3, n4)
  testRes <- NULL
  # According to the size of data, different types of contingency table test
  if (sum(testTable) > 40) {
    if (!any(theorFreq < 5)) {
      testRes <- chisq.test(testTable, correct = FALSE)
    }
    if (sum(theorFreq < 5) == 1) {
      testRes <- chisq.test(testTable, correct = TRUE)
    }
    if (sum(theorFreq < 5) > 1 || any(theorFreq < 1)) {
      testRes <- fisher.test(testTable)
    }
  } else {
    testRes <- fisher.test(testTable)
  }

  # Calculate odds ratio
  oddsratioRes <- fmsb::oddsratio(testTable)

  resList <- list(oddsratioRes$estimate, testRes$p.value, oddsratioRes$conf.int, oddsratioRes$method, testRes$method)
  names(resList) <- c("oddsratio", "pvalue", "oddsratio_95CI", "oddsratio_method", "test_method")
  return(resList)
}
