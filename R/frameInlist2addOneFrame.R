#' Convert list to data.frame, the new data frame adds one column to
#' each element's data frame
#'
#' @param  specific_list A special list in which each element is a data
#' frame with the same structure
#' @param  new_colname Specifies the name of the column to be added
#'
#' @return a data frame

frameInlist2addOneFrame <- function(specific_list, new_colname) {
  combine_df <- do.call(rbind.data.frame, specific_list)
  combine_df[, new_colname] <- rep(names(specific_list), times = sapply(specific_list, nrow))

  return(combine_df)
}
