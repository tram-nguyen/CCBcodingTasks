#' Convert count data to data frame
#'
#' @param df
#'
#' @return data frame in long format
#' @export
#'
#' @examples
#' LongCountDF <- edit_df(df)
#'
edit_df <- function(df) {

  library(dplyr)
  df %>%
    rownames_to_column(var = "GeneID") %>%
    melt(id.vars = c("Condition", "GeneID"), variable.name = "Sample", value.name = "ReadCount")

}
