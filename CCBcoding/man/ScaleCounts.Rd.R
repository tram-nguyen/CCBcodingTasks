#' Title
#'
#' @param mat
#'
#' @return Scaled matrix
#' @export
#'
#' @examples
#' ScaledMat <- scaleCounts(mat)
#'
scaleCounts <- function(mat) {
  col_sums <- colSums(mat)
  normalized_mat <- (mat / col_sums) * 10^6
  return(normalized_mat)
}
