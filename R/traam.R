#' The traam class
#'
#' @param nmf_model An NMF model from `singlet`
#' @param source A SpatialExperiment or SingleCellExperiment object
#' @param targets A list of SpatialExperiment or SingleCellExperiment object(s)
#' @param multinom An object of class `multinom` from `nnet`
#' @param cors A matrix of correlations between NMF factors and annotations
#'
#' @return An object of class traam
#' @export
traam <- function(nmf_model, source, targets, multinom, cors){
  res <- list(nmf=nmf_model, source=source, targets=targets, multinom_mod=multinom, correlations=cors)
  class(res) <- "traam"
  return(res)
}
