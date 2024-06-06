#' The traam class
#'
#' @param nmf_model
#' @param source
#' @param targets
#' @param multinom
#' @param cors
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
traam <- function(nmf_model, source, targets, multinom, cors, ...){
  res <- list(nmf=nmf_model, source=source, targets=targets, multinom_mod=multinom, correlations=cors)
  class(res) <- "traam"
  return(res)
}
