
#' Compute correlation of each NMF factor with annotations in the source dataset
#'
#' @param factors A Data.Frame of NMF factors
#' @param source_annotations A vector of per-cell annotations for your source dataset. For example, cell-type labels or sample ID.
#'
#' @import stats
#'
#' @return A NMF factor by unique annotations matrix of correlations between each factor and each unique annotation.
#' @export
compute_factor_correlations <- function(factors, source_annotations) {

  labels <- unique(source_annotations)
  n_labels <-length(labels)
  n_obs <- length(source_annotations)

  factors.ind <- matrix(,nrow=n_obs, ncol=n_labels)
  for (i in 1:n_labels){
    label <- labels[[i]]
    factors.ind[,i] <- as.integer(source_annotations == label)
  }

  colnames(factors.ind) <- labels
  factors.ind <- as.data.frame(factors.ind)
  cor.mat <- cbind(factors, factors.ind)

  M <- stats::cor(cor.mat, use="complete.obs")
  M <- M[grepl("NMF", rownames(M)), !grepl("NMF", colnames(M))]

  return(M)
}

identify_factors_representing_annotations <- function(correlations, n_factors_per_annot=3){
  n_domains <- ncol(correlations)
  selected_factors <- c()
  for (i in 1:n_domains){
    domain <- correlations[,i]
    domain <- domain[! names(domain) %in% selected_factors]
    ord_domain <- domain[order(abs(domain), decreasing=TRUE)]
    top_fcts <- names(ord_domain[1:n_factors_per_annot])

    new.fcts <- setdiff(top_fcts, selected_factors) # keep the factor only if it's not already been included
    selected_factors <- c(selected_factors, new.fcts)
  }
  selected_factors <- selected_factors[!is.na(selected_factors)]
  return(selected_factors)
}
