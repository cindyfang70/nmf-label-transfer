
#' Compute correlation of each NMF factor with annotations in the source dataset
#'
#' @param factors A Data.Frame of NMF factors
#' @param cor_variable A vector of per-cell annotations or discrete technical levels for your source dataset. For example, cell-type labels or sample ID.
#'
#' @import stats
#'
#' @return A NMF factor by unique annotations matrix of correlations between each factor and each unique annotation.
#' @export
compute_factor_correlations <- function(factors, cor_variable) {

  cor_variable <- as.factor(cor_variable)
  print('computing correlation')
  factors.ind <- model.matrix(~0+cor_variable, na.action=na.exclude)
  print(dim(factors.ind))
  colnames(factors.ind) <- unlist(lapply(strsplit(colnames(factors.ind), split = "cor_variable"), "[", 2))
  factors.ind <- as.data.frame(factors.ind)
  print(dim(factors.ind))
  print(dim(factors))
  cor.mat <- cbind(factors, factors.ind)
  print(head(cor.mat))

  M <- stats::cor(cor.mat, use="complete.obs")
  M <- M[grepl("NMF", rownames(M)), !grepl("NMF", colnames(M))]

  return(M)
}


identify_factors_representing_annotations <- function(correlations, technical_correlations, n_factors_per_annot=3, cor_cutoff=0.3){
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

  print(selected_factors)

  selected_factors_keep <- c()
  for (i in 1:length(selected_factors)){
    abs_fct_sample_cor <- abs(technical_correlations[selected_factors[[i]],])
    if (!any(abs_fct_sample_cor > cor_cutoff)){
      print(selected_factors[[i]])
      selected_factors_keep <- c(selected_factors_keep, selected_factors[[i]])
    }
  }
  print(selected_factors_keep)
  return(selected_factors_keep)
}
