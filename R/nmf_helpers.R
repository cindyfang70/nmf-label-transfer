#' Run NMF on a dataset
#'
#' @param data SingleCellExperiment or SpatialExperiment object
#' @param assay string indicating the assay to run NMF on
#' @param k integer indicating the number of factors for NMF
#' @param seed a random seed
#'
#' @return NMF model object
#'
#' @import SingleCellExperiment
#' @import RcppML
#' @export
run_nmf<- function(data, assay, k, seed=1237){
  A <- assays(data, assay)

  if(missing(k)){
    warning("Number of factors for NMF not specified. Using cross-validation to idenitfy optimal number of clusters")
      k <- find_num_factors(A)
  }


  model <- RcppML::nmf(A, k = k, seed=seed)

  return(model)
}

find_num_factors <- function(A, ranks = c(50,100,200)){
  if(any(ranks >= ncol(A))){
    stop("ranks must be less than the number of columns in A")
    }
  cv <- singlet::cross_validate_nmf(A, ranks = ranks,
                                    n_replicates = 3,
                                    verbose=3)
  num_factors <- singlet::GetBestRank(cv)

  return(num_factors)
}

project_factors <- function(source, target, assay, nmf_model){
  loadings <- nmf_model$w
  #subset to the genes shared between the source and the target
  i<-intersect(rownames(target),rownames(rowData(source)))
  loadings<-loadings[rownames(loadings) %in% i,]
  loadings<-loadings[match(rownames(target),rownames(loadings)),]

  A <- assays(target, assay)
  proj<-project(A=A, w=loadings, L1=0)

  factors <- t(proj)/nmf_model$d #scale by the constant factor
  colnames(factors) <- paste0("NMF", ncol(factors))

  return(factors)
}
