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
  A <- as.matrix(assay(data, assay))
  print(k)

  if(is.null(k)){
    warning("Number of factors for NMF not specified. Using cross-validation to idenitfy optimal number of clusters")
      k <- find_num_factors(A)
  }


  model <- singlet::run_nmf(A, rank=k)

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
  if(is(target, "SpatialExperiment")){
    stopifnot(assay %in% assayNames(target))
  }

  loadings <- nmf_model$w
  #subset to the genes shared between the source and the target
  i<-intersect(rownames(target), rowData(source)$gene_name) # need to check for gene names in source object later
  rownames(loadings) <- rowData(source)$gene_name
  loadings<-loadings[rownames(loadings) %in% i,]
  loadings<-loadings[match(rownames(target),rownames(loadings)),]


  A <- assay(target, assay)
  print(head(loadings))
  options(RcppML.threads = 0) #line below doesn't work otherwise
  proj<-RcppML::project(data=as.matrix(A), w=loadings, threads=0, L1=0, mask=NULL)

  factors <- t(proj)/nmf_model$d #scale by the constant factor
  print(head(factors))
  colnames(factors) <- paste0("NMF", 1:ncol(factors))

  return(factors)
}
