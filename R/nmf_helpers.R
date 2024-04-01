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
run_nmf<- function(data, assay, k=NULL, seed=1237){

  #add in checks for k
  message("Running NMF")

  if(is.null(k)){
    warning("Number of factors for NMF not specified. Using cross-validation to idenitfy optimal number of factors.", immediate. = TRUE)
      #k <- find_num_factors(A)
    model <- run_rank_determination_nmf(data, assay)
  }else{
    A <- as.matrix(assay(data, assay))
    model <- singlet::run_nmf(A, rank=k)
  }

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

#' run_rank_determination_nmf
#'
#' @param data A SingleCellExperiment or SpatialExperiment object
#' @param assay
#'
#' @return
#' @import RcppML singlet SummarizedExperiment SingleCellExperiment S4Vectors, SpatialExperiment
#' @examples
run_rank_determination_nmf <- function(data, assay){
  requireNamespace("SingleCellExperiment")
  requireNamespace("singlet")
  requireNamespace("SummarizedExperiment")
  requireNamespace("S4Vectors")
  data_nmf <- singlet::RunNMF(data)
  nmf_mod <- S4Vectors::metadata(data_nmf)$nmf_model
  return(nmf_mod)
}

project_factors <- function(source, target, assay, nmf_model){
  if(is(target, "SpatialExperiment")){
    if (!(assay %in% assayNames(source))){
      stop(sprintf("Assay %s not found in the source dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }

    if(!("gene_name" %in% colnames(rowData(target)))){
      stop("Please provide gene symbols in your target dataset as a column named 'gene_name' in rowData.")
    }
  }

  if(is(source, "SpatialExperiment")){
    if (!(assay %in% assayNames(source))){
      stop(sprintf("Assay %s not found in the source dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }

    if(!("gene_name" %in% colnames(rowData(source)))){
      stop("Please provide gene symbols in your source dataset as a column named 'gene_name' in rowData.")
    }
  }



  loadings <- nmf_model$w
  #subset to the genes shared between the source and the target
  i<-intersect(rowData(target)$gene_name, rowData(source)$gene_name) # need to check for gene names in source object later

  if(length(i) == 0){
    stop("No intersecting genes between target and source dataset.")
  }

  rownames(loadings) <- rowData(source)$gene_name
  loadings<-loadings[rownames(loadings) %in% i,]
  loadings<-loadings[match(rowData(target)$gene_name,rownames(loadings)),]


  A <- assay(target, assay)
  options(RcppML.threads = 0) #line below doesn't work otherwise
  proj<-RcppML::project(data=as.matrix(A), w=loadings, threads=0, L1=0, mask=NULL)

  factors <- t(proj)/nmf_model$d #scale by the constant factor
  colnames(factors) <- paste0("NMF", 1:ncol(factors))

  return(factors)
}
