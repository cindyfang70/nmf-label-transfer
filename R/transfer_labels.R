#' Transfer labels from the source dataset to the target dataset
#'
#' @param source A SingleCellExperiment or SpatialExperiment object with cell-type or spatial domain labels
#' @param target A list of SingleCellExperiment or SpatialExperiment objects to transfer labels into
#' @param assay Name of the assay in `source` to use for NMF
#' @param annotationsName Name of the annotations in `source` as found in the `colData`
#' @param seed A random seed
#' @param ... Additional parameters passed to `run_nmf`
#'
#' @return A list of predicted labels for each dataset in `target`
#'
#' @import SpatialExperiment
#' @import SummarizedExperiment
#' @import methods
#' @export
transfer_labels <- function(source, target, assay="logcounts", annotationsName, seed=123, save_nmf=TRUE,...) {

  if(is(source, "SpatialExperiment")){
    if (!(assay %in% assayNames(source))){
      stop(sprintf("Assay %s not found in the source dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }
    if(!(annotationsName %in% colnames(colData(source)))){
      stop(sprintf("%s is not found in the colData of the source dataset.", annotationsName))
    }
    if(!("gene_name" %in% colnames(rowData(source)))){
      stop("Please provide gene symbols in your source dataset as a column named 'gene_name' in rowData.")
    }

  }

  for (i in 1:length(target)){
    check_targets_validity(target)
  }



  # 1: run NMF on the source dataset
  source_nmf_mod <- run_nmf(data=source, assay=assay, seed=seed, ...)

  if(save_nmf){
    saveRDS(source_nmf_mod, "source_nmf_mod.RDS")
  }

  source_factors <- t(source_nmf_mod$h)
  colnames(source_factors) <- paste0("NMF", 1:ncol(source_factors))

  # 2: select the important factors based on correlation
  annots <- colData(source)[[annotationsName]]
  factor_annot_cors <- compute_factor_correlations(source_factors, annots)
  factors_use_names <- identify_factors_representing_annotations(factor_annot_cors)


  # 3: fit multinomial model on source factors
  factors_use <- source_factors[,factors_use_names]

  # compute neighbourhood mean of each selected NMF factor
  multinom_mod <- fit_multinom_model(as.data.frame(factors_use), annots)

  all_preds <- list()
  for (i in 1:length(target)){
      preds <- project_and_predict_on_one_target(target[[i]])
      all_preds[[i]] <- preds
  }

  return(all_preds)
}

check_targets_validity <- function(target){
  if (is(target,"SpatialExperiment")){
    if (!(assay %in% assayNames(target))){
      stop(sprintf("Assay %s is not found in the target dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }

    if(!("gene_name" %in% colnames(rowData(target)))){
      stop("Please provide gene symbols in your target dataset as a column named 'gene_name' in rowData.")
    }

  }
}


project_and_predict_on_one_target <- function(target){
  # 4: project patterns onto target dataset
  projections <- project_factors(source, target, assay, source_nmf_mod)

  #print(head(projections))
  projections <- projections[,factors_use_names]
  # 5: predict annotations on target factors

  probs <- predict(multinom_mod, newdata=projections, type='probs',
                   na.action=na.exclude)

  preds <- unlist(lapply(1:nrow(probs), function(xx){
    colnames(probs)[which.max(probs[xx,])]
  }))

  return(preds)
}
