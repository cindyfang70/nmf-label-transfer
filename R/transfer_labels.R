#' @export
transfer_labels.list <- function(targets, source, assay="logcounts", annotationsName, seed=123, save_nmf=TRUE, nmf_path="nmf_mod.RDS",...) {

  check_source_validity(source, assay, annotationsName)

  for (i in 1:length(targets)){
    check_targets_validity(assay, targets[[i]])
  }

  source_outputs <- source_nmf_and_model_fitting(source, assay, seed,
                                                 save_nmf, nmf_path,
                                                 annotationsName,...)

  source_nmf_mod <- source_outputs$source_nmf
  factors_use_names <- source_outputs$factors_use_names
  multinom_mod <- source_outputs$multinom

  for (i in 1:length(targets)){
    target <- targets[[i]]
    # 4: project patterns onto target dataset
    projections <- project_factors(source, target, assay, source_nmf_mod)
    reducedDim(target, "nmf_projections") <- projections
    projections <- projections[,factors_use_names]

    # 5: predict on the projected factors using the multinomial model
    probs <- predict(multinom_mod, newdata=projections, type='probs',
                     na.action=na.exclude)

    preds <- unlist(lapply(1:nrow(probs), function(xx){
      colnames(probs)[which.max(probs[xx,])]
    }))


    colData(target)$nmf_preds <- preds
    targets[[i]] <- target
  }

  traam.out <- traam(nmf_model=source_nmf_mod,
                     source=source,
                     targets=targets,
                     multinom=multinom_mod,
                     cors=source_outputs$factor_annot_cors)

  return(traam.out)
}


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
#' @export
transfer_labels <- function(targets, source, assay="logcounts", annotationsName, seed=123, save_nmf=TRUE,...){
  UseMethod("transfer_labels")
}

#' @export
transfer_labels.SpatialExperiment <- function(targets, source, assay="logcounts", annotationsName, seed=123, save_nmf=TRUE,...){

  check_source_validity(source, assay, annotationsName)
  check_targets_validity(assay, targets)

  source_outputs <- source_nmf_and_model_fitting(source, assay, seed,
                                                 save_nmf, nmf_path,
                                                 annotationsName,...)

  source_nmf_mod <- source_outputs$source_nmf
  factors_use_names <- source_outputs$factors_use_names
  multinom_mod <- source_outputs$multinom

  # 4: project patterns onto target dataset
  projections <- project_factors(source, targets, assay, source_nmf_mod)
  reducedDim(targets, "nmf_projections") <- projections
  projections <- projections[,factors_use_names]

  # 5: predict on the projected factors using the multinomial model
  probs <- predict(multinom_mod, newdata=projections, type='probs',
                     na.action=na.exclude)

  preds <- unlist(lapply(1:nrow(probs), function(xx){
      colnames(probs)[which.max(probs[xx,])]
    }))


  colData(targets)$nmf_preds <- preds


  traam.out <- traam(nmf_model=source_nmf_mod,
                     source=source,
                     targets=targets,
                     multinom=multinom_mod,
                     cors=source_outputs$factor_annot_cors)

  return(traam.out)
}


