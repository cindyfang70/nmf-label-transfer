#' @export
transfer_labels.list <- function(targets, source, assay="logcounts", annotationsName, technicalVarName, seed=123, nmf_path="nmf_mod.RDS", save_nmf=TRUE,...) {

  check_source_validity(source, assay, annotationsName)

  for (i in 1:length(targets)){
    check_targets_validity(assay, targets[[i]])
  }
  if(any(is.na(colData(source)[[annotationsName]]))){
    warning("NAs in annotations, removing from source dataset")
    source <- source[,which(!is.na(colData(source)[[annotationsName]]))]
  }
  source_outputs <- source_nmf_and_model_fitting(source, assay, seed,
                                                 save_nmf, nmf_path,
                                                 annotationsName, technicalVarName,...)

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
#' @param targets A list of SingleCellExperiment or SpatialExperiment objects to transfer labels into
#' @param assay Name of the assay in `source` to use for NMF
#' @param annotationsName Name of the annotations in `source` as found in the `colData`
#' @param technicalVarName Name of a technical variable in `source` as found in the `colData`, used to identify NMF factors representing batch effects.
#' @param seed A random seed
#' @param save_nmf TRUE/FALSE specifying whether to save the NMF model to disk
#' @param nmf_path if saving the NMF model to disk, the file path to save it to.
#' @param ... Additional parameters passed to `run_nmf`

#'
#' @return A list of predicted labels for each dataset in `target`
#'
#' @import SpatialExperiment
#' @import SummarizedExperiment
#' @import methods
#' @export
#'
#' @examples
#' library(spatialLIBD)
#' library(SpatialExperiment)
#'
#' # fetch the data
#' ehub <- ExperimentHub::ExperimentHub()
#' layer_labs <- "BayesSpace_harmony_09"
#' data_type <- "spatialDLPFC_Visium"
#' vis_anno <- spatialLIBD::fetch_data(type = data_type, eh = ehub)
#'
#' # select reference datasets
#' vis_anno_ref <- vis_anno[,which(vis_anno$sample_id %in% unique(vis_anno$sample_id)[c(1,2)])]
#'
#' # select a target dataset
#' vis_anno_target <- vis_anno[,which(vis_anno$sample_id %in% unique(vis_anno$sample_id)[3])]
#'
#' # transfer labels
#' res <- transfer_labels(targets=vis_anno_target,
#'             source=vis_anno_ref,
#'             assay="logcounts",
#'             annotationsName=layer_labs,
#'             technicalVarName="sample_id",
#'             save_nmf=FALSE, k=10, tol=1e-5)
#'
#' target_with_preds <- res$targets
#' print(target_with_preds)
#' # display the results
#' table(target_with_preds$nmf_preds)
#'
transfer_labels <- function(targets, source, assay="logcounts", annotationsName, technicalVarName, seed=123, nmf_path="nmf_mod.RDS", save_nmf=TRUE,...){
  UseMethod("transfer_labels")
}

#' @export
transfer_labels.SpatialExperiment <- function(targets, source, assay="logcounts", annotationsName, technicalVarName, seed=123, nmf_path="nmf_mod.RDS", save_nmf=TRUE,...){

  check_source_validity(source, assay, annotationsName)
  check_targets_validity(assay, targets)
  if(any(is.na(colData(source)[[annotationsName]]))){
    warning("NAs in annotations, removing from source dataset")
    source <- source[,which(!is.na(colData(source)[[annotationsName]]))]
  }

  source_outputs <- source_nmf_and_model_fitting(source, assay, seed,
                                                 save_nmf, nmf_path,
                                                 annotationsName,
                                                 technicalVarName,...)

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


