#' Transfer labels from the source dataset to the target dataset
#'
#' @param source A SingleCellExperiment or SpatialExperiment object with cell-type or spatial domain labels
#' @param target A SingleCellExperiment or SpatialExperiment object to transfer labels into
#' @param assay Name of the assay in `source` to use for NMF
#' @param annotationsName Name of the annotations in `source` as found in the `colData`
#' @param seed A random seed
#'
#' @return A SingleCellExperiment or SpatialExperiment object
#'
#' @import SpatialExperiment
#' @import SummarizedExperiment
#' @import methods
#' @export
transfer_labels <- function(source, target, assay, annotationsName, seed=123) {

  if(is(source, "SpatialExperiment")){
    stopifnot(assay %in% assays(source))
    stopifnot(annotationsName %in% colnames(colData(source)))
  }

  # 1: run NMF on the source dataset
  source_nmf_mod <- run_nmf(data=source, assay="logcounts", k=100,
                            seed=seed)

  # 2: select the important factors based on correlation
  annots<- colData(source)[[annotationsName]]
  factor_annot_cors <- compute_factor_correlations(source_nmf_mod, annots)
  factors_use <- identify_factors_representing_annotations(factor_annot_cors)
  # 3: fit multinomial model on source factors
  multinom_mod <- fit_multinom_model(factors_use, annots)
  # 4: project patterns onto target dataset
  projections <- project_factors(source, target, assay, source_nmf_mod)
  # 5: predict annotations on target factors

  probs <- predict(multinom, newdata=projections, type='probs',
                   na.action=na.exclude)

  preds <- unlist(lapply(1:nrow(probs), function(xx){
    colnames(probs)[which.max(probs[xx,])]
  }))

  return(preds)
}

