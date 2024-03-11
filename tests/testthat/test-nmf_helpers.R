test_that("multiplication works", {
  # get the annotated source data
  ehub <- ExperimentHub::ExperimentHub()
  layer_labs <- "BayesSpace_harmony_09"
  data_type <- "spatialDLPFC_Visium"
  vis_anno <- spatialLIBD::fetch_data(type = data_type, eh = ehub)

  vis_anno_sub <- vis_anno[,which(vis_anno$sample_id %in% unique(vis_anno$sample_id)[1:2])]

  # get the unannotated target data
  data(sfe)

  source_nmf_mod <- run_nmf(data=vis_anno_sub, assay="logcounts", k=10,
                            seed=seed)

  source_factors <- t(source_nmf_mod$h)
  colnames(source_factors) <- paste0("NMF", 1:ncol(source_factors))

  annots<- colData(vis_anno_sub)[[layer_labs]]
  factor_annot_cors <- compute_factor_correlations(source_factors, annots)
  factors_use <- identify_factors_representing_annotations(factor_annot_cors)
  print(factors_use)

  factors_use <- source_factors[,factors_use]
  print(head(factors_use))
  multinom_mod <- fit_multinom_model(factors_use, colData(vis_anno_sub)[[layer_labs]])


  sfe <- scuttle::logNormCounts(sfe)
  projections <- project_factors(vis_anno_sub, sfe, assay="logcounts", source_nmf_mod)

  probs <- predict(multinom_mod, newdata=projections, type='probs',
                   na.action=na.exclude)

  preds <- unlist(lapply(1:nrow(probs), function(xx){
    colnames(probs)[which.max(probs[xx,])]
  }))
})
