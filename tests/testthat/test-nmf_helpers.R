test_that("outputted NMF factors are all nonnegative", {
  # get the annotated source data
  ehub <- ExperimentHub::ExperimentHub()
  layer_labs <- "BayesSpace_harmony_09"
  data_type <- "spatialDLPFC_Visium"
  vis_anno <- spatialLIBD::fetch_data(type = data_type, eh = ehub)

  vis_anno_sub <- vis_anno[,which(vis_anno$sample_id %in% unique(vis_anno$sample_id)[1:2])]

  # get the unannotated target data
  data(sfe)

  source_nmf_mod <- run_nmf(data=vis_anno_sub, assay="logcounts", k=10,
                            seed=123)

  source_factors <- t(source_nmf_mod$h)

  expect_equal(sum(source_factors >= 0), length(source_factors))
})

test_that("predicted probabilities are in [0,1]",{
  # get the annotated source data
  ehub <- ExperimentHub::ExperimentHub()
  layer_labs <- "BayesSpace_harmony_09"
  data_type <- "spatialDLPFC_Visium"
  vis_anno <- spatialLIBD::fetch_data(type = data_type, eh = ehub)

  vis_anno_sub <- vis_anno[,which(vis_anno$sample_id %in% unique(vis_anno$sample_id)[1:2])]

  # get the unannotated target data
  data(spe)
  rowData(spe)$gene_name <- rownames(rowData(spe))

  source_nmf_mod <- run_nmf(data=vis_anno_sub, assay="logcounts", k=10,
                            seed=1237)

  source_factors <- t(source_nmf_mod$h)
  colnames(source_factors) <- paste0("NMF", 1:ncol(source_factors))

  annots<- colData(vis_anno_sub)[[layer_labs]]
  factor_annot_cors <- compute_factor_correlations(source_factors, annots)
  technical_cors <- compute_factor_correlations(source_factors, colData(vis_anno_sub)[["sample_id"]])
  factors_use <- identify_factors_representing_annotations(factor_annot_cors, technical_cors)

  factors_use <- source_factors[,factors_use]
  multinom_mod <- fit_multinom_model(factors_use, colData(vis_anno_sub)[[layer_labs]])


  spe <- scuttle::logNormCounts(spe)
  projections <- project_factors(vis_anno_sub, spe, assay="logcounts", source_nmf_mod)

  probs <- predict(multinom_mod, newdata=projections, type='probs',
                   na.action=na.exclude)

  expect_equal(sum(probs >= 0 & probs <= 1), length(probs))
  expect_equal(as.vector(rowSums(probs)), rep(1, nrow(probs)))
})
