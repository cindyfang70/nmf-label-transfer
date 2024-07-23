test_that("single sample label transfer works", {
  # get the annotated source data
  ehub <- ExperimentHub::ExperimentHub()
  layer_labs <- "BayesSpace_harmony_09"
  data_type <- "spatialDLPFC_Visium"
  vis_anno <- spatialLIBD::fetch_data(type = data_type, eh = ehub)

  vis_anno_sub <- vis_anno[,which(vis_anno$sample_id %in% unique(vis_anno$sample_id)[c(1,2)])]

  # get the unannotated target data
  data(spe)
  rowData(spe)$gene_name <- rownames(spe)
  spe <- scuttle::logNormCounts(spe)

  preds <- transfer_labels(targets=spe, source=vis_anno_sub, assay="logcounts", annotationsName=layer_labs,
                           technicalVarName="sample_id", save_nmf=FALSE, k=10,
                           cv_tol=1e-4, tol=1e-5)

  #print(preds)

  expect_true(length(preds$targets$nmf_preds)==ncol(spe))
})
