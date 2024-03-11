test_that("label transfer works", {
  # get the annotated source data
  ehub <- ExperimentHub::ExperimentHub()
  layer_labs <- "BayesSpace_harmony_09"
  data_type <- "spatialDLPFC_Visium"
  vis_anno <- spatialLIBD::fetch_data(type = data_type, eh = ehub)

  vis_anno_sub <- vis_anno[,which(vis_anno$sample_id %in% unique(vis_anno$sample_id)[1:2])]

  # get the unannotated target data
  data(sfe)

  sfe_annot <- transfer_labels(vis_anno_sub, sfe, "logcounts", layer_labs)
})
