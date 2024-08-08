plot_factor_correlations <- function(traam.out,...){
  cors <- traam.out$correlations
  ComplexHeatmap::Heatmap(cors,...)
}

plot_nmf_factor_on_tissue <- function(sce, fcts, fct_name, annotationsName){

  colData(sce)[fct_name] <- fcts[fct_name]
  escheR::make_escheR(sce) %>%
    escheR::add_fill(fct_name) %>%
    escheR::add_ground(annotationsName)

}


#' Plot clusters by facet
#'
#' @param spe A SpatialExperiment or SpatialFeatureExperiment object
#' @param annotationsName the name of the `colData(spe)` column that stores the annotations to plot
#'
#' @return A list of `ggplot` objects with each item in the list being a plot of a level in `annotationsName`
#' @export
plot_faceted_clusters <- function(spe, annotationsName){
  annots <- unique(colData(spe)[[annotationsName]])
  annot_plts <- list()
  for (l in 1:length(annots)){
    spe$isAnnot <- (annots[[l]] == colData(spe)[[annotationsName]])
    annot_plts[[l]] <- escheR::make_escheR(spe, y_reverse=FALSE)%>%
      escheR::add_fill("isAnnot", point_size=0.75)+
      ggplot2::scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"))+
      ggplot2::ggtitle(annots[[l]])+
      ggplot2::theme(legend.text=element_text(size=20),
            plot.title=element_text(size=20))
  }
  return(annot_plts)
}
