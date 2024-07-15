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
