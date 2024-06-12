find_neighbours <- function(spe, k=10){
  spe_coords <- as.data.frame(spatialCoords(spe))
  nns <- FNN::get.knn(data=spe_coords, k=k)
  #print(nns)
  reducedDim(spe, "KNN") <- nns$nn.index
  return(spe)
}


compute_nbhd_factors_mean <- function(spe, factors){
  spe <- find_neighbours(spe)
  nbhd_fcts_mean <- matrix(,nrow=nrow(factors), ncol=ncol(factors))
  for(i in 1:ncol(factors)){
    fct <- factors[,i]
    fct_nbhd_mean <- unlist(lapply(1:ncol(spe), function(xx){
      mean(fct[reducedDim(spe, "KNN")[xx,]])
    }))
    nbhd_fcts_mean[,i] <- fct_nbhd_mean
  }
  colnames(nbhd_fcts_mean) <- paste0(colnames(factors), "nns-mean")
  fcts_with_nbhd_means <- cbind(factors, nbhd_fcts_mean)
  return(fcts_with_nbhd_means)
}
