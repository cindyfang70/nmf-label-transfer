smooth_label <- function(spe, annotationsName, label, n_neighours=10, cutoff=0.2){
  spe_coords <- as.data.frame(spatialCoords(spe))

}


#' @importFrom dbscan kNN
#' this is taken from BANKSY: https://github.com/prabhakarlab/Banksy/blob/bioc/R/cluster.R
#' @export
smoother <-
  function(labels_curr,
           locs,
           k = 10L,
           props = NULL,
           max_iter = 10,
           verbose = TRUE) {

    # Set the proportions to be 0.5 for each class if the user did not specify
    if (is.null(props)){
      props <- rep(0.5, length(unique(labels_curr)))
      names(props) <- names(unique(labels_curr))
    }

    # Get neighbors
    knn <- kNN(locs, k = k)$id
    N <- nrow(knn)

    labels_raw <- labels_curr
    labels_update <- labels_curr

    if (max_iter == -1) max_iter <- Inf

    iter <- 1
    while (iter < max_iter) {
      if (verbose) {
        message("Iteration ", iter)
      }

      # Iterate across cells
      for (i in seq(N)) {
        # Get neighbors
        neighbor_labels <- labels_curr[c(i, knn[i, ])]
        neighbor_props <- table(neighbor_labels) /
          length(neighbor_labels)

        cell_type_i <- labels_curr[[i]]
        prop_thres <- props[cell_type_i]
        # Change label based on condition
        if (any(neighbor_props > prop_thres)) {
          labels_update[i] <- as.numeric(
            names(which.max(neighbor_props))
          )
        } else {
          labels_update[i] <- labels_curr[i]
        }
      }

      change <- sum(labels_update != labels_curr)
      if (verbose) {
        message("Change: ", change)
      }

      if (change == 0) {
        break
      }

      labels_curr <- labels_update
      iter <- iter + 1
    }
    labels_update
  }
