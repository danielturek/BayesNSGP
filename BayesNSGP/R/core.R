

#' Just created a single function.  Ok to delete.
#' @export
ok_to_delete_this <- function() {
    print('placeholder; just delete this')
}


## Script #1: nsgpCore.R (core bayesNSGP functionality)
## 
## - inverseEigen: calculate covariance elements based on eigendecomposition components
## - nsCorr: calculate a nonstationary Matern correlation matrix
## - nsCrosscorr: calculate a nonstationary Matern cross-correlation matrix
## - nsDist: calculate coordinate-specific distance matrices
## - nsCrossdist: calculate coordinate-specific cross-distance matrices
## - nsDist3d (formerly ns_dist_3d)
## - nsCrossdist3d (TODO)
## - nsgpModel: NIMBLE code for a generic nonstationary GP model 
## - TODO: add SGV likelihood calculation
## - nsgpPredict: posterior prediction for the NSGP TODO
## 
## 
## Script #2: nsgpOrderingNN.R (functions for ordering and finding nearest neighbors)
## 
## - orderCoordinatesMMD: order coordinates by maxmin distance
## - determineNeighbors: identify k nearest neighbors


