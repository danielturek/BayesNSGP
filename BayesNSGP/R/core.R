

#================================================
# Bayesian nonstationary Gaussian process 
# modeling in NIMBLE
# Mark Risser and Daniel Turek
# Lawrence Berkeley National Laboratory
# January, 2019
#================================================




#================================================
# Functions for ordering coordinates and finding
# nearest neighbors
#================================================

## Script #2: nsgpOrderingNN.R (functions for ordering and finding nearest neighbors)
## 
## - orderCoordinatesMMD: order coordinates by maxmin distance
## - determineNeighbors: identify k nearest neighbors


#==============================================================================
# Maximum-minimum distance (MMD) coordinate ordering
#==============================================================================

# ROxygen comments ----
#' Order coordinates according to a maximum-minimum distance criterion.
#' 
#' \code{orderCoordinatesMMD} orders an array of (x,y) spatial coordinates 
#' according to the "maximum minimum distance" (MMD), as described in Guinness, 
#' 2018. (Points are selected to maximize their minimum distance to already-
#' selected points).
#' 
#' @param s N x 2 array of N 2-dimensional (x,y) spatial coordinates.
#' @param exact Logical; \code{FALSE} uses a fast approximation to MMD ordering 
#' (and is almost always recommended), while \code{TRUE} uses exact MMD 
#' ordering but is infeasible for large number of locations.
#' 
#' @return A list with two components: (1) an N x 2 array containing the 
#' same spatial coordinates, ordered by MMD, and (2) the same thing, but with 
#' any NA values removed.
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
orderCoordinatesMMD <- function(s, exact = FALSE) {
  ## input s: an Nx2 array of spatial coordinates
  N <- dim(s)[1]
  if(N < 3) return(s)
  if(!exact) {       ## approximate MMD ordering
    initialOrdering <- sample(1:N)
    orderedIndices <- c(initialOrdering, rep(NA, 3*N))
    indexLookupVector <- order(initialOrdering)
    maxNeighbors <- floor(sqrt(N))
    NN <- get.knn(s, k = maxNeighbors)$nn.index
    nextSpot <- N+1
    cycleCheckIndex <- -1
    for(i in 2:(3*N)) {
      (targetIndex <- orderedIndices[i])
      if(cycleCheckIndex == targetIndex) break
      if(cycleCheckIndex == -1) cycleCheckIndex <- targetIndex
      targetNeighbors <- NN[targetIndex, 1:min(maxNeighbors, round(N/(i+N-nextSpot)))]
      targetNeighborLocs <- indexLookupVector[targetNeighbors]
      if(min(targetNeighborLocs) < i) {   ## relocate this index to the back
        orderedIndices[nextSpot] <- targetIndex
        orderedIndices[i] <- NA
        indexLookupVector[targetIndex] <- nextSpot
        nextSpot <- nextSpot + 1
      } else cycleCheckIndex <- -1
    }
    orderedIndicesNoNA <- orderedIndices[!is.na(orderedIndices)]
    orderedS <- s[orderedIndicesNoNA,]
  } else {           ## exact MMD ordering
    availableIndices <- 1:N
    orderedS <- array(NA, c(N,2))
    sbar <- apply(s, 2, mean)   ## group centroid
    iNext <- which.min(sapply(1:N, function(i) sum((s[i,] - sbar)^2)))
    orderedS[1,] <- s[iNext,]
    availableIndices <- setdiff(availableIndices, iNext)
    for(i in 2:N) {
      aIndNext <- which.max(    ## this indexes the availableIndices vector
        sapply(1:(N-i+1), function(j) {
          min(sapply(1:(i-1), function(k) sum((s[availableIndices[j],] - orderedS[k,])^2)))
        }))
      iNext <- availableIndices[aIndNext]   ## this indexes rows of the original s[] array
      orderedS[i,] <- s[iNext,]
      availableIndices <- setdiff(availableIndices, iNext)
    }
    orderedIndicesNoNA <- NULL
  }
  return(list(orderedS = orderedS, orderedIndicesNoNA = orderedIndicesNoNA))
}

#==============================================================================
# Determine the k-nearest neighbors
#==============================================================================

# ROxygen comments ----
#' Determine the k-nearest neighbors for each spatial coordinate.
#' 
#' \code{determineNeighbors} returns an N x k matrix of the nearest neighbors 
#' for spatial locations s, with the ith row giving indices of the k nearest 
#' neighbors to the ith location, which are selected from among the 1,...(i-1) 
#' other spatial locations. The first row is -1's, since the first location has 
#' no neighbors. The i=2 through i=(k+1) rows each necessarily contain 1:i.
#' 
#' @param s N x 2 array of N 2-dimensional (x,y) spatial coordinates.
#' @param k Scalar; number of neighbors
#' 
#' @return An N x k matrix of nearest neighbor indices
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
determineNeighbors <- function(s, k) {
  N <- dim(s)[1]
  d <- dim(s)[2]
  if(k+2 > N) stop()
  nID <- array(-1, c(N,k))     ## populate unused values with -1, to prevent a warning from NIMBLE
  for(i in 2:(k+1))   nID[i, 1:(i-1)] <- as.numeric(1:(i-1))
  if(d == 2){
    for(i in (k+2):N)   nID[i, 1:k] <- as.numeric(order((s[1:(i-1),1] - s[i,1])^2 + (s[1:(i-1),2] - s[i,2])^2)[1:k])
  } else{
    for(i in (k+2):N){
      disti <- 0
      for(j in 1:d){
        disti <- disti + (s[1:(i-1),j] - s[i,j])^2
      }
      nID[i, 1:k] <- as.numeric(order(disti)[1:k])
    }   
  }
  return(nID)
}








#================================================
# Core package functionality
#================================================

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
## - nsgpPredict: posterior prediction for the NSGP TODO

# TODO: 
# -- Code up prediction functions (full GP, NNGP, SGV)
# -- nsCrossdist3d for prediction


#==============================================================================
# Inverse eigendecomposition 
#==============================================================================

# ROxygen comments ----
#' Calculate covariance elements based on eigendecomposition components
#'
#' \code{inverseEigen} calculates the inverse eigendecomposition -- in other
#' words, the covariance elements based on the eigenvalues and vectors (see
#' Paciorek and Schervish, 2006, for details on the parameterization). The 
#' function is coded as a \code{nimbleFunction} (see the \code{nimble} package) 
#' but can also be used as a regular R function.
#' 
#' @param eigen_comp1 N-vector; contains values of the log of the second
#' anisotropy eigenvalue for a set of locations.
#' @param eigen_comp2 N-vector; contains values of the first eigenvector
#' component for a set of locations.
#' @param eigen_comp3 N-vector; contains values of the second eigenvector 
#' component for a set of locations.
#' @param which_Sigma Scalar; one of \code{(1,2,3)}, corresponding to which
#' covariance component should be calculated (Sigma11, Sigma22, or Sigma12,
#' respectively).
#' 
#' @return A correlation matrix for a fixed set of stations and fixed
#' parameter values.
#'
#' @examples
#' # TODO
#'
#' @export
#'
# inverseEigen <- nimbleFunction(
#   run = function( eigen_comp1 = double(1), eigen_comp2 = double(1),
#                   eigen_comp3 = double(1), which_Sigma = double(0) ) {
#     returnType(double(1))
#     
#     Gam_denom <- sqrt( eigen_comp2^2 + eigen_comp3^2 )
#     Gam11 <- eigen_comp2/Gam_denom
#     Gam22 <- eigen_comp2/Gam_denom
#     Gam12 <- -eigen_comp3/Gam_denom
#     Gam21 <- eigen_comp3/Gam_denom
#     
#     Lam2 <- exp(eigen_comp1)
#     Lam1 <- eigen_comp2^2 + eigen_comp3^2 
#     
#     if( which_Sigma == 1 ){ # Return Sigma11
#       return( Gam11^2*Lam1 + Gam12^2*Lam2 )
#     }
#     if( which_Sigma == 2 ){ # Return Sigma22
#       return( Gam21^2*Lam1 + Gam22^2*Lam2 )
#     }
#     if( which_Sigma == 3 ){ # Return Sigma12
#       return( Gam11*Gam21*Lam1 + Gam12*Gam22*Lam2 )
#     }
# 
#     stop('Error in inverseEigen function')  ## prevent compiler warning
#     return(numeric(10))                     ## prevent compiler warning
#     
#   }
# )
#
# Alternatively, parameterize in terms of the two log eigenvalues and rotation parameter
# eigen_comp1 = log(lambda1)
# eigen_comp2 = log(lambda2)
# eigen_comp3 = logit(2*gamma/pi)
inverseEigen <- nimbleFunction(
  run = function( eigen_comp1 = double(1), eigen_comp2 = double(1),
                  eigen_comp3 = double(1), which_Sigma = double(0) ) {
    returnType(double(1))
    
    rotAngle <- (3.141593/2)*exp(eigen_comp3)/(1 + exp(eigen_comp3)) # pi = 3.141593
    Gam11 <- cos(rotAngle)
    Gam22 <- cos(rotAngle)
    Gam12 <- -sin(rotAngle)
    Gam21 <- sin(rotAngle)
    
    Lam1 <- exp(eigen_comp1)
    Lam2 <- exp(eigen_comp2)
    
    if( which_Sigma == 1 ){ # Return Sigma11
      return( Gam11^2*Lam1 + Gam12^2*Lam2 )
    }
    if( which_Sigma == 2 ){ # Return Sigma22
      return( Gam21^2*Lam1 + Gam22^2*Lam2 )
    }
    if( which_Sigma == 3 ){ # Return Sigma12
      return( Gam11*Gam21*Lam1 + Gam12*Gam22*Lam2 )
    }
    
    stop('Error in inverseEigen function')  ## prevent compiler warning
    return(numeric(10))                     ## prevent compiler warning
    
  }
)


#==============================================================================
# Compiled besselK function 
#==============================================================================

# ROxygen comments ----
# Compiled besselK function 
# 
# \code{RbesselK} and \code{CbesselK} calculates the modified Bessel function
# of the third kind.
# 
# @param dst Matrix; contains distances for the besselK function
# @param nu Scalar; smoothness.
# 
# @return A matrix with values of the corresponding Bessel function.
# 
# RbesselK <- nimbleFunction(
#   run = function(dst = double(2), nu = double(0)) {
#     xVector <- besselK(dst, nu)
#     xMatrix <- matrix(xVector, dim(dst)[1], dim(dst)[2])
#     returnType(double(2))
#     return(xMatrix)
#   }
# )


#==============================================================================
# Calculate a nonstationary Matern correlation matrix 
#==============================================================================



# ROxygen comments ----
#' Calculate a nonstationary Matern correlation matrix
#' 
#' \code{nsCorr} calculates a nonstationary correlation matrix for a 
#' fixed set of locations, based on vectors of the unique anisotropy 
#' parameters for each station. Since the correlation function uses a 
#' spatially-varying Mahalanobis distance, this function requires coordinate-
#' specific distance matrices (see below). The function is coded as a 
#' \code{nimbleFunction} (see the \code{nimble} package) but can also be 
#' used as a regular R function.
#' 
#' @param dist1_sq N x N matrix; contains values of pairwise squared distances
#' in the x-coordinate.
#' @param dist2_sq N x N matrix; contains values of pairwise squared distances
#' in the y-coordinate.
#' @param dist12 N x N matrix; contains values of pairwise signed cross-
#' distances between the x- and y-coordinates. The sign of each element is
#' important; see \code{nsDist} function for the details of this calculation.
#' in the x-coordinate.
#' @param Sigma11 Vector of length N; contains the 1-1 element of the 
#' anisotropy process for each station. 
#' @param Sigma22 Vector of length N; contains the 2-2 element of the 
#' anisotropy process for each station. 
#' @param Sigma12 Vector of length N; contains the 1-2 element of the 
#' anisotropy process for each station.
#' @param nu Scalar; Matern smoothness parameter. \code{nu = 0.5} corresponds 
#' to the Exponential correlation; \code{nu = Inf} corresponds to the Gaussian
#' correlation function.
#' @param d TODO
#' 
#' @return A correlation matrix for a fixed set of stations and fixed
#' parameter values.
#' 
#' @examples
#' # TODO
#' 
#' @export
#' 
nsCorr <- nimbleFunction(
  run = function( dist1_sq = double(2), dist2_sq = double(2), dist12 = double(2), 
                  Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1), 
                  nu = double(0), d = double(0) ) {
    
    returnType(double(2))
    N <- length(Sigma11)
    if( dist2_sq[1,1] == -1 ){ # Isotropic case
      # Calculate the scale matrix 
      if(N == 1){
        det1 <- Sigma11^d
        diagSqrtSqrtDet1 <- matrix(sqrt(sqrt(det1)), N, N)
      } else{
        det1 <- Sigma11^d
        diagSqrtSqrtDet1 <- diag(sqrt(sqrt(det1)))
      }
      mat11_a <- matrix(Sigma11, nrow = N, ncol = N)
      mat11 <- 0.5*(mat11_a + t(mat11_a))
      det12 <- mat11^d
      oneOverDet12 <- 1/det12
      Scale.mat <- diagSqrtSqrtDet1 %*% sqrt(oneOverDet12) %*% diagSqrtSqrtDet1
      # Calculate the distance matrix
      Dist.mat <- sqrt( dist1_sq/mat11 )
    } else{
      # Calculate the scale matrix 
      if(N == 1){
        det1 <- Sigma11*Sigma22 - Sigma12^2
        diagSqrtSqrtDet1 <- matrix(sqrt(sqrt(det1)), N, N)
      } else{
        det1 <- Sigma11*Sigma22 - Sigma12^2 
        diagSqrtSqrtDet1 <- diag(sqrt(sqrt(det1)))
      }
      mat11_a <- matrix(Sigma11, nrow = N, ncol = N)
      mat22_a <- matrix(Sigma22, nrow = N, ncol = N)
      mat12_a <- matrix(Sigma12, nrow = N, ncol = N)
      mat11 <- 0.5*(mat11_a + t(mat11_a))
      mat22 <- 0.5*(mat22_a + t(mat22_a))
      mat12 <- 0.5*(mat12_a + t(mat12_a))
      det12 <- mat11*mat22 - mat12^2
      oneOverDet12 <- 1/det12
      Scale.mat <- diagSqrtSqrtDet1 %*% sqrt(oneOverDet12) %*% diagSqrtSqrtDet1
      # Calculate the distance matrix
      inv11 <-  mat22 * oneOverDet12
      inv22 <-  mat11 * oneOverDet12
      inv12 <- -mat12 * oneOverDet12
      Dist.mat <- sqrt( inv11*dist1_sq + 2*inv12*dist12 + inv22*dist2_sq )
    }
    
    # Combine 
    if( nu == 0.5 ){ # Exponential correlation
      Unscl.corr <- exp(-Dist.mat) 
    } else{
      if( nu == Inf ){ # Gaussian (squared exponential) correlation
        Unscl.corr <- exp(-(Dist.mat^2)) 
      } else{ # Else: Matern with smoothness nu
        ##Unscl.corr <- (exp(lgamma(nu)) * 2^(nu - 1))^(-1) * (Dist.mat)^nu * besselK( x = Dist.mat, nu = nu )
        xVector <- besselK(Dist.mat, nu)
        xMatrix <- matrix(xVector, dim(Dist.mat)[1], dim(Dist.mat)[2])
        Unscl.corr <- (exp(lgamma(nu)) * 2^(nu - 1))^(-1) * (Dist.mat)^nu * xMatrix
        diag(Unscl.corr) <- 1
      } 
    }
    nsCorr <- Scale.mat*Unscl.corr
    return(nsCorr)
  }
)



#==============================================================================
# Calculate a stationary Matern correlation matrix 
#==============================================================================

# ROxygen comments ----
#' Calculate a stationary Matern correlation matrix
#' 
#' \code{matern_corr} calculates a stationary Matern correlation matrix for a 
#' fixed set of locations, based on a range and smoothness parameter. This 
#' function is primarily used for the "npGP" and "approxGP" models. The 
#' function is coded as a \code{nimbleFunction} (see the \code{nimble} package) 
#' but can also be used as a regular R function.
#' 
#' @param dist N x N matrix; contains values of pairwise Euclidean distances in 
#' the x-y plane.
#' @param rho Scalar; "range" parameter used to rescale distances
#' @param nu Scalar; Matern smoothness parameter. \code{nu = 0.5} corresponds 
#' to the Exponential correlation; \code{nu = Inf} corresponds to the Gaussian
#' correlation function.
#' 
#' @return A correlation matrix for a fixed set of stations and fixed
#' parameter values.
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
matern_corr <- nimbleFunction(
  run = function( dist = double(2), rho = double(0), nu = double(0) ) {
    returnType(double(2))
    
    Nr <- dim(dist)[1]
    Nc <- dim(dist)[2]
    if( nu == 0.5 ){ # Exponential correlation
      return(exp(-dist/rho))
    }
    if( nu == Inf ){ # Gaussian (squared exponential) correlation
      return(exp(-(dist/rho)^2))
    } 
    
    # Else: Matern with smoothness nu
    ##temp <- (exp(lgamma(nu)) * 2^(nu - 1))^(-1) * (dist/rho)^nu * besselK( x = dist/rho, nu = nu )
    xVector <- besselK(dist/rho, nu)
    xMatrix <- matrix(xVector, dim(dist)[1], dim(dist)[2])
    temp <- (exp(lgamma(nu)) * 2^(nu - 1))^(-1) * (dist/rho)^nu * xMatrix
    if(Nr == Nc){
      diag(temp) <- 1
    }
    return(temp)
  }
)



#==============================================================================
# Calculate a nonstationary Matern cross-correlation matrix 
#==============================================================================

# ROxygen comments ----
#' Calculate a nonstationary Matern cross-correlation matrix
#' 
#' \code{nsCrosscorr} calculates a nonstationary cross-correlation matrix 
#' between two fixed sets of locations (a prediction set with M locations, and
#' the observed set with N locations), based on vectors of the unique anisotropy 
#' parameters for each station. Since the correlation function uses a 
#' spatially-varying Mahalanobis distance, this function requires coordinate-
#' specific distance matrices (see below). The function is coded as a 
#' \code{nimbleFunction} (see the \code{nimble} package) but can also be 
#' used as a regular R function.
#' 
#' @param Xdist1_sq M x N matrix; contains values of pairwise squared cross-distances
#' in the x-coordinate.
#' @param Xdist2_sq M x N matrix; contains values of pairwise squared cross-distances
#' in the y-coordinate.
#' @param Xdist12 M x N matrix; contains values of pairwise signed cross/cross-
#' distances between the x- and y-coordinates. The sign of each element is
#' important; see \code{nsDist} function for the details of this calculation.
#' in the x-coordinate.
#' @param Sigma11 Vector of length N; contains the 1-1 element of the 
#' anisotropy process for each observed location. 
#' @param Sigma22 Vector of length N; contains the 2-2 element of the 
#' anisotropy process for each observed location. 
#' @param Sigma12 Vector of length N; contains the 1-2 element of the 
#' anisotropy process for each observed location.
#' @param PSigma11 Vector of length N; contains the 1-1 element of the 
#' anisotropy process for each prediction location. 
#' @param PSigma22 Vector of length N; contains the 2-2 element of the 
#' anisotropy process for each prediction location. 
#' @param PSigma12 Vector of length N; contains the 1-2 element of the 
#' anisotropy process for each prediction location.
#' @param nu Scalar; Matern smoothness parameter. \code{nu = 0.5} corresponds 
#' to the Exponential correlation; \code{nu = Inf} corresponds to the Gaussian
#' correlation function.
#' @param d TODO
#' 
#' @return A cross-correlation matrix for two fixed sets of stations and fixed
#' parameter values.
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
nsCrosscorr <- nimbleFunction(
  run = function( Xdist1_sq = double(2), Xdist2_sq = double(2), Xdist12 = double(2), 
                  Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1),
                  PSigma11 = double(1), PSigma22 = double(1), PSigma12 = double(1), 
                  nu = double(0), d = double(0) ) {
    
    returnType(double(2))
    N <- length(Sigma11)
    M <- length(PSigma11)
 
    if( Xdist2_sq[1,1] == -1 ){ # Isotropic case
      # Calculate the scale matrix 
      if(N == 1){
        det1 <- Sigma11^d
        diagSqrtSqrtDet1 <- matrix(sqrt(sqrt(det1)), N, N)
      } else{
        det1 <- Sigma11^d 
        diagSqrtSqrtDet1 <- diag(sqrt(sqrt(det1)))
      }
      if(M == 1){
        Pdet1 <- PSigma11^d 
        diagSqrtSqrtPDet1 <- matrix(sqrt(sqrt(Pdet1)), M, M)
      } else{
        Pdet1 <- PSigma11^d
        diagSqrtSqrtPDet1 <- diag(sqrt(sqrt(Pdet1)))
      }
      mat11_1 <- t(matrix(Sigma11, nrow = N, ncol = M))
      mat11_2 <- matrix(PSigma11, nrow = M, ncol = N)
      mat11 <- 0.5*(mat11_1 + mat11_2)
      det12 <- mat11^d
      oneOverDet12 <- 1/det12
      Scale.mat <- diagSqrtSqrtPDet1 %*% sqrt(oneOverDet12) %*% diagSqrtSqrtDet1
      
      # Calculate the distance matrix
      Dist.mat <- sqrt( Xdist1_sq/mat11 )
    } else{
      # Calculate the scale matrix 
      if(N == 1){
        det1 <- Sigma11*Sigma22 - Sigma12^2
        diagSqrtSqrtDet1 <- matrix(sqrt(sqrt(det1)), N, N)
      } else{
        det1 <- Sigma11*Sigma22 - Sigma12^2 
        diagSqrtSqrtDet1 <- diag(sqrt(sqrt(det1)))
      }
      if(M == 1){
        Pdet1 <- PSigma11*PSigma22 - PSigma12^2 
        diagSqrtSqrtPDet1 <- matrix(sqrt(sqrt(Pdet1)), M, M)
      } else{
        Pdet1 <- PSigma11*PSigma22 - PSigma12^2 
        diagSqrtSqrtPDet1 <- diag(sqrt(sqrt(Pdet1)))
      }
      mat11_1 <- t(matrix(Sigma11, nrow = N, ncol = M))
      mat11_2 <- matrix(PSigma11, nrow = M, ncol = N)
      mat22_1 <- t(matrix(Sigma22, nrow = N, ncol = M))
      mat22_2 <- matrix(PSigma22, nrow = M, ncol = N)
      mat12_1 <- t(matrix(Sigma12, nrow = N, ncol = M))
      mat12_2 <- matrix(PSigma12, nrow = M, ncol = N)
      mat11 <- 0.5*(mat11_1 + mat11_2)
      mat22 <- 0.5*(mat22_1 + mat22_2)
      mat12 <- 0.5*(mat12_1 + mat12_2)
      det12 <- mat11*mat22 - mat12^2
      oneOverDet12 <- 1/det12
      Scale.mat <- diagSqrtSqrtPDet1 %*% sqrt(oneOverDet12) %*% diagSqrtSqrtDet1
      
      # Calculate the distance matrix
      inv11 <-  mat22 * oneOverDet12
      inv22 <-  mat11 * oneOverDet12
      inv12 <- -mat12 * oneOverDet12
      Dist.mat <- sqrt( inv11*Xdist1_sq + 2*inv12*Xdist12 + inv22*Xdist2_sq )
    }
 
    # Combine 
    if( nu == 0.5 ){ # Exponential correlation
      Unscl.corr <- exp(-Dist.mat) 
    } else{
      if( nu == Inf ){ # Gaussian (squared exponential) correlation
        Unscl.corr <- exp(-(Dist.mat^2)) 
      } else{ # Else: Matern with smoothness nu
        ##Unscl.corr <- (exp(lgamma(nu)) * 2^(nu - 1))^(-1) * (Dist.mat)^nu * besselK( x = Dist.mat, nu = nu )
        xVector <- besselK(Dist.mat, nu)
        xMatrix <- matrix(xVector, dim(Dist.mat)[1], dim(Dist.mat)[2])
        Unscl.corr <- (exp(lgamma(nu)) * 2^(nu - 1))^(-1) * (Dist.mat)^nu * xMatrix
        ## this line will not fly:  Unscl.corr[Unscl.corr == Inf] <- 1
        ## inelegant, but accomplishes the same:
        if(min(Dist.mat) == 0) {
          for(i in 1:dim(Unscl.corr)[1]) {
            for(j in 1:dim(Unscl.corr)[2]) {
              if(Unscl.corr[i,j] == Inf)   Unscl.corr[i,j] <- 1
            }
          }
        }
        # diag(Unscl.corr) <- 1
      } 
    }
    nsCrosscorr <- Scale.mat*Unscl.corr
    return(nsCrosscorr)
  }
)



#==============================================================================
# Calculate coordinate-specific distance matrices 
#==============================================================================

# ROxygen comments ----
#' Calculate coordinate-specific distance matrices
#' 
#' \code{nsDist} calculates x, y, and x-y distances for use in the 
#' nonstationary correlation calculation. The sign of the cross-distance
#' is important. The function contains an optional argument for re-scaling
#' the distances such that the coordinates lie in a square.
#' 
#' @param coords N x 2 matrix; contains the x-y coordinates of stations
#' @param scale_factor Scalar; optional argument for re-scaling the distances.
#' @param isotropic Logical; indicates whether distances should be calculated
#' separately for each coordinate dimension (FALSE) or simultaneously for all
#' coordinate dimensions (TRUE). \code{isotropic = TRUE} can only be used for
#' two-dimensional coordinate systems.
#' 
#' @return A list of distances matrices, with the following components:
#' \item{dist1_sq}{N x N matrix; contains values of pairwise squared distances
#' in the x-coordinate.}
#' \item{dist2_sq}{N x N matrix; contains values of pairwise squared distances
#' in the y-coordinate.}
#' \item{dist12}{N x N matrix; contains values of pairwise signed cross-
#' distances between the x- and y-coordinates.}
#' \item{scale_factor}{Value of the scale factor used to rescale distances.}
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
nsDist <- function( coords, scale_factor = NULL, isotropic = FALSE ){
  
  N <- nrow(coords)
  d <- ncol(coords)
  if( !isotropic & d != 2 ) stop("Anisotropy (isotropic = FALSE) only available for 2-dimensional coordinate systems.")
  
  if(!isotropic){
    # Calculate distances
    dists1 <- as.matrix(dist(coords[,1], upper = T, diag = T))
    dists2 <- as.matrix(dist(coords[,2], upper = T, diag = T))
    
    temp1 <- matrix(coords[,1], nrow = N, ncol = N) 
    temp2 <- matrix(coords[,2], nrow = N, ncol = N) 
    
    sgn_mat1 <- ( temp1 - t(temp1) >= 0 )
    sgn_mat1[sgn_mat1 == FALSE] <- -1 
    sgn_mat2 <- ( temp2 - t(temp2) >= 0 )
    sgn_mat2[sgn_mat2 == FALSE] <- -1 
    
    dist1_sq <- dists1^2
    dist2_sq <- dists2^2
    dist12 <- sgn_mat1*dists1*sgn_mat2*dists2
  } else{
    dist1_sq <- as.matrix(dist(coords, upper = T, diag = T))^2
    dist2_sq <- matrix(-1, N, N)
    dist12 <- matrix(0, N, N)
  }
  # Rescale if needed
  if( !is.null(scale_factor) ){
    dist1_sq <- dist1_sq/scale_factor  
    dist2_sq <- dist2_sq/scale_factor  
    dist12 <- dist12/scale_factor     
  }   
  
  return(list( 
    dist1_sq = dist1_sq, dist2_sq = dist2_sq, 
    dist12 = dist12, scale_factor = scale_factor ))
}

#==============================================================================
# Coordinate-specific distance matrices, only for NN
#==============================================================================

# ROxygen comments ----
#' Calculate coordinate-specific distance matrices, only for nearest neighbors
#' and store in an array
#' 
#' \code{nsDist3d} generates and returns new 3-dimensional arrays containing
#' the former dist1_sq, dist2_s1, and dist12 matrices, but
#' only as needed for the k nearest-neighbors of each location.
#' these 3D matrices (dist1_3d, dist2_3d, and dist12_3d)
#' are used in the new implementation of calculateAD_ns().
#' 
#' @param coords N x 2 matrix; contains the x-y coordinates of stations.
#' @param nID N x k matrix; contains indices of nearest neighbors.
#' @param scale_factor Scalar; optional argument for re-scaling the distances.
#' @param isotropic Logical; indicates whether distances should be calculated
#' separately for each coordinate dimension (FALSE) or simultaneously for all
#' coordinate dimensions (TRUE). \code{isotropic = TRUE} can only be used for
#' two-dimensional coordinate systems.
#' 
#' @return Arrays with nearest neighbor distances in each coordinate 
#' direction.
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
nsDist3d <- function(coords, nID, scale_factor = NULL, isotropic = FALSE) {
  N <- nrow(coords)
  
  d <- ncol(coords)
  if( !isotropic & d != 2 ) stop("Anisotropy (isotropic = FALSE) only available for 2-dimensional coordinate systems.")
  
  k <- ncol(nID)
  dist1_3d <- array(0, c(N, k+1, k+1))
  dist2_3d <- array(0, c(N, k+1, k+1))
  dist12_3d <- array(0, c(N, k+1, k+1))
  
  if(!isotropic){
    for(i in 2:N) {
      if(i<=k)     nNei <- i-1      else      nNei <- k
      ind <- c( nID[i,1:nNei], i )
      thisN <- nNei + 1
      theseCoords <- coords[ind, ]
      dists1 <- as.matrix(dist(theseCoords[,1]))
      dists2 <- as.matrix(dist(theseCoords[,2]))
      temp1 <- matrix(theseCoords[,1], nrow = thisN, ncol = thisN)
      temp2 <- matrix(theseCoords[,2], nrow = thisN, ncol = thisN)
      sgn_mat1 <- ( temp1 - t(temp1) >= 0 )
      sgn_mat1[sgn_mat1 == FALSE] <- -1
      sgn_mat2 <- ( temp2 - t(temp2) >= 0 )
      sgn_mat2[sgn_mat2 == FALSE] <- -1
      dist1_3d[i, 1:thisN, 1:thisN] <- dists1^2
      dist2_3d[i, 1:thisN, 1:thisN] <- dists2^2
      dist12_3d[i, 1:thisN, 1:thisN] <- sgn_mat1*dists1*sgn_mat2*dists2
    }  
  } else{
    dist2_3d[1,,] <- -1
    for(i in 2:N) {
      if(i<=k)     nNei <- i-1      else      nNei <- k
      ind <- c( nID[i,1:nNei], i )
      thisN <- nNei + 1
      theseCoords <- coords[ind, ]
      dists1 <- as.matrix(dist(theseCoords))
      dist1_3d[i, 1:thisN, 1:thisN] <- dists1^2
      dist2_3d[i,,] <- -1
    }
  }
  
  if(!is.null(scale_factor)) {
    dist1_3d  <- dist1_3d  / scale_factor
    dist2_3d  <- dist2_3d  / scale_factor
    dist12_3d <- dist12_3d / scale_factor
  }
  return(list(dist1_3d  = dist1_3d,
              dist2_3d  = dist2_3d,
              dist12_3d = dist12_3d,
              scale_factor = scale_factor))
}



#==============================================================================
# Calculate coordinate-specific cross-distance matrices 
#==============================================================================

# ROxygen comments ----
#' Calculate coordinate-specific cross-distance matrices
#' 
#' \code{nsCrossdist} calculates coordinate-specific cross distances in x, y,
#' and x-y for use in the nonstationary cross-correlation calculation. This 
#' function is useful for calculating posterior predictions.
#' 
#' @param coords N x 2 matrix; contains x-y coordinates of station (observed)
#' locations.
#' @param Pcoords M x 2 matrix; contains x-y coordinates of prediction
#' locations.
#' @param scale_factor Scalar; optional argument for re-scaling the distances.
#' @param isotropic TODO
#' 
#' @return A list of distances matrices, with the following components:
#' \item{dist1_sq}{M x N matrix; contains values of pairwise squared cross-
#' distances in the x-coordinate.}
#' \item{dist2_sq}{M x N matrix; contains values of pairwise squared cross-
#' distances in the y-coordinate.}
#' \item{dist12}{M x N matrix; contains values of pairwise signed cross-
#' distances between the x- and y-coordinates.}
#' \item{scale_factor}{Value of the scale factor used to rescale distances.}
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
nsCrossdist <- function(coords, Pcoords, scale_factor = NULL, isotropic = FALSE ){
  
  N <- nrow(coords)
  M <- nrow(Pcoords)
  d <- ncol(coords)
  if( !isotropic & d != 2 ) stop("Anisotropy (isotropic = FALSE) only available for 2-dimensional coordinate systems.")
  
  if(!isotropic){
    # Calculate distances
    dists1 <- mahalanobis.dist(data.x = Pcoords[,1], data.y = coords[,1], vc = diag(1))
    dists2 <- mahalanobis.dist(data.x = Pcoords[,2], data.y = coords[,2], vc = diag(1))
    
    temp1a <- matrix(coords[,1], nrow = M, ncol = N, byrow = TRUE) 
    temp1b <- matrix(Pcoords[,1], nrow = M, ncol = N) 
    temp2a <- matrix(coords[,2], nrow = M, ncol = N, byrow = TRUE) 
    temp2b <- matrix(Pcoords[,2], nrow = M, ncol = N) 
    
    sgn_mat1 <- ( temp1a - temp1b >= 0 )
    sgn_mat1[sgn_mat1 == FALSE] <- -1 
    sgn_mat2 <- ( temp2a - temp2b >= 0 )
    sgn_mat2[sgn_mat2 == FALSE] <- -1 
    
    dist1_sq <- dists1^2
    dist2_sq <- dists2^2
    dist12 <- sgn_mat1*dists1*sgn_mat2*dists2
  } else{
    dist1_sq <- mahalanobis.dist(data.x = Pcoords, data.y = coords, vc = diag(d))^2
    dist2_sq <- matrix(-1, M, N)
    dist12 <- matrix(0, M, N)
  }
  # Rescale if needed
  if( !is.null(scale_factor) ){
    dist1_sq <- dist1_sq/scale_factor  
    dist2_sq <- dist2_sq/scale_factor  
    dist12 <- dist12/scale_factor     
  }   
  
  return(list( 
    dist1_sq = dist1_sq, dist2_sq = dist2_sq, 
    dist12 = dist12, scale_factor = scale_factor ))
}   

#==============================================================================
# Coordinate-specific distance matrices, only for NN
#==============================================================================

# ROxygen comments ----
#' Calculate coordinate-specific distance matrices, only for nearest neighbors
#' and store in an array
#' 
#' \code{nsCrossdist3d} generates and returns new 3-dimensional arrays containing
#' the former dist1_sq, dist2_s1, and dist12 matrices, but
#' only as needed for the k nearest-neighbors of each location.
#' these 3D matrices (dist1_3d, dist2_3d, and dist12_3d)
#' are used in the new implementation of calculateAD_ns().
#' 
#' @param coords N x 2 matrix; contains the x-y coordinates of stations.
#' @param predCoords TODO
#' @param P_nID N x k matrix; contains indices of nearest neighbors.
#' @param scale_factor Scalar; optional argument for re-scaling the distances.
#' @param isotropic Logical; indicates whether distances should be calculated
#' separately for each coordinate dimension (FALSE) or simultaneously for all
#' coordinate dimensions (TRUE). \code{isotropic = TRUE} can only be used for
#' two-dimensional coordinate systems.
#' 
#' @return Arrays with nearest neighbor distances in each coordinate 
#' direction.
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
nsCrossdist3d <- function(coords, predCoords, P_nID, scale_factor = NULL, isotropic = FALSE) {
  N <- nrow(coords)
  M <- nrow(predCoords)

  d <- ncol(coords)
  if( !isotropic & d != 2 ) stop("Anisotropy (isotropic = FALSE) only available for 2-dimensional coordinate systems.")
  
  k <- ncol(P_nID)
  thisN <- k+1
  dist1_3d <- array(0, c(M, k+1, k+1))
  dist2_3d <- array(0, c(M, k+1, k+1))
  dist12_3d <- array(0, c(M, k+1, k+1))
  
  if(!isotropic){
    for(i in 1:M){
      theseCoords <- rbind(coords[P_nID[i,], ], predCoords[i,])
      dists1 <- as.matrix(dist(theseCoords[,1]))
      dists2 <- as.matrix(dist(theseCoords[,2]))
      temp1 <- matrix(theseCoords[,1], nrow = thisN, ncol = thisN)
      temp2 <- matrix(theseCoords[,2], nrow = thisN, ncol = thisN)
      sgn_mat1 <- ( temp1 - t(temp1) >= 0 )
      sgn_mat1[sgn_mat1 == FALSE] <- -1
      sgn_mat2 <- ( temp2 - t(temp2) >= 0 )
      sgn_mat2[sgn_mat2 == FALSE] <- -1
      dist1_3d[i, 1:thisN, 1:thisN] <- dists1^2
      dist2_3d[i, 1:thisN, 1:thisN] <- dists2^2
      dist12_3d[i, 1:thisN, 1:thisN] <- sgn_mat1*dists1*sgn_mat2*dists2
    }
  } else{
    for(i in 1:M){
      theseCoords <- rbind(coords[P_nID[i,], ], predCoords[i,])
      dists1 <- as.matrix(dist(theseCoords))
      dist1_3d[i, 1:thisN, 1:thisN] <- dists1^2
      dist2_3d[i,,] <- -1
    }
  }
  
  if(!is.null(scale_factor)) {
    dist1_3d  <- dist1_3d  / scale_factor
    dist2_3d  <- dist2_3d  / scale_factor
    dist12_3d <- dist12_3d / scale_factor
  }
  return(list(dist1_3d  = dist1_3d,
              dist2_3d  = dist2_3d,
              dist12_3d = dist12_3d,
              scale_factor = scale_factor))
}



#==============================================================================
# NIMBLE code for a generic nonstationary GP model
#==============================================================================

# ROxygen comments ----
#' NIMBLE code for a generic nonstationary GP model
#' 
#' TODO: add documentation
#' 
#' @param tau_model Character; specifies the model to be used for the log(tau) 
#' process. Options are \code{"logLinReg"} (log-linear regression), 
#' \code{"mixComp"} (mixture component representation), \code{"GP"} (stationary
#' Gaussian process), and \code{"approxGP"} (approximation to a Gaussian 
#' process). 
#' @param sigma_model Character; specifies the model to be used for the 
#' log(sigma) process. See \code{tau_model} for options.
#' @param Sigma_model Character; specifies the model to be used for the 
#' Sigma anisotropy process. Options are \code{"covReg"} (covariance 
#' regression), \code{"compReg"} (componentwise regression), 
#' \code{"npMixComp"} (nonparameteric regression via the mixture component
#' approach), \code{"npGP"} (nonparameteric regression via a stationary 
#' Gaussian process), or \code{"npApproxGP"} (nonparameteric regression via an
#' approximation to a stationary Gaussian process).
#' @param mu_model TODO
#' @param likelihood TODO
#' @param coords TODO
#' @param z TODO
#' @param constants TODO
#' @param returnModelComponents TODO
#' @param ... TODO
#' 
#' @return A \code{nimbleCode} object.
#' 
#' @examples
#' # TODO
#'
#' @export
#' 
nsgpModel <- function( tau_model   = "constant",
                       sigma_model = "constant",
                       Sigma_model = "constant",
                       mu_model    = "constant",
                       likelihood  = "fullGP",
                       coords,
                       z,
                       constants = list(),
                       returnModelComponents = FALSE,
                       ... ) {
  
  ##============================================
  ## Models for tau
  ##============================================
  tau_model_list <- list(
    constant = list(
      ## 1. tau_HP1          Standard deviation for the log-linear standard deviation
      ## 2. delta            Scalar; represents the standard deviation (constant over the domain)
      ## 3. ones             N-vector of 1's
      code = quote({
        log_tau_vec[1:N] <- log(sqrt(delta))*ones[1:N]
        delta ~ dunif(0, tau_HP1)
      }),
      constants_needed = c("ones", "tau_HP1"),
      inits = list(delta = quote(tau_HP1/10))
    ),
    
    logLinReg = list(
      ## 1. X_tau            N x p_tau design matrix; leading column of 1's with (p_tau - 1) other covariates
      ## 2. tau_HP1          Standard deviation for the log-linear regression coefficients
      ## 3. p_tau            Number of design columns
      ## 4. delta            Vector of length p_tau; represents log-linear regression coefficients
      code = quote({
        log_tau_vec[1:N] <- X_tau[1:N,1:p_tau] %*% delta[1:p_tau]
        for(l in 1:p_tau){
          delta[l] ~ dnorm(0, sd = tau_HP1)
        }
      }),
      constants_needed = c("X_tau", "p_tau", "tau_HP1"),
      inits = list(delta = quote(rep(0, p_tau)))
    ),
    
    # GP = list(
    #   ## 1. tau_HP1         Gaussian process standard deviation 
    #   ## 2. tau_HP2         Gaussian process mean
    #   ## 3. tau_HP3         Gaussian process range
    #   ## 4. tau_HP4         Gaussian process smoothness
    #   ## 5. ones            N-vector of 1's
    #   ## 6. dist            N x N matrix of inter-point Euclidean distances
    #   code = quote({
    #     log_tau_vec[1:N] ~ dmnorm(mean = log_tau_mn[1:N], cov = log_tau_C[1:N,1:N])
    #     log_tau_mn[1:N] <- tauGP_mu*ones[1:N] 
    #     log_tau_C[1:N,1:N] <- tauGP_sigma^2 * matern_corr(dist[1:N,1:N], tauGP_phi, tau_HP2)
    #     
    #     # Hyperparameters
    #     tauGP_mu ~ dnorm(0, sd = tau_HP1)
    #     tauGP_phi ~ dunif(0, tau_HP3) # Range parameter, GP
    #     tauGP_sigma ~ dunif(0, tau_HP4) # SD parameter, GP
    #     
    #   }),
    #   constants_needed = c("ones", "dist", "tau_HP1", "tau_HP2", "tau_HP3", "tau_HP4"),
    #   inits = list(
    #     log_tau_vec = quote(rep(0, N)),
    #     tauGP_mu = quote(0),
    #     tauGP_phi = quote(tau_HP3/2),
    #     tauGP_sigma = quote(tau_HP4/2)
    #   )
    # ),
    
    approxGP = list(
      ## 1. tau_HP1          Gaussian process standard deviation
      ## 2. tau_HP2          Gaussian process mean
      ## 3. tau_HP3          Gaussian process range
      ## 4. tau_HP4          Gaussian process smoothness
      ## 5. ones             N-vector of 1's
      ## 6. tau_cross_dist   N x p_tau matrix of inter-point Euclidean distances, obs. coords vs. knot locations
      ## 7. tau_knot_dist    p_tau x p_tau matrix of inter-point Euclidean distances, knot locations
      ## 8. p_tau            Number of knot locations
      code = quote({
        
        log_tau_vec[1:N][1:N] <- tauGP_mu*ones[1:N] + tauGP_sigma*Pmat_tau[1:N,1:p_tau] %*% w_tau[1:p_tau]
        Pmat_tau[1:N,1:p_tau] <- matern_corr(tau_cross_dist[1:N,1:p_tau], tauGP_phi, tau_HP2)
        Vmat_tau[1:p_tau,1:p_tau] <- matern_corr(tau_knot_dist[1:p_tau,1:p_tau], tauGP_phi, tau_HP2)
        w_tau_mean[1:p_tau] <- 0*ones[1:p_tau]
        w_tau[1:p_tau] ~ dmnorm( mean = w_tau_mean[1:p_tau], prec = Vmat_tau[1:p_tau,1:p_tau] )
        
        # Hyperparameters
        tauGP_mu ~ dnorm(0, sd = tau_HP1)
        tauGP_phi ~ dunif(0, tau_HP3) # Range parameter, GP
        tauGP_sigma ~ dunif(0, tau_HP4) # SD parameter, GP
        
      }),
      constants_needed = c("ones", "tau_cross_dist", "tau_knot_dist", "p_tau", "tau_HP1", "tau_HP2", "tau_HP3", "tau_HP4"),
      inits = list(
        w_tau = quote(rep(0, p_tau)),
        tauGP_mu = quote(0),
        tauGP_phi = quote(tau_HP3/2),
        tauGP_sigma = quote(tau_HP4/2)
      )
      
    )
  )
  
  ##============================================
  ## Models for sigma
  ##============================================
  
  sigma_model_list <- list(
    constant = list(
      ## 1. sigma_HP1        Standard deviation for the log-linear standard deviation
      ## 2. alpha            Scalar; represents the standard deviation (constant over the domain)
      ## 3. ones             N-vector of 1's
      code = quote({
        log_sigma_vec[1:N] <- log(sqrt(alpha))*ones[1:N]
        alpha ~ dunif(0, sigma_HP1)
      }),
      constants_needed = c("ones", "sigma_HP1"),
      inits = list(alpha = quote(sigma_HP1/10))
    ),
    
    logLinReg = list(
      ## 1. X_sigma          N x p_sigma design matrix; leading column of 1's with (p_sigma - 1) other covariates
      ## 2. sigma_HP1        Standard deviation for the log-linear regression coefficients
      ## 3. p_sigma          Number of design columns
      ## 4. alpha            Vector of length p_sigma; represents log-linear regression coefficients
      code = quote({
        log_sigma_vec[1:N] <- X_sigma[1:N,1:p_sigma] %*% alpha[1:p_sigma]
        for(l in 1:p_sigma){
          alpha[l] ~ dnorm(0, sd = sigma_HP1)
        }
      }),
      constants_needed = c("X_sigma", "p_sigma", "sigma_HP1"),
      inits = list(alpha = quote(rep(0, p_sigma)))
    ),
    
    # GP = list(
    #   ## 1. sigma_HP1        Gaussian process standard deviation 
    #   ## 2. sigma_HP2        Gaussian process mean
    #   ## 3. sigma_HP3        Gaussian process range
    #   ## 4. sigma_HP4        Gaussian process smoothness
    #   ## 5. ones             N-vector of 1's
    #   ## 6. dist             N x N matrix of inter-point Euclidean distances
    #   code = quote({
    #     log_sigma_vec[1:N] ~ dmnorm(mean = log_sigma_mn[1:N], cov = log_sigma_C[1:N,1:N])
    #     log_sigma_mn[1:N] <- sigmaGP_mu*ones[1:N] 
    #     log_sigma_C[1:N,1:N] <- sigmaGP_sigma^2 * matern_corr(dist[1:N,1:N], sigmaGP_phi, sigma_HP2)
    #     
    #     # Hyperparameters
    #     sigmaGP_mu ~ dnorm(0, sd = sigma_HP1)
    #     sigmaGP_phi ~ dunif(0, sigma_HP3) # Range parameter, GP
    #     sigmaGP_sigma ~ dunif(0, sigma_HP4) # SD parameter, GP
    #   }),
    #   constants_needed = c("ones", "dist", "sigma_HP1", "sigma_HP2", "sigma_HP3", "sigma_HP4"),
    #   inits = list(
    #     log_sigma_vec = quote(rep(0, N)),
    #     sigmaGP_mu = quote(0),
    #     sigmaGP_phi = quote(sigma_HP3/2),
    #     sigmaGP_sigma = quote(sigma_HP4/2)
    #   )
    # ),
    
    approxGP = list(
      ## 1. sigma_HP1        Gaussian process standard deviation
      ## 2. sigma_HP2        Gaussian process mean
      ## 3. sigma_HP3        Gaussian process range
      ## 4. sigma_HP4        Gaussian process smoothness
      ## 5. ones             N-vector of 1's
      ## 6. sigma_cross_dist N x p_sigma matrix of inter-point Euclidean distances, obs. coords vs. knot locations
      ## 7. sigma_knot_dist  p_sigma x p_sigma matrix of inter-point Euclidean distances, knot locations
      ## 8. p_sigma          Number of knot locations
      code = quote({
        
        log_sigma_vec[1:N][1:N] <- sigmaGP_mu*ones[1:N] + sigmaGP_sigma*Pmat_sigma[1:N,1:p_sigma] %*% w_sigma[1:p_sigma]
        Pmat_sigma[1:N,1:p_sigma] <- matern_corr(sigma_cross_dist[1:N,1:p_sigma], sigmaGP_phi, sigma_HP2)
        Vmat_sigma[1:p_sigma,1:p_sigma] <- matern_corr(sigma_knot_dist[1:p_sigma,1:p_sigma], sigmaGP_phi, sigma_HP2)
        w_sigma_mean[1:p_sigma] <- 0*ones[1:p_sigma]
        w_sigma[1:p_sigma] ~ dmnorm( mean = w_sigma_mean[1:p_sigma], prec = Vmat_sigma[1:p_sigma,1:p_sigma] )
        
        # Hyperparameters
        sigmaGP_mu ~ dnorm(0, sd = sigma_HP1)
        sigmaGP_phi ~ dunif(0, sigma_HP3) # Range parameter, GP
        sigmaGP_sigma ~ dunif(0, sigma_HP4) # SD parameter, GP
      }),
      constants_needed = c("ones", "sigma_cross_dist", "sigma_knot_dist", "p_sigma", "sigma_HP1", "sigma_HP2", "sigma_HP3", "sigma_HP4"),
      inits = list(
        w_sigma = quote(rep(0, p_sigma)),
        sigmaGP_mu = quote(0),
        sigmaGP_phi = quote(sigma_HP3/2),
        sigmaGP_sigma = quote(sigma_HP4/2)
      )
    )
  )
  
  ##============================================
  ## Models for Sigma
  ##============================================
  
  Sigma_model_list <- list(
    
    constant = list(
      ## 1. ones                N-vector of 1's
      ## 2. Sigma_HP1           Upper bound for the eigenvalues
      ## 3. Sigma_coef{1,2,3}   Vectors of length p_Sigma; represents the anisotropy components
      code = quote({
        
        Sigma11[1:N] <- ones[1:N]*(Sigma_coef1*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef2*sin(Sigma_coef3)*sin(Sigma_coef3))
        Sigma22[1:N] <- ones[1:N]*(Sigma_coef2*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef1*sin(Sigma_coef3)*sin(Sigma_coef3))
        Sigma12[1:N] <- ones[1:N]*(Sigma_coef1*cos(Sigma_coef3)*sin(Sigma_coef3) - Sigma_coef2*cos(Sigma_coef3)*sin(Sigma_coef3))
        
        Sigma_coef1 ~ dunif(0, Sigma_HP1) # phi1
        Sigma_coef2 ~ dunif(0, Sigma_HP1) # phi2
        Sigma_coef3 ~ dunif(0, 1.570796)  # eta --> 1.570796 = pi/2
        
      }),
      constants_needed = c("ones", "Sigma_HP1"),
      inits = list(
        Sigma_coef1 = quote(Sigma_HP1/2),
        Sigma_coef2 = quote(Sigma_HP1/2),
        Sigma_coef3 = 0.7853982 # pi/4
      )
    ),
    constantIso = list( # Isotropic version of 
      ## 1. ones                 N-vector of 1's
      ## 2. Sigma_HP1            Standard deviation for the anisotropy components
      ## 3. Sigma_coef{1,2,3}    Vectors of length p_Sigma; represents the anisotropy components
      code = quote({
        Sigma11[1:N] <- ones[1:N]*Sigma_coef1
        Sigma22[1:N] <- ones[1:N]*Sigma_coef1
        Sigma12[1:N] <- ones[1:N]*0
        
        Sigma_coef1 ~ dunif(0, Sigma_HP1) # phi1
      }),
      constants_needed = c("ones", "Sigma_HP1"),
      inits = list( Sigma_coef1 = quote(Sigma_HP1[1]/2) )
    ),
    
    covReg = list(
      code = quote({
        ## 1. X_Sigma                N x p_Sigma design matrix; leading column of 1's with (p_Sigma - 1) other covariates
        ## 2. Sigma_HP1              Standard deviation for the covariance regression coefficients
        ## 3. p_Sigma                Number of design columns
        ## 4. gamma1, gamma2         Vectors of length p_Sigma; represents covariance regression coefficients
        ## 5. psi11, psi22, rho      Baseline covariance regression parameters
        ## 6. Sigma_HP2              Upper bound for the baseline covariance regression variances
        Sigma11[1:N] <- psi11*ones[1:N] + (X_Sigma[1:N,1:p_Sigma] %*% gamma1[1:p_Sigma])^2
        Sigma12[1:N] <- rho*sqrt(psi11*psi22)*ones[1:N] + (X_Sigma[1:N,1:p_Sigma]%*%gamma1[1:p_Sigma])*(X_Sigma[1:N,1:p_Sigma]%*%gamma2[1:p_Sigma])
        Sigma22[1:N] <- psi22*ones[1:N] + (X_Sigma[1:N,1:p_Sigma] %*% gamma2[1:p_Sigma])^2
        psi11 ~ dunif(0, Sigma_HP2[1])
        psi22 ~ dunif(0, Sigma_HP2[1])
        rho ~ dunif(-1, 1)
        for(j in 1:p_Sigma){
          gamma1[j] ~ dnorm(0, sd = Sigma_HP1[1])
          gamma2[j] ~ dnorm(0, sd = Sigma_HP1[1])
        }
        # Constraints: upper limits on eigen_comp1 and eigen_comp2
        constraint1 ~ dconstraint( max(Sigma11[1:N]) < Sigma_HP5 )
        constraint2 ~ dconstraint( max(Sigma22[1:N]) < Sigma_HP5 )
    }),
      constants_needed = c("ones", "X_Sigma", "p_Sigma", "Sigma_HP1", "Sigma_HP2", "Sigma_HP5"),
      inits = list(
        psi11 = quote(Sigma_HP2[1]/2),
        psi22 = quote(Sigma_HP2[1]/2),
        rho = 0,
        gamma1 = quote(rep(0, p_Sigma)),
        gamma2 = quote(rep(0, p_Sigma)),
        constraint1 = 1, constraint2 = 1
      )
    ), 
    compReg = list(
      code = quote({
        ## 1. X_Sigma                N x p_Sigma design matrix; leading column of 1's with (p_Sigma - 1) other covariates
        ## 2. Sigma_HP1              Standard deviation for the component regression coefficients
        ## 3. p_Sigma                Number of design columns
        ## 4. Sigma_coef{1,2,3}      Vectors of length p_Sigma; represents component regression coefficients
        eigen_comp1[1:N] <- X_Sigma[1:N,1:p_Sigma] %*% Sigma_coef1[1:p_Sigma]
        eigen_comp2[1:N] <- X_Sigma[1:N,1:p_Sigma] %*% Sigma_coef2[1:p_Sigma]
        eigen_comp3[1:N] <- X_Sigma[1:N,1:p_Sigma] %*% Sigma_coef3[1:p_Sigma]
        Sigma11[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 1)
        Sigma12[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 3) 
        Sigma22[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 2)
        for(j in 1:p_Sigma){
          Sigma_coef1[j] ~ dnorm(0, sd = Sigma_HP1)
          Sigma_coef2[j] ~ dnorm(0, sd = Sigma_HP1)
          Sigma_coef3[j] ~ dnorm(0, sd = Sigma_HP1)
        }
        # Constraints: upper limits on eigen_comp1 and eigen_comp2
        constraint1 ~ dconstraint( max(eigen_comp1[1:N]) < log(Sigma_HP5) )
        constraint2 ~ dconstraint( max(eigen_comp2[1:N]) < log(Sigma_HP5) )
      }),
      constants_needed = c("X_Sigma", "p_Sigma", "Sigma_HP1", "Sigma_HP5"),
      inits = list(
        Sigma_coef1 = quote(rep(0, p_Sigma)),
        Sigma_coef2 = quote(rep(0, p_Sigma)),
        Sigma_coef3 = quote(rep(0, p_Sigma)),
        constraint1 = 1, constraint2 = 1
      )
    ),
    
    compRegIso = list( # Isotropic version of compReg
      code = quote({
        ## 1. X_Sigma                N x p_Sigma design matrix; leading column of 1's with (p_Sigma - 1) other covariates
        ## 2. Sigma_HP1              Standard deviation for the component regression coefficients
        ## 3. p_Sigma                Number of design columns
        ## 4. Sigma_coef{1,2,3}      Vectors of length p_Sigma; represents component regression coefficients
        eigen_comp1[1:N] <- X_Sigma[1:N,1:p_Sigma] %*% Sigma_coef1[1:p_Sigma]
        Sigma11[1:N] <- exp(eigen_comp1[1:N])
        Sigma22[1:N] <- exp(eigen_comp1[1:N])
        Sigma12[1:N] <- ones[1:N]*0
        for(j in 1:p_Sigma){
          Sigma_coef1[j] ~ dnorm(0, sd = Sigma_HP1)
        }
        # Constraints: upper limits on eigen_comp1
        constraint1 ~ dconstraint( max(eigen_comp1[1:N]) < log(Sigma_HP5) )
      }),
      constants_needed = c("X_Sigma", "p_Sigma", "Sigma_HP1", "Sigma_HP5"),
      inits = list(
        Sigma_coef1 = quote(rep(0, p_Sigma)),
        constraint1 = 1
      )
    ),
    
    # npGP = list( 
    #   code = quote({
    #     ## 1. Sigma_HP1          3-vector; Gaussian process mean
    #     ## 2. Sigma_HP2          3-vector; Gaussian process smoothness
    #     ## 3. ones                   N-vector of 1's
    #     ## 4. dist                   N x N matrix of inter-point Euclidean distances
    #     
    #     Sigma11[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 1)
    #     Sigma12[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 3) 
    #     Sigma22[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 2)
    #     
    #     # GP1
    #     eigen_comp1[1:N] ~ dmnorm(mean = eigen1_mn[1:N], cov = eigen1_C[1:N,1:N])
    #     eigen1_mn[1:N] <- SigmaGP_mu[1]*ones[1:N] 
    #     eigen1_C[1:N,1:N] <- SigmaGP_sigma[1]^2 * matern_corr(dist[1:N,1:N], SigmaGP_phi[1], Sigma_HP2[1])
    #     
    #     # GP2
    #     eigen_comp2[1:N] ~ dmnorm(mean = eigen2_mn[1:N], cov = eigen2_C[1:N,1:N])
    #     eigen2_mn[1:N] <- SigmaGP_mu[1]*ones[1:N] 
    #     eigen2_C[1:N,1:N] <- SigmaGP_sigma[1]^2 * matern_corr(dist[1:N,1:N], SigmaGP_phi[1], Sigma_HP2[2])
    #     
    #     # GP3
    #     eigen_comp3[1:N] ~ dmnorm(mean = eigen3_mn[1:N], cov = eigen3_C[1:N,1:N])
    #     eigen3_mn[1:N] <- SigmaGP_mu[2]*ones[1:N] 
    #     eigen3_C[1:N,1:N] <- SigmaGP_sigma[2]^2 * matern_corr(dist[1:N,1:N], SigmaGP_phi[2], Sigma_HP2[3])
    #     
    #     # Hyperparameters
    #     for(w in 1:2){
    #       SigmaGP_mu[w] ~ dnorm(0, sd = Sigma_HP1[w])
    #       SigmaGP_phi[w] ~ dunif(0, Sigma_HP3[w]) # Range parameter, GP
    #       SigmaGP_sigma[w] ~ dunif(0, Sigma_HP4[w]) # SD parameter, GP
    #     }
    #     
    #   }),
    #   constants_needed = c("ones", "dist", "Sigma_HP1", "Sigma_HP2", "Sigma_HP3", "Sigma_HP4"),    
    #   inits = list(
    #     eigen_comp1 = quote(rep(0,N)),
    #     eigen_comp2 = quote(rep(0,N)),
    #     eigen_comp3 = quote(rep(0,N)),
    #     SigmaGP_mu = quote(rep(0,2)),
    #     SigmaGP_phi = quote(rep(Sigma_HP3/2,2)),
    #     SigmaGP_sigma = quote(rep(Sigma_HP4/2,2))
    #   )
    # ),
    # 
    # npGPIso = list( 
    #   code = quote({
    #     ## 1. Sigma_HP1          3-vector; Gaussian process mean
    #     ## 2. Sigma_HP2          3-vector; Gaussian process smoothness
    #     ## 3. ones                   N-vector of 1's
    #     ## 4. dist                   N x N matrix of inter-point Euclidean distances
    #     
    #     Sigma11[1:N] <- exp(eigen_comp1[1:N])
    #     Sigma22[1:N] <- exp(eigen_comp1[1:N])
    #     Sigma12[1:N] <- ones[1:N]*0
    #     
    #     # GP1
    #     eigen_comp1[1:N] ~ dmnorm(mean = eigen1_mn[1:N], cov = eigen1_C[1:N,1:N])
    #     eigen1_mn[1:N] <- SigmaGP_mu[1]*ones[1:N] 
    #     eigen1_C[1:N,1:N] <- SigmaGP_sigma[1]^2 * matern_corr(dist[1:N,1:N], SigmaGP_phi[1], Sigma_HP2[1])
    #     
    #     # Hyperparameters
    #     for(w in 1){
    #       SigmaGP_mu[w] ~ dnorm(0, sd = Sigma_HP1[w])
    #       SigmaGP_phi[w] ~ dunif(0, Sigma_HP3[w]) # Range parameter, GP
    #       SigmaGP_sigma[w] ~ dunif(0, Sigma_HP4[w]) # SD parameter, GP
    #     }
    #     
    #   }),
    #   constants_needed = c("ones", "dist", "Sigma_HP1", "Sigma_HP2", "Sigma_HP3", "Sigma_HP4"),    
    #   inits = list(
    #     eigen_comp1 = quote(rep(0,N)),
    #     SigmaGP_mu = quote(rep(0,1)),
    #     SigmaGP_phi = quote(rep(Sigma_HP3/2,1)),
    #     SigmaGP_sigma = quote(rep(Sigma_HP4/2,1))
    #   )
    # ),
    
    npApproxGP = list( 
      code = quote({
        ## 1. Sigma_HP1          3-vector; Gaussian process mean
        ## 2. Sigma_HP2          3-vector; Gaussian process smoothness
        ## 5. ones                   N-vector of 1's
        ## 6. dist                   N x N matrix of inter-point Euclidean distances
        ## 7. Sigma_cross_dist       N x p_Sigma matrix of inter-point Euclidean distances, obs. coords vs. knot locations
        ## 8. Sigma_knot_dist        p_Sigma x p_Sigma matrix of inter-point Euclidean distances, knot locations
        ## 9. p_Sigma                Number of knot locations
        
        Sigma11[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 1)
        Sigma12[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 3) 
        Sigma22[1:N] <- inverseEigen(eigen_comp1[1:N], eigen_comp2[1:N], eigen_comp3[1:N], 2)
        
        # approxGP1, approxGP2
        eigen_comp1[1:N] <- SigmaGP_mu[1]*ones[1:N] + SigmaGP_sigma[1] * Pmat12_Sigma[1:N,1:p_Sigma] %*% w1_Sigma[1:p_Sigma]
        eigen_comp2[1:N] <- SigmaGP_mu[1]*ones[1:N] + SigmaGP_sigma[1] * Pmat12_Sigma[1:N,1:p_Sigma] %*% w2_Sigma[1:p_Sigma]
        
        Pmat12_Sigma[1:N,1:p_Sigma] <- matern_corr(Sigma_cross_dist[1:N,1:p_Sigma], SigmaGP_phi[1], Sigma_HP2[1])
        Vmat12_Sigma[1:p_Sigma,1:p_Sigma] <- matern_corr(Sigma_knot_dist[1:p_Sigma,1:p_Sigma], SigmaGP_phi[1], Sigma_HP2[1])
        w12_Sigma_mean[1:p_Sigma] <- 0*ones[1:p_Sigma]
        
        w1_Sigma[1:p_Sigma] ~ dmnorm( mean = w12_Sigma_mean[1:p_Sigma], prec = Vmat12_Sigma[1:p_Sigma,1:p_Sigma] )
        w2_Sigma[1:p_Sigma] ~ dmnorm( mean = w12_Sigma_mean[1:p_Sigma], prec = Vmat12_Sigma[1:p_Sigma,1:p_Sigma] )
        
        # approxGP3
        eigen_comp3[1:N] <- SigmaGP_mu[2]*ones[1:N] + SigmaGP_sigma[2] * Pmat3_Sigma[1:N,1:p_Sigma] %*% w3_Sigma[1:p_Sigma]
        Pmat3_Sigma[1:N,1:p_Sigma] <- matern_corr(Sigma_cross_dist[1:N,1:p_Sigma], SigmaGP_phi[2], Sigma_HP2[2])
        Vmat3_Sigma[1:p_Sigma,1:p_Sigma] <- matern_corr(Sigma_knot_dist[1:p_Sigma,1:p_Sigma], SigmaGP_phi[2], Sigma_HP2[2])
        w3_Sigma_mean[1:p_Sigma] <- 0*ones[1:p_Sigma]
        w3_Sigma[1:p_Sigma] ~ dmnorm( mean = w3_Sigma_mean[1:p_Sigma], prec = Vmat3_Sigma[1:p_Sigma,1:p_Sigma] )
        
        # Hyperparameters
        for(w in 1:2){
          SigmaGP_mu[w] ~ dnorm(0, sd = Sigma_HP1[w])
          SigmaGP_phi[w] ~ dunif(0, Sigma_HP3[w]) # Range parameter, GP
          SigmaGP_sigma[w] ~ dunif(0, Sigma_HP4[w]) # SD parameter, GP
        }
        
        # Constraints: upper limits on eigen_comp1 and eigen_comp2
        constraint1 ~ dconstraint( max(eigen_comp1[1:N]) < log(Sigma_HP5) )
        constraint2 ~ dconstraint( max(eigen_comp2[1:N]) < log(Sigma_HP5) )
        
      }),
      constants_needed = c("ones", "Sigma_HP1", "Sigma_HP2", "Sigma_HP3", "Sigma_HP4",
                           "Sigma_HP5", "Sigma_cross_dist", "Sigma_knot_dist", "p_Sigma"),    
      inits = list(
        w1_Sigma = quote(rep(0,p_Sigma)),
        w2_Sigma = quote(rep(0,p_Sigma)),
        w3_Sigma = quote(rep(0,p_Sigma)),
        SigmaGP_mu = quote(rep(0,2)),
        SigmaGP_phi = quote(rep(Sigma_HP3/2,2)),
        SigmaGP_sigma = quote(rep(Sigma_HP4/2,2)),
        constraint1 = 1,
        constraint2 = 1
      )
      
    ),
    
    npApproxGPIso = list( 
      code = quote({
        ## 1. Sigma_HP1          3-vector; Gaussian process mean
        ## 2. Sigma_HP2          3-vector; Gaussian process smoothness
        ## 5. ones                   N-vector of 1's
        ## 6. dist                   N x N matrix of inter-point Euclidean distances
        ## 7. Sigma_cross_dist       N x p_Sigma matrix of inter-point Euclidean distances, obs. coords vs. knot locations
        ## 8. Sigma_knot_dist        p_Sigma x p_Sigma matrix of inter-point Euclidean distances, knot locations
        ## 9. p_Sigma                Number of knot locations
        
        Sigma11[1:N] <- exp(eigen_comp1[1:N])
        Sigma22[1:N] <- exp(eigen_comp1[1:N])
        Sigma12[1:N] <- ones[1:N]*0
        
        # approxGP1
        eigen_comp1[1:N] <- SigmaGP_mu[1]*ones[1:N] + SigmaGP_sigma[1] * Pmat12_Sigma[1:N,1:p_Sigma] %*% w1_Sigma[1:p_Sigma]
        
        Pmat12_Sigma[1:N,1:p_Sigma] <- matern_corr(Sigma_cross_dist[1:N,1:p_Sigma], SigmaGP_phi[1], Sigma_HP2[1])
        Vmat12_Sigma[1:p_Sigma,1:p_Sigma] <- matern_corr(Sigma_knot_dist[1:p_Sigma,1:p_Sigma], SigmaGP_phi[1], Sigma_HP2[1])
        w12_Sigma_mean[1:p_Sigma] <- 0*ones[1:p_Sigma]
        
        w1_Sigma[1:p_Sigma] ~ dmnorm( mean = w12_Sigma_mean[1:p_Sigma], prec = Vmat12_Sigma[1:p_Sigma,1:p_Sigma] )
        
        # Hyperparameters
        for(w in 1){
          SigmaGP_mu[w] ~ dnorm(0, sd = Sigma_HP1[w])
          SigmaGP_phi[w] ~ dunif(0, Sigma_HP3[w]) # Range parameter, GP
          SigmaGP_sigma[w] ~ dunif(0, Sigma_HP4[w]) # SD parameter, GP
        }
        
        # Constraints: upper limits on eigen_comp1 and eigen_comp2
        constraint1 ~ dconstraint( max(eigen_comp1[1:N]) < log(Sigma_HP5) )
        
      }),
      constants_needed = c("ones", "Sigma_HP1", "Sigma_HP2", "Sigma_HP3", "Sigma_HP4",
                           "Sigma_HP5", "Sigma_cross_dist", "Sigma_knot_dist", "p_Sigma"),    
      inits = list(
        w1_Sigma = quote(rep(0,p_Sigma)),
        SigmaGP_mu = quote(rep(0,1)),
        SigmaGP_phi = quote(rep(Sigma_HP3/2,1)),
        SigmaGP_sigma = quote(rep(Sigma_HP4/2,1)),
        constraint1 = 1
      )
      
    )
  )
  
  ##============================================
  ## Models for mu
  ##============================================
  
  mu_model_list <- list(
    constant = list(
      ## 1. sigma_HP1          Standard deviation for the log-linear standard deviation
      ## 2. alpha                  Scalar; represents log-linear standard deviation (constant over the domain)
      ## 3. ones                   N-vector of 1's
      code = quote({
        mu[1:N] <-beta*ones[1:N]
        beta ~ dnorm(0, sd = mu_HP1)
      }),
      constants_needed = c("ones", "mu_HP1"),
      inits = list(beta = 0)),
    linReg = list(
      ## 1. X_mu                N x p_mu design matrix; leading column of 1's with (p_mu - 1) other covariates
      ## 2. p_mu                Number of design columns
      ## 3. beta                Vector of length p_mu; represents regression coefficients
      code = quote({
        mu[1:N] <- X_mu[1:N,1:p_mu] %*% beta[1:p_mu]
        for(l in 1:p_mu){
          beta[l] ~ dnorm(0, sd = mu_HP1)
        }
      }),
      constants_needed = c("X_mu", "p_mu", "mu_HP1"),
      inits = list(beta = quote(rep(0, p_mu)))),
    zero = list(
      ## 1. zeros               N-vector of 0's
      code = quote({
        mu[1:N] <- zeros[1:N]
      }),
      constants_needed = c("zeros"),
      inits = list()
    )
  )
  
  ##============================================
  ## Models for likelihood
  ##============================================
  
  likelihood_list <- list(
    fullGP = list(
      code = quote({
        Cor[1:N,1:N] <- nsCorr(dist1_sq[1:N,1:N], dist2_sq[1:N,1:N], dist12[1:N,1:N],
                               Sigma11[1:N], Sigma22[1:N], Sigma12[1:N], nu, d)
        sigmaMat[1:N,1:N] <- diag(exp(log_sigma_vec[1:N]))
        Cov[1:N, 1:N] <- sigmaMat[1:N,1:N] %*% Cor[1:N,1:N] %*% sigmaMat[1:N,1:N]
        C[1:N,1:N] <- Cov[1:N, 1:N] + diag(exp(log_tau_vec[1:N])^2)
        z[1:N] ~ dmnorm(mean = mu[1:N], cov = C[1:N,1:N])
      }),
      constants_needed = c("N", "coords", "d", "dist1_sq", "dist2_sq", "dist12", "nu"),                ## keep N, coords, d here
      inits = list()
    ),
    NNGP = list(
      code = quote({
        AD[1:N,1:(k+1)] <- calculateAD_ns(dist1_3d[1:N,1:(k+1),1:(k+1)],
                                          dist2_3d[1:N,1:(k+1),1:(k+1)],
                                          dist12_3d[1:N,1:(k+1),1:(k+1)],
                                          Sigma11[1:N], Sigma22[1:N], Sigma12[1:N],
                                          log_sigma_vec[1:N], log_tau_vec[1:N],
                                          nID[1:N,1:k], N, k, nu, d)
        z[1:N] ~ dmnorm_nngp(mu[1:N], AD[1:N,1:(k+1)], nID[1:N,1:k], N, k)
      }),
      constants_needed = c("N", "coords", "d", "dist1_3d", "dist2_3d", "dist12_3d", "nID", "k", "nu"),    ## keep N, coords, d here
      inits = list()
    ),
    SGV = list(
      code = quote({
        U[1:num_NZ,1:3] <- calculateU_ns( dist1_3d[1:N,1:(k+1),1:(k+1)], 
                                          dist2_3d[1:N,1:(k+1),1:(k+1)],
                                          dist12_3d[1:N,1:(k+1),1:(k+1)],
                                          Sigma11[1:N], Sigma22[1:N], Sigma12[1:N],
                                          log_sigma_vec[1:N], log_tau_vec[1:N], 
                                          nu, nID[1:N,1:k], cond_on_y[1:N,1:k], N, k, d )
        z[1:N] ~ dmnorm_sgv(mu[1:N], U[1:num_NZ,1:3], N, k)
      }),
      constants_needed = c("N", "coords", "d", "dist1_3d", "dist2_3d", "dist12_3d", "nID", "k", "nu", "cond_on_y", "num_NZ"),    ## keep N, coords, d here
      inits = list()
    )
  )
  
  tau_model_list$mixComp <- tau_model_list$logLinReg        ## add duplicated "mixComp" model option
  sigma_model_list$mixComp <- sigma_model_list$logLinReg    ## add duplicated "mixComp" model option
  Sigma_model_list$npMixComp <- Sigma_model_list$compReg    ## add duplicated "npMixComp" model option
  
  if(is.null(  tau_model_list[[  tau_model]])) stop("unknown specification for tau_model")
  if(is.null(sigma_model_list[[sigma_model]])) stop("unknown specification for sigma_model")
  if(is.null(Sigma_model_list[[Sigma_model]])) stop("unknown specification for Sigma_model")
  if(is.null(   mu_model_list[[   mu_model]])) stop("unknown specification for mu_model")
  if(is.null( likelihood_list[[ likelihood]])) stop("unknown specification for likelihood")
  
  model_selections_list <- list(
    tau        = tau_model_list  [[tau_model]],
    sigma      = sigma_model_list[[sigma_model]],
    Sigma      = Sigma_model_list[[Sigma_model]],
    mu         = mu_model_list   [[mu_model]],
    likelihood = likelihood_list [[likelihood]]
  )
  
  ## code
  
  code_template <- quote({
    SIGMA_MODEL             ## Log variance
    TAU_MODEL               ## Log nugget -- gotta respect the nugget
    CAP_SIGMA_MODEL         ## Anisotropy
    MU_MODEL                ## Mean
    LIKELIHOOD_MODEL        ## Likelihood
  })
  
  code <-
    eval(substitute(substitute(
      CODE,
      list(TAU_MODEL        = model_selections_list$tau$code,
           SIGMA_MODEL      = model_selections_list$sigma$code,
           CAP_SIGMA_MODEL  = model_selections_list$Sigma$code,
           CAP_SIGMA_MODEL  = model_selections_list$Sigma$code,
           MU_MODEL         = model_selections_list$mu$code,
           LIKELIHOOD_MODEL = model_selections_list$likelihood$code)),
      list(CODE = code_template)))
  
  if(missing(z)) stop("must provide data as 'z' argument")
  N <- length(z)
  
  if(missing(coords)) stop("must provide 'coords' argument, array of spatial coordinates")
  d <- ncol(coords)

  sd_default <- 100
  mu_default <- 0
  matern_rho_default <- 1
  matern_nu_default <- 0.5     ## Mark: is this 0.5 good default choice for 'nu'?
  
  constants_defaults_list <- list(
    N = N,
    coords = coords,
    d = d,
    zeros = rep(0, N),
    ones = rep(1, N),
    tau_HP1 = sd_default,            ## standard deviation
    tau_HP2 = mu_default,            ## mean
    tau_HP3 = matern_rho_default,    ## matern_corr 'rho' parameter
    tau_HP4 = matern_nu_default,     ## matern_corr 'nu'  parameter
    sigma_HP1 = sd_default,            ## standard deviation
    sigma_HP2 = mu_default,            ## mean
    sigma_HP3 = matern_rho_default,    ## matern_corr 'rho' parameter
    sigma_HP4 = matern_nu_default,     ## matern_corr 'nu'  parameter
    Sigma_HP1 = 10,    ## standard deviation
    Sigma_HP2 = 10,    ## uniform upper bound for covReg 'psi' parameters
    mu_HP1 = sd_default,            ## standard deviation
    nu = matern_nu_default         ## matern_corr 'nu'  parameter
  )

  ## initialize constants_to_use with constants_defaults_list
  constants_to_use <- constants_defaults_list
  
  ## update constants_to_use with those arguments provided via ...
  dotdotdot <- list(...)
  ## make sure all ... arguments were provided with names
  if(length(dotdotdot) > 0 && (is.null(names(dotdotdot)) || any(names(dotdotdot) == ""))) stop("Only named arguemnts should be provided through ... argument")
  constants_to_use[names(dotdotdot)] <- dotdotdot

  ## add 'constants' argument to constants_to_use list:
  ## if provided, make sure 'constants' argument is a named list
  if(!missing(constants)) {
    if(length(constants) > 0 && (is.null(names(constants)) || any(names(constants) == ""))) stop("All elements in constants list argument must be named")
    constants_to_use[names(constants)] <- constants
  }

  ## generate and add dist1_sq, dist2_sq, and dist12 arrays to constants_to_use list
  if(likelihood == 'fullGP') {
      dist_list <- nsDist(coords)
  } else { ## likelihood is NNGP, or SGV:
      if(is.null(constants_to_use$k)) stop(paste0('missing k constants argument for ', likelihood, ' likelihood'))
      nID <- determineNeighbors(coords, constants_to_use$k)
      constants_to_use$nID <- nID
      dist_list <- nsDist(coords, nID)
  }
  constants_to_use[names(dist_list)] <- dist_list

  ## add the following (derived numbers of columns) to constants_to_use:
  ## p_tau:
  if(!is.null(constants_to_use$X_tau)) constants_to_use$p_tau <- ncol(constants_to_use$X_tau)
  if(!is.null(constants_to_use$tau_cross_dist)) constants_to_use$p_tau <- ncol(constants_to_use$tau_cross_dist)
  ## p_sigma:
  if(!is.null(constants_to_use$X_sigma)) constants_to_use$p_sigma <- ncol(constants_to_use$X_sigma)
  if(!is.null(constants_to_use$sigma_cross_dist)) constants_to_use$p_sigma <- ncol(constants_to_use$sigma_cross_dist)
  ## p_Sigma:
  if(!is.null(constants_to_use$X_Sigma)) constants_to_use$p_Sigma <- ncol(constants_to_use$X_Sigma)
  if(!is.null(constants_to_use$Sigma_cross_dist)) constants_to_use$p_Sigma <- ncol(constants_to_use$Sigma_cross_dist)
  ## p_mu:
  if(!is.null(constants_to_use$X_mu)) constants_to_use$p_mu <- ncol(constants_to_use$X_mu)

         
  ## get a vector of all the constants we need for this model
  constants_needed <- unique(unlist(lapply(model_selections_list, function(x) x$constants_needed), use.names = FALSE))
  ## check if we're missing any constants we need, and throw an error if any are missing
  constants_missing <- setdiff(constants_needed, names(constants_to_use))
  if(length(constants_missing) > 0) {
    stop(paste0("Missing values for the following model constants: ",
                paste0(constants_missing, collapse = ", "),
                ".\nThese values should be provided as named arguments, or named elements in the constants list argument"))
}
  
  ## generate the constants list
  constants <- constants_to_use[constants_needed]
  
  ## data
  data <- list(z = z)
  
  ## inits
  inits_uneval <- do.call("c", unname(lapply(model_selections_list, function(x) x$inits)))
  inits <- lapply(inits_uneval, function(x) eval(x, envir = constants))
  
  if(returnModelComponents) return(list(code=code, constants=constants, data=data, inits=inits))
  
  ## generate the "name" for the nimble model object, containing which submodels were used
  thisName <- paste0(
    'tau='       , tau_model  , '_',
    'sigma='     , sigma_model, '_',
    'Sigma='     , Sigma_model, '_',
    'mu='        , mu_model   , '_',
    'likelihood=', likelihood
  )
  
  ## NIMBLE model object
  Rmodel <- nimbleModel(code, constants, data, inits, name = thisName)
  lp <- Rmodel$getLogProb()
  if(is(lp, 'try-error') || is.nan(lp) || is.na(lp) || abs(lp) == Inf) stop('model not properly initialized')

  ## store 'constants' list into Rmodel$isDataEnv
  Rmodel$isDataEnv$.BayesNSGP_constants_list <- constants
  
  return(Rmodel)
}



#==============================================================================
# Posterior prediction for the NSGP
#==============================================================================

# ROxygen comments ----
#' Posterior prediction for the NSGP
#' 
#' \code{nsgpPredict} conducts posterior prediction for MCMC samples generated
#' using nsgpModel.
#' 
#' @param model TODO
#' @param samples TODO
#' @param coords.predict TODO
#' @param predict.y TODO
#' @param constants TODO
#' @param seed TODO
#' @param ... TODO
#' 
#' @return  TODO
#' 
#' @examples
#' # TODO
#' 
#' @export
#' 
nsgpPredict <- function(model, samples, coords.predict, predict.y = TRUE, constants, seed = 0, ... ) {

  if(!is.model(model)) stop('first argument must be NSGP NIMBLE model object')
  Rmodel <- if(is.Rmodel(model)) model else model$Rmodel
  if(!is.Rmodel(Rmodel)) stop('something went wrong')
  
  model_constants <- Rmodel$isDataEnv$.BayesNSGP_constants_list
  coords <- model_constants$coords
  
  mcmc_samples <- samples
  predCoords <- coords.predict

  ## extract the "submodel" information from the nimble model object "name"
  thisName <- Rmodel$getModelDef()$name
  modelsList <- lapply(strsplit(thisName, '_')[[1]], function(x) strsplit(x, '=')[[1]][2])
  names(modelsList) <- sapply(strsplit(thisName, '_')[[1]], function(x) strsplit(x, '=')[[1]][1])
  ## avaialble for use:
  ## modelsList$tau
  ## modelsList$sigma
  ## modelsList$Sigma
  ## modelsList$mu
  ## modelsList$likelihood


  constants_to_use <- list()
  ## update constants_to_use with those arguments provided via ...
  dotdotdot <- list(...)
  ## make sure all ... arguments were provided with names
  if(length(dotdotdot) > 0 && (is.null(names(dotdotdot)) || any(names(dotdotdot) == ""))) stop("Only named arguemnts should be provided through ... argument")
  constants_to_use[names(dotdotdot)] <- dotdotdot
  ## add 'constants' argument to constants_to_use list:
  ## if provided, make sure 'constants' argument is a named list
  if(!missing(constants)) {
    if(length(constants) > 0 && (is.null(names(constants)) || any(names(constants) == ""))) stop("All elements in constants list argument must be named")
    constants_to_use[names(constants)] <- constants
  }
  ## check for discrepancies in any duplicates, between
  ## constants_to_use provided here, and model_constants from Rmodel
  duplicatedNames <- intersect(names(constants_to_use), names(model_constants))
  discrepancies <- character()
  if(length(duplicatedNames) > 0) {
      for(name in duplicatedNames) {
          if(!identical(constants_to_use[[name]], model_constants[[name]]))
              discrepancies <- c(discrepancies, name)
      }
  }
  if(length(discrepancies) > 0) stop(paste0('Inconsistent values were provided for the following constants: ', paste0(discrepancies, collapse = ', ')))
  ## move original model_constants from nsgpModel into constants_to_use:
  constants_to_use[names(model_constants)] <- model_constants
  ## determine the constants needed
  predictConstantsNeeded <- list(
      tau = list(
          constant = character(),
          logLinReg = c('X_tau', 'PX_tau'),
          approxGP = c('p_tau', 'tau_cross_dist_obs', 'tau_cross_dist_pred', 'tau_HP2')),
      sigma = list(
          constant = character(),
          logLinReg = c('X_sigma', 'PX_sigma'),
          approxGP = c('p_sigma', 'sigma_cross_dist_obs', 'sigma_cross_dist_pred', 'sigma_HP2')),
      Sigma = list(
          constant = character(),
          constantIso = character(),
          covReg = c('X_Sigma', 'PX_Sigma'),
          compReg = c('X_Sigma', 'PX_Sigma'),
          compRegIso = c('X_Sigma', 'PX_Sigma'),
          npApproxGP = c('p_Sigma', 'Sigma_cross_dist_obs', 'Sigma_cross_dist_pred', 'Sigma_HP2'),
          npApproxGPIso = c('p_Sigma', 'Sigma_cross_dist_obs', 'Sigma_cross_dist_pred', 'Sigma_HP2')),
      mu = list(
          constant = character(),
          linReg = c('X_mu', 'PX_mu'),
          zero = character()),
      likelihood = list(
          fullGP = character(),
          NNGP = character(),
          SGV = character())
  )
  constants_needed <- unique(unlist(lapply(1:length(modelsList), function(i) predictConstantsNeeded[[names(modelsList)[i]]][[modelsList[[i]]]] )))
  ## check if we're missing any constants we need, and throw an error if any are missing
  constants_missing <- setdiff(constants_needed, names(constants_to_use))
  if(length(constants_missing) > 0) {
    stop(paste0("Missing values for the following model constants: ",
                paste0(constants_missing, collapse = ", "),
                ".\nThese values should be provided as named arguments, or named elements in the constants list argument"))
}
  ## generate the constants list
  ## do NOT truncate constants list like this:
  ##constants <- constants_to_use[constants_needed]
  constants <- constants_to_use
  
  
  d <- constants$d
  if( modelsList$likelihood == "fullGP" ){ # Predictions for the full GP likelihood
    # Extract needed variables from constants
    dist1_sq <- constants$dist1_sq
    dist2_sq <- constants$dist2_sq
    dist12 <- constants$dist12
    N <- constants$N # number of observed locations
    nu <- constants$nu
    # Prediction distances
    if(dist2_sq[1,1] == -1){ # Isotropic
      Pdist <- nsDist(predCoords, isotropic = TRUE)
      Xdist <- nsCrossdist(coords, predCoords, isotropic = TRUE)
    } else{
      Pdist <- nsDist(predCoords)
      Xdist <- nsCrossdist(coords, predCoords)
    }
    Pdist1_sq <- Pdist$dist1_sq
    Pdist2_sq <- Pdist$dist2_sq
    Pdist12 <- Pdist$dist12
    M <- nrow(Pdist1_sq) # number of prediction locations
    # Cross distances
    Xdist1_sq <- Xdist$dist1_sq
    Xdist2_sq <- Xdist$dist2_sq
    Xdist12 <- Xdist$dist12
    postPredDrawsCols <- M
  } else if( modelsList$likelihood == "NNGP" ){ # Predictions for the NNGP likelihood
    # "Local kriging" only possible for NNGP
    # Extract needed variables from constants
    dist1_3d <- constants$dist1_3d
    dist2_3d <- constants$dist2_3d
    dist12_3d <- constants$dist12_3d
    N <- constants$N # number of observed locations
    k <- constants$k # number of neighbors
    nu <- constants$nu
    # Prediction/cross distances
    P_nID <- get.knnx(coords, predCoords, k = k)$nn.index # Prediction NN
    if(dist2_3d[1,1,1] == -1){
      Pdist <- nsCrossdist3d(coords, predCoords, P_nID, isotropic = TRUE)
    } else{
      Pdist <- nsCrossdist3d(coords, predCoords, P_nID)
    }
    Pdist1_3d <- Pdist$dist1_3d
    Pdist2_3d <- Pdist$dist2_3d
    Pdist12_3d <- Pdist$dist12_3d
    M <- dim(Pdist1_3d)[1] # number of prediction locations
    postPredDrawsCols <- M
  } else if( modelsList$likelihood == "SGV" ){ # Predictions for the SGV likelihood
    if(predict.y == FALSE){   ## TODO: is this right??  when predict.y=FALSE, tell the user that "predicton for Z is not available" ??   -DT
      stop("Prediction for Z(.) not available with SGV.")
    }
    # Extract needed variables from constants
    z <- Rmodel$z
    dist1_3d <- constants$dist1_3d
    dist2_3d <- constants$dist2_3d
    dist12_3d <- constants$dist12_3d
    N <- constants$N # number of observed locations
    k <- constants$k # number of neighbors
    nu <- constants$nu
    # Prediction setup
    predSGV_setup <- sgvSetup(coords, predCoords, k, seed = seed)##SGV_setup$seed)
    prednID_SGV <- predSGV_setup$nID_ord
    if(dist2_3d[1,1,1] == -1){
      preddist_SGV <- nsDist3d( predSGV_setup$locs_ord, prednID_SGV, isotropic = TRUE )
    } else{
      preddist_SGV <- nsDist3d( predSGV_setup$locs_ord, prednID_SGV )
    }
    Alldist1_3d <- preddist_SGV$dist1_3d
    Alldist2_3d <- preddist_SGV$dist2_3d
    Alldist12_3d <- preddist_SGV$dist12_3d
    M <- dim(predCoords)[1] # number of prediction locations
    # REORDER INPUTS (for model != "constant") =================
    obs_ord <- predSGV_setup$ord
    pred_ord <- predSGV_setup$ord_pred
    if( modelsList$tau == "logLinReg" ){
      X_tau_ord <- constants$X_tau[obs_ord,]
      PX_tau_ord <- constants$PX_tau[pred_ord,]
    }
    if( modelsList$tau == "approxGP" ){
      tau_cross_dist_obs_ord <- constants$tau_cross_dist_obs[obs_ord,]
      tau_cross_dist_pred_ord <- constants$tau_cross_dist_pred[pred_ord,]
    }
    if( modelsList$sigma == "logLinReg" ){
      X_sigma_ord <- constants$X_sigma[obs_ord,]
      PX_sigma_ord <- constants$PX_sigma[pred_ord,]
    }
    if( modelsList$sigma == "approxGP" ){
      sigma_cross_dist_obs_ord <- constants$sigma_cross_dist_obs[obs_ord,]
      sigma_cross_dist_pred_ord <- constants$sigma_cross_dist_pred[pred_ord,]
    }
    if( modelsList$Sigma == "covReg" | modelsList$Sigma == "compReg" ){
      X_Sigma_ord <- constants$X_Sigma[obs_ord,]
      PX_Sigma_ord <- constants$PX_Sigma[pred_ord,]
    }
    if( modelsList$Sigma == "npApproxGP" ){
      Sigma_cross_dist_obs_ord <- constants$Sigma_cross_dist_obs[obs_ord,]
      Sigma_cross_dist_pred_ord <- constants$Sigma_cross_dist_pred[pred_ord,]
    }
    if( modelsList$mu == "linReg" ){
      X_mu_ord <- constants$X_mu[obs_ord,]
      PX_mu_ord <- constants$PX_mu[pred_ord,]
    }
    postPredDrawsCols <- M+N
  } else stop('')
  
  J <- nrow(mcmc_samples)
  
  # Posterior draws - storage
  postPredDraws <- matrix(NA, nrow = J, ncol = postPredDrawsCols)
  
  cat("|-------------|-------------|-------------|-------------|")
  cat("\n|")
  for(j in 1:J){ # Loop over MCMC samples
    
    samp_j <- mcmc_samples[j,] 
    
    # Calculate log_tau_vec and Plog_tau_vec ======================
    if( modelsList$tau == "constant" ){
      # Required constants: none
      log_tau_vec_j <- log(sqrt(samp_j["delta"]))*rep(1,N)
      Plog_tau_vec_j <- log(sqrt(samp_j["delta"]))*rep(1,M)
    }
    if( modelsList$tau == "logLinReg" ){
      # Required constants: X_tau, PX_tau
      X_tau <- constants$X_tau
      PX_tau <- constants$PX_tau
      log_tau_vec_j <- as.numeric(X_tau %*% samp_j[paste("delta[",1:ncol(X_tau),"]",sep = "")])
      Plog_tau_vec_j <- as.numeric(PX_tau %*% samp_j[paste("delta[",1:ncol(PX_tau),"]",sep = "")])
    }
    if( modelsList$tau == "approxGP" ){
      # Required constants: p_tau, tau_cross_dist_obs, tau_cross_dist_pred, tau_HP2
      p_tau <- constants$p_tau
      tau_cross_dist_obs <- constants$tau_cross_dist_obs
      tau_cross_dist_pred <- constants$tau_cross_dist_pred
      tau_HP2 <- constants$tau_HP2
      w_tau_j <- as.numeric(samp_j[paste("w_tau[",1:p_tau,"]",sep = "")])
      
      # Obs locations
      Pmat_tau_obs_j <- matern_corr(tau_cross_dist_obs, samp_j["tauGP_phi"], tau_HP2)
      log_tau_vec_j <- as.numeric(samp_j["tauGP_mu"]*rep(1,N) + samp_j["tauGP_sigma"] * Pmat_tau_obs_j %*% w_tau_j)
      
      # Pred locations
      Pmat_tau_pred_j <- matern_corr(tau_cross_dist_pred, samp_j["tauGP_phi"], tau_HP2)
      Plog_tau_vec_j <- as.numeric(samp_j["tauGP_mu"]*rep(1,M) + samp_j["tauGP_sigma"] * Pmat_tau_pred_j %*% w_tau_j)
    }
    
    # Calculate log_sigma_vec and Plog_sigma_vec ==================
    if( modelsList$sigma == "constant" ){
      # Required constants: none
      log_sigma_vec_j <- log(sqrt(samp_j["alpha"]))*rep(1,N)
      Plog_sigma_vec_j <- log(sqrt(samp_j["alpha"]))*rep(1,M)
    }
    if( modelsList$sigma == "logLinReg" ){
      # Required constants: X_sigma, PX_sigma
      X_sigma <- constants$X_sigma
      PX_sigma <- constants$PX_sigma
      log_sigma_vec_j <- as.numeric(X_sigma %*% samp_j[paste("alpha[",1:ncol(X_sigma),"]",sep = "")])
      Plog_sigma_vec_j <- as.numeric(PX_sigma %*% samp_j[paste("alpha[",1:ncol(PX_sigma),"]",sep = "")])
    }
    if( modelsList$sigma == "approxGP" ){
      # Required constants: p_sigma, sigma_cross_dist_obs, sigma_cross_dist_pred, sigma_HP2
      p_sigma <- constants$p_sigma
      sigma_cross_dist_obs <- constants$sigma_cross_dist_obs
      sigma_cross_dist_pred <- constants$sigma_cross_dist_pred
      sigma_HP2 <- constants$sigma_HP2
      w_sigma_j <- as.numeric(samp_j[paste("w_sigma[",1:p_sigma,"]",sep = "")])
      
      # Obs locations
      Pmat_sigma_obs_j <- matern_corr(sigma_cross_dist_obs, samp_j["sigmaGP_phi"], sigma_HP2)
      log_sigma_vec_j <- as.numeric(samp_j["sigmaGP_mu"]*rep(1,N) + samp_j["sigmaGP_sigma"] * Pmat_sigma_obs_j %*% w_sigma_j)
      
      # Pred locations
      Pmat_sigma_pred_j <- matern_corr(sigma_cross_dist_pred, samp_j["sigmaGP_phi"], sigma_HP2)
      Plog_sigma_vec_j <- as.numeric(samp_j["sigmaGP_mu"]*rep(1,M) + samp_j["sigmaGP_sigma"] * Pmat_sigma_pred_j %*% w_sigma_j)
    }
    
    # Calculate SigmaXX and PSigmaXX ==============================
    if( modelsList$Sigma == "constant" ){
      # Required constants: none
      Sigma_coef1 <- samp_j["Sigma_coef1"]
      Sigma_coef2 <- samp_j["Sigma_coef2"]
      Sigma_coef3 <- samp_j["Sigma_coef3"]
      Sigma11_j <- rep(1,N)*(Sigma_coef1*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef2*sin(Sigma_coef3)*sin(Sigma_coef3))
      Sigma22_j <- rep(1,N)*(Sigma_coef2*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef1*sin(Sigma_coef3)*sin(Sigma_coef3))
      Sigma12_j <- rep(1,N)*(Sigma_coef1*cos(Sigma_coef3)*sin(Sigma_coef3) - Sigma_coef2*cos(Sigma_coef3)*sin(Sigma_coef3))
      PSigma11_j <- rep(1,M)*(Sigma_coef1*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef2*sin(Sigma_coef3)*sin(Sigma_coef3))
      PSigma22_j <- rep(1,M)*(Sigma_coef2*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef1*sin(Sigma_coef3)*sin(Sigma_coef3))
      PSigma12_j <- rep(1,M)*(Sigma_coef1*cos(Sigma_coef3)*sin(Sigma_coef3) - Sigma_coef2*cos(Sigma_coef3)*sin(Sigma_coef3))
    }
    if( modelsList$Sigma == "constantIso" ){
      # Required constants: none
      Sigma_coef1 <- samp_j["Sigma_coef1"]
      Sigma_coef2 <- samp_j["Sigma_coef1"]
      Sigma_coef3 <- 0
      Sigma11_j <- rep(1,N)*(Sigma_coef1*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef2*sin(Sigma_coef3)*sin(Sigma_coef3))
      Sigma22_j <- rep(1,N)*(Sigma_coef2*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef1*sin(Sigma_coef3)*sin(Sigma_coef3))
      Sigma12_j <- rep(1,N)*(Sigma_coef1*cos(Sigma_coef3)*sin(Sigma_coef3) - Sigma_coef2*cos(Sigma_coef3)*sin(Sigma_coef3))
      PSigma11_j <- rep(1,M)*(Sigma_coef1*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef2*sin(Sigma_coef3)*sin(Sigma_coef3))
      PSigma22_j <- rep(1,M)*(Sigma_coef2*cos(Sigma_coef3)*cos(Sigma_coef3) + Sigma_coef1*sin(Sigma_coef3)*sin(Sigma_coef3))
      PSigma12_j <- rep(1,M)*(Sigma_coef1*cos(Sigma_coef3)*sin(Sigma_coef3) - Sigma_coef2*cos(Sigma_coef3)*sin(Sigma_coef3))
    }
    if( modelsList$Sigma == "covReg" ){
      # Required constants: X_Sigma, PX_Sigma
      X_Sigma <- constants$X_Sigma
      PX_Sigma <- constants$PX_Sigma
      Sigma11_j <- as.numeric(samp_j["psi11"]*rep(1,N) + (X_Sigma %*% samp_j[paste("gamma1[",1:ncol(X_Sigma),"]",sep = "")])^2)
      Sigma12_j <- as.numeric(samp_j["rho"]*sqrt(samp_j["psi11"]*samp_j["psi22"])*rep(1,N) + (X_Sigma %*% samp_j[paste("gamma1[",1:ncol(X_Sigma),"]",sep = "")])*(X_Sigma %*% samp_j[paste("gamma2[",1:ncol(X_Sigma),"]",sep = "")]))
      Sigma22_j <- as.numeric(samp_j["psi22"]*rep(1,N) + (X_Sigma %*% samp_j[paste("gamma2[",1:ncol(X_Sigma),"]",sep = "")])^2)
      PSigma11_j <- as.numeric(samp_j["psi11"]*rep(1,M) + (PX_Sigma %*% samp_j[paste("gamma1[",1:ncol(PX_Sigma),"]",sep = "")])^2)
      PSigma12_j <- as.numeric(samp_j["rho"]*sqrt(samp_j["psi11"]*samp_j["psi22"])*rep(1,M) + (PX_Sigma %*% samp_j[paste("gamma1[",1:ncol(PX_Sigma),"]",sep = "")])*(PX_Sigma %*% samp_j[paste("gamma2[",1:ncol(PX_Sigma),"]",sep = "")]))
      PSigma22_j <- as.numeric(samp_j["psi22"]*rep(1,M) + (PX_Sigma %*% samp_j[paste("gamma2[",1:ncol(PX_Sigma),"]",sep = "")])^2)
    }
    if( modelsList$Sigma == "compReg" ){
      # Required constants: X_Sigma, PX_Sigma
      X_Sigma <- constants$X_Sigma
      PX_Sigma <- constants$PX_Sigma
      eigen_comp1_j <- X_Sigma %*% samp_j[paste("Sigma_coef1[",1:ncol(X_Sigma),"]",sep = "")]
      eigen_comp2_j <- X_Sigma %*% samp_j[paste("Sigma_coef2[",1:ncol(X_Sigma),"]",sep = "")]
      eigen_comp3_j <- X_Sigma %*% samp_j[paste("Sigma_coef3[",1:ncol(X_Sigma),"]",sep = "")]
      Sigma11_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 1))
      Sigma12_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 3))
      Sigma22_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 2))
      
      Peigen_comp1_j <- PX_Sigma %*% samp_j[paste("Sigma_coef1[",1:ncol(X_Sigma),"]",sep = "")]
      Peigen_comp2_j <- PX_Sigma %*% samp_j[paste("Sigma_coef2[",1:ncol(X_Sigma),"]",sep = "")]
      Peigen_comp3_j <- PX_Sigma %*% samp_j[paste("Sigma_coef3[",1:ncol(X_Sigma),"]",sep = "")]
      PSigma11_j <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 1))
      PSigma12_j <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 3)) 
      PSigma22_j <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 2))
    }
    if( modelsList$Sigma == "compRegIso" ){
      # Required constants: X_Sigma, PX_Sigma
      X_Sigma <- constants$X_Sigma
      PX_Sigma <- constants$PX_Sigma
      eigen_comp1_j <- X_Sigma %*% samp_j[paste("Sigma_coef1[",1:ncol(X_Sigma),"]",sep = "")]
      Sigma11_j <- as.numeric(exp(eigen_comp1_j))
      Sigma12_j <- as.numeric(exp(eigen_comp1_j))
      Sigma22_j <- rep(0,N)
 
      Peigen_comp1_j <- PX_Sigma %*% samp_j[paste("Sigma_coef1[",1:ncol(X_Sigma),"]",sep = "")]
      PSigma11_j <- as.numeric(exp(Peigen_comp1_j))
      PSigma12_j <- as.numeric(exp(Peigen_comp1_j))
      PSigma22_j <- rep(0,N)
    }
    if( modelsList$Sigma == "npApproxGP" ){
      # Required constants: p_Sigma, Sigma_cross_dist_obs, Sigma_cross_dist_pred, Sigma_HP2
      p_Sigma <- constants$p_Sigma
      Sigma_cross_dist_obs <- constants$Sigma_cross_dist_obs
      Sigma_cross_dist_pred <- constants$Sigma_cross_dist_pred
      Sigma_HP2 <- constants$Sigma_HP2
      w1_Sigma_j <- samp_j[paste("w1_Sigma[",1:p_Sigma,"]",sep = "")]
      w2_Sigma_j <- samp_j[paste("w2_Sigma[",1:p_Sigma,"]",sep = "")]
      w3_Sigma_j <- samp_j[paste("w3_Sigma[",1:p_Sigma,"]",sep = "")]
      
      # Obs locations
      Pmat12_Sigma_obs_j <- matern_corr(Sigma_cross_dist_obs, samp_j["SigmaGP_phi[1]"], Sigma_HP2[1])
      Pmat3_Sigma_obs_j <- matern_corr(Sigma_cross_dist_obs, samp_j["SigmaGP_phi[2]"], Sigma_HP2[2])
      eigen_comp1_j <- samp_j["SigmaGP_mu[1]"]*rep(1,N) + samp_j["SigmaGP_sigma[1]"] * Pmat12_Sigma_obs_j %*% w1_Sigma_j
      eigen_comp2_j <- samp_j["SigmaGP_mu[1]"]*rep(1,N) + samp_j["SigmaGP_sigma[1]"] * Pmat12_Sigma_obs_j %*% w2_Sigma_j
      eigen_comp3_j <- samp_j["SigmaGP_mu[2]"]*rep(1,N) + samp_j["SigmaGP_sigma[2]"] * Pmat3_Sigma_obs_j %*% w3_Sigma_j
      Sigma11_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 1))
      Sigma12_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 3))
      Sigma22_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 2))
      
      # Pred locations
      Pmat12_Sigma_pred_j <- matern_corr(Sigma_cross_dist_pred, samp_j["SigmaGP_phi[1]"], Sigma_HP2[1])
      Pmat3_Sigma_pred_j <- matern_corr(Sigma_cross_dist_pred, samp_j["SigmaGP_phi[2]"], Sigma_HP2[2])
      Peigen_comp1_j <- samp_j["SigmaGP_mu[1]"]*rep(1,M) + samp_j["SigmaGP_sigma[1]"] * Pmat12_Sigma_pred_j %*% w1_Sigma_j
      Peigen_comp2_j <- samp_j["SigmaGP_mu[1]"]*rep(1,M) + samp_j["SigmaGP_sigma[1]"] * Pmat12_Sigma_pred_j %*% w2_Sigma_j
      Peigen_comp3_j <- samp_j["SigmaGP_mu[2]"]*rep(1,M) + samp_j["SigmaGP_sigma[2]"] * Pmat3_Sigma_pred_j %*% w3_Sigma_j
      PSigma11_j <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 1))
      PSigma12_j <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 3))
      PSigma22_j <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 2))
    }
    if( modelsList$Sigma == "npApproxGPIso" ){
      # Required constants: p_Sigma, Sigma_cross_dist_obs, Sigma_cross_dist_pred, Sigma_HP2
      p_Sigma <- constants$p_Sigma
      Sigma_cross_dist_obs <- constants$Sigma_cross_dist_obs
      Sigma_cross_dist_pred <- constants$Sigma_cross_dist_pred
      Sigma_HP2 <- constants$Sigma_HP2
      w1_Sigma_j <- samp_j[paste("w1_Sigma[",1:p_Sigma,"]",sep = "")]
 
      # Obs locations
      Pmat12_Sigma_obs_j <- matern_corr(Sigma_cross_dist_obs, samp_j["SigmaGP_phi[1]"], Sigma_HP2[1])
      eigen_comp1_j <- samp_j["SigmaGP_mu[1]"]*rep(1,N) + samp_j["SigmaGP_sigma[1]"] * Pmat12_Sigma_obs_j %*% w1_Sigma_j
      Sigma11_j <- as.numeric(exp(eigen_comp1_j))
      Sigma12_j <- as.numeric(exp(eigen_comp1_j))
      Sigma22_j <- rep(0,N)
      
      # Pred locations
      Pmat12_Sigma_pred_j <- matern_corr(Sigma_cross_dist_pred, samp_j["SigmaGP_phi[1]"], Sigma_HP2[1])
      Peigen_comp1_j <- samp_j["SigmaGP_mu[1]"]*rep(1,M) + samp_j["SigmaGP_sigma[1]"] * Pmat12_Sigma_pred_j %*% w1_Sigma_j
      PSigma11_j <- as.numeric(exp(Peigen_comp1_j))
      PSigma12_j <- as.numeric(exp(Peigen_comp1_j))
      PSigma22_j <- rep(0,N)
    }
    
    # Calculate mu and Pmu ========================================
    if( modelsList$mu == "constant" ){
      mu <- samp_j[paste("beta")]*rep(1,N)
      Pmu <- samp_j[paste("beta")]*rep(1,M)
    }
    if( modelsList$mu == "linReg" ){
      X_mu <- constants$X_mu
      PX_mu <- constants$PX_mu
      mu <- as.numeric(X_mu %*% samp_j[paste("beta[",1:ncol(X_mu),"]",sep = "")])
      Pmu <- as.numeric(PX_mu %*% samp_j[paste("beta[",1:ncol(X_mu),"]",sep = "")])
    }
    if( modelsList$mu == "zero" ){
      mu <- 0*rep(1,N)
      Pmu <- 0*rep(1,M)
    }
    
    if( modelsList$likelihood == "fullGP" ){ # Predictions for the full GP likelihood
      # Posterior predictive draw ===================================
      # Obs covariance
      Cor <- nsCorr(dist1_sq, dist2_sq, dist12, Sigma11_j, Sigma22_j, Sigma12_j, nu, d)
      sigmaMat <- diag(exp(log_sigma_vec_j))
      Cov <- sigmaMat %*% Cor %*% sigmaMat
      C <- Cov + diag(exp(log_tau_vec_j)^2)
      C_chol <- chol(Cov)
      # Prediction covariance
      PCor <- nsCorr(Pdist1_sq, Pdist2_sq, Pdist12, PSigma11_j, PSigma22_j, PSigma12_j, nu, d)
      PsigmaMat <- diag(exp(Plog_sigma_vec_j))
      PCov <- PsigmaMat %*% PCor %*% PsigmaMat
      if(predict.y){ # Do not include the nugget variance
        PC <- PCov
      } else{ # Include nugget variance
        PC <- PCov + diag(exp(Plog_tau_vec_j)^2)
      }
      # Cross-covariance
      XCor <- nsCrosscorr(Xdist1_sq, Xdist2_sq, Xdist12, 
                          Sigma11_j, Sigma22_j, Sigma12_j,
                          PSigma11_j, PSigma22_j, PSigma12_j, nu, d)
      XCov <- PsigmaMat %*% XCor %*% sigmaMat
      # Conditional mean/covariance
      crscov_covinv <- t(backsolve(C_chol, backsolve(C_chol, t(XCov), transpose = TRUE)))
      condMean <- Pmu + crscov_covinv %*% (z - mu)
      condCov <- PC - crscov_covinv %*% t(XCov)
      condCov_chol <- chol(condCov)
      # Store
      postPredDraws[j,] <- condMean + t(condCov_chol) %*% rnorm(M)
    } else if( modelsList$likelihood == "NNGP" ){ # Predictions for the NNGP likelihood
      # Local kriging ===============================================
      for(m in 1:M){
        # Obs covariance -- for nearest neighbors
        Cor <- nsCorr(Pdist1_3d[m,1:k,1:k],
                      Pdist2_3d[m,1:k,1:k],
                      Pdist12_3d[m,1:k,1:k],
                      Sigma11_j[P_nID[m,]], Sigma22_j[P_nID[m,]], 
                      Sigma12_j[P_nID[m,]], nu, d)
        sigmaMat <- diag(exp(log_sigma_vec_j[P_nID[m,]]))
        Cov <- sigmaMat %*% Cor %*% sigmaMat
        C <- Cov + diag(exp(log_tau_vec_j[P_nID[m,]])^2)
        C_chol <- chol(Cov)
        # Prediction variance
        # PCor <- nsCorr(Pdist1_sq, Pdist2_sq, Pdist12, PSigma11_j, PSigma22_j, PSigma12_j, nu)
        PsigmaMat <- exp(Plog_sigma_vec_j[m])
        # PCov <- PsigmaMat %*% Cor %*% PsigmaMat
        if(predict.y){ # Do not include the nugget variance
          PC <- exp(Plog_sigma_vec_j[m])^2
        } else{ # Include nugget variance
          PC <- exp(Plog_tau_vec_j[m])^2 + exp(Plog_sigma_vec_j[m])^2
        }
        # Cross-covariance
        XCor <- nsCrosscorr(matrix(Pdist1_3d[m,k+1,1:k], nrow = 1, ncol = k),
                            matrix(Pdist2_3d[m,k+1,1:k], nrow = 1, ncol = k),
                            matrix(Pdist12_3d[m,k+1,1:k], nrow = 1, ncol = k),
                            Sigma11_j[P_nID[m,]], 
                            Sigma22_j[P_nID[m,]], 
                            Sigma12_j[P_nID[m,]],
                            PSigma11_j[m], PSigma22_j[m], PSigma12_j[m], nu, d)
        XCov <- PsigmaMat %*% XCor %*% sigmaMat
        # Conditional mean/covariance
        crscov_covinv <- t(backsolve(C_chol, backsolve(C_chol, t(XCov), transpose = TRUE)))
        condMean <- Pmu[m] + crscov_covinv %*% (z[P_nID[m,]] - mu[P_nID[m,]])
        condCov <- PC - crscov_covinv %*% t(XCov)
        condCov_chol <- chol(condCov)
        pred_sample <- condMean + t(condCov_chol) %*% rnorm(1)  
        # Store
        postPredDraws[j,m] <- pred_sample
      }
    } else if( modelsList$likelihood == "SGV" ){ # Predictions for the SGV likelihood
      # SGV prediction ==============================================
      U_j <- calculateU_ns(
        dist1_3d = Alldist1_3d, dist2_3d = Alldist2_3d, dist12_3d = Alldist12_3d,
        Sigma11 = c(Sigma11_j, PSigma11_j), 
        Sigma22 = c(Sigma22_j, PSigma22_j), 
        Sigma12 = c(Sigma12_j, PSigma12_j),
        log_sigma_vec = c(log_sigma_vec_j, Plog_sigma_vec_j), 
        log_tau_vec = c(log_tau_vec_j, Plog_tau_vec_j), nu = nu, 
        nID = prednID_SGV, cond_on_y = predSGV_setup$condition_on_y_ord, 
        N = N, k = k, M = M, d = d )
      Usm <- sparseMatrix(i = U_j[,1], j = U_j[,2], x = U_j[,3])
      Asm <- Usm[c(seq(from = 1, to = 2*N, by = 2), 2*N + 1:M),]
      Bsm <- Usm[seq(from = 2, to = 2*N, by = 2),]
      # Calculate V = rchol(W)
      Aoo_sm <- Asm[1:N,1:(2*N)]
      Woo_sm <- tcrossprod(Aoo_sm)[N:1,N:1]
      Voo_sm <- t(chol(Woo_sm))[N:1,N:1]
      Vsm <- Asm[,-(1:N)]
      Vsm[1:N,1:N] <- Voo_sm
      # Kriging predictor (mean zero)
      ABtz <- crossprod(t(Asm), crossprod(Bsm, z - mu))
      krigPredictor <- -as.numeric(solve(tcrossprod(Vsm), ABtz))
      # Draw mean zero
      pred_sample <- solve(t(Vsm), rnorm(N+M), system = "L")
      # Combine (with mean)
      # postpred_draw <- Pmu + pred_sample[-(1:N)] + krigPredictor[-(1:N)]
      postpred_draw <- c(mu,Pmu) + as.numeric(pred_sample + krigPredictor)
      # Store
      postPredDraws[j,] <- postpred_draw
      PPD_obs_ord <- postPredDraws[,1:N]
      PPD_pred_ord <- postPredDraws[,-(1:N)]
      PPD_obs_orig <- PPD_obs_ord[,order(predSGV_setup$ord)]
      PPD_pred_orig <- PPD_pred_ord[,order(predSGV_setup$ord_pred)]
    } else stop('')
    
    # Progress
    if( j %% ceiling(J/56) == 0 ){cat("-")}
  }
  cat("|\n")
  
  if( modelsList$likelihood == "fullGP" ){ # Predictions for the full GP likelihood
    output <- list(obs = NULL, pred = postPredDraws)
  } else if( modelsList$likelihood == "NNGP" ){ # Predictions for the NNGP likelihood
    output <- list(obs = NULL, pred = postPredDraws)
  } else if( modelsList$likelihood == "SGV" ){ # Predictions for the SGV likelihood
    output <- list(obs = PPD_obs_orig, pred = PPD_pred_orig)
  } else stop('')
  
  return(output)
  
}












