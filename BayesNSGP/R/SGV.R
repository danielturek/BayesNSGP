

#================================================
# Bayesian nonstationary Gaussian process 
# modeling in NIMBLE
# Mark Risser and Daniel Turek
# Lawrence Berkeley National Laboratory
# January, 2019
#================================================

#================================================
# Functions for the SGV approximation
#================================================

## Script #4: nsgpSGV.R (functions for the SGV approximation)
## 
## - conditionLatentObs
## - sgpSetup (formerly sgv_setup)
## - calculateU_ns
## - dmnorm_sgv (formerly sgv_loglikelihood)
## - rmnorm_sgv

#==============================================================================
# Assign conditioning sets for the SGV approximation
#==============================================================================

# ROxygen comments ----
#' Assign conditioning sets for the SGV approximation
#'
#' \code{conditionLatentObs} assigns q_y(i) vs q_z(i) following Section 5.1 
#' in Katzfuss and Guinness (2018). This function only needs to be run once 
#' per SGV analysis.
#' 
#' @param nID N x k matrix of neighbor indices.
#' @param locs_ord N x 2 matrix of locations.
#' @param N Scalar; number of locations (observed only!).
#' 
#' @return A matrix indicating whether the conditioning set for each location 
#' is on the latent process (y, \code{1}) or the observed values (z, \code{0}).
#'
#' @examples
#' # TODO
#'
#' @export
#' 
conditionLatentObs <- function( nID, locs_ord, N ){
  
  Nall <- nrow(nID)
  k <- ncol(nID)
  d <- ncol(locs_ord)
  
  cond_on_y <- matrix(0, Nall, k) 
  cond_on_y[nID == -1] <- -1 ## populate unused values with -1, to prevent a warning from NIMBLE
  cond_on_y[2,1] <- 1
  for(i in 3:N){ # i = 1 has no conditioning set; i = 2 automatically conditions on y_1
    q_i <- (nID[i,])[nID[i,] != -1] 
    size_intrsct_qyj_qi <- rep(NA, length(q_i))
    for(j in 1:length(q_i)){
      q_y_j <- which(cond_on_y[q_i[j],] == 1)
      size_intrsct_qyj_qi[j] <- sum(q_y_j %in% q_i)
    }
    ind_h_i <- which(size_intrsct_qyj_qi == max(size_intrsct_qyj_qi))
    h_i <- q_i[ind_h_i]
    ind_k_i <- which.min( as.numeric(mahalanobis.dist(data.x = matrix(locs_ord[i,], ncol = d, byrow = TRUE),
                                                                 data.y = matrix(locs_ord[h_i,], ncol = d, byrow = TRUE),
                                                                 vc = diag(d))) ) 
    k_i <- h_i[ind_k_i]
    
    q_y_k_i <- nID[k_i,which(cond_on_y[k_i,] == 1)]
    cond_on_y[i, which(c(q_y_k_i[which( q_y_k_i %in% q_i )], k_i) %in% (nID[i,])[nID[i,] != -1])] <- 1
  }
  if(Nall > N) cond_on_y[N:Nall,] <- 1 # This case involves prediction locations
  
  return(cond_on_y)
}

#==============================================================================
# One-time setup wrapper function for the SGV approximation
#==============================================================================

# ROxygen comments ----
#' One-time setup wrapper function for the SGV approximation
#'
#' \code{sgvSetup} is a wrapper function that sets up the SGV approximation. 
#' Three objects are required: (1) ordering the locations, (2) identify nearest
#' neighbors, and (3) determine the conditioning set. This function only needs
#' to be run once per SGV analysis.
#' 
#' @param locs Matrix of observed locations.
#' @param locs_pred Optional matrix of prediction locations.
#' @param k Number of neighbors.
#' @param seed TODO
#' 
#' @return A list with the following components:
#' \item{ord}{A vector of ordering position for the observed locations.}
#' \item{ord_pred}{A vector of ordering position for the prediction 
#' locations (if \code{locs_pred} is provided).}
#' \item{ord_all}{A concatenated vector of \code{ord} and \code{ord_pred}.}
#' \item{locs_ord}{A matrix of ordered locations (observed and prediction),
#' included for convenience.}
#' \item{nID_ord}{A matrix of (ordered) neighbor indices.}
#' \item{condition_on_y_ord}{A matrix indicating whether the conditioning
#' set for each (ordered) location is on the latent process (y, \code{1}) or
#' the observed values (z, \code{0}).}
#'
#' @examples
#' # TODO
#'
#' @export
#' 
sgvSetup <- function( locs, locs_pred = NULL, k = 15, seed = NULL ){
  
  if(is.null(seed)) seed <- sample(1e5, 1) # Set seed for reproducibility (randomness in orderCoordinatesMMD function)
  
  d <- ncol(locs) # Spatial dimension
  n <- nrow(locs) # Number of (observed) locations
  
  #--------------------------------------------------------
  # Task 1: Order the locations
  #--------------------------------------------------------
  set.seed(seed)
  if(is.null(locs_pred)){ # If no prediction
    locs_mmd <- orderCoordinatesMMD(locs )
    locs_pred_mmd <- NULL
    ord <- locs_mmd$orderedIndicesNoNA
    ord_pred <- NULL
    ord_all <- ord
    locs_ord <- locs[ord_all,]
    
  } else{
    n_pred <- nrow(locs_pred) # Number of prediction locations
    # locs_ord <- orderCoordinatesMMD(locs)
    # locs_pred_ord <- orderCoordinatesMMD(locs_pred)
    
    locs_mmd <- orderCoordinatesMMD(locs)
    locs_pred_mmd <- orderCoordinatesMMD(locs_pred)
    ord <- locs_mmd$orderedIndicesNoNA
    ord_pred <- locs_pred_mmd$orderedIndicesNoNA
    ord_all <- c(ord, n+ord_pred)
    locs_ord <- rbind(locs, locs_pred)[ord_all,]
  }
  
  #--------------------------------------------------------
  # Task 2: Get nearest neighbors
  #--------------------------------------------------------
  nID_ord <- determineNeighbors(locs_ord, k)
  
  #--------------------------------------------------------
  # Task 3: Conditioning on y or z?
  #--------------------------------------------------------
  condition_on_y_ord <- conditionLatentObs( nID_ord, locs_ord, n )
  
  return(list( seed = seed, 
               ord = ord, ord_pred = ord_pred, ord_all = ord_all, 
               locs_ord = locs_ord, nID_ord = nID_ord, 
               condition_on_y_ord = condition_on_y_ord ))
}

#==============================================================================
# Calculate the (sparse) matrix U 
#==============================================================================

# ROxygen comments ----
#' Calculate the (sparse) matrix U 
#'
#' \code{calculateU_ns} calculates the (sparse) matrix U (i.e., the Cholesky 
#' of the inverse covariance matrix) using a nonstationary covariance function.
#' The output only contains non-zero values and is stored as three vectors: 
#' (1) the row indices, (2) the column indices, and (3) the non-zero values.
#' NOTE: this code assumes the all inputs correspond to the ORDERED locations.
#' 
#' @param dist1_3d N x (k+1) x (k+1) array of distances in the x-coordinate 
#' direction.
#' @param dist2_3d N x (k+1) x (k+1) array of distances in the y-coordinate 
#' direction.
#' @param dist12_3d N x (k+1) x (k+1) array of cross-distances.
#' @param Sigma11 N-vector; 1-1 element of the Sigma() process.
#' @param Sigma12 N-vector; 1-2 element of the Sigma() process.
#' @param Sigma22 N-vector; 2-2 element of the Sigma() process.
#' @param log_sigma_vec N-vector; process standard deviation values.
#' @param log_tau_vec N-vector; nugget standard deviation values.
#' @param nu Scalar; Matern smoothness parameter.
#' @param nID N x k matrix of (ordered) neighbor indices.
#' @param cond_on_y A matrix indicating whether the conditioning set for each 
#' (ordered) location is on the latent process (y, \code{1}) or the observed 
#' values (z, \code{0}). Calculated in \code{sgvSetup}.
#' @param N Scalar; number of data measurements.
#' @param k Scalar; number of nearest neighbors.
#' @param d TODO
#' @param M TODO
#' 
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @export
#' 
calculateU_ns <- nimbleFunction(  # Create the sparse U matrix for specific theta
  run = function(
    dist1_3d = double(3), dist2_3d = double(3), dist12_3d = double(3),
    Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1),
    log_sigma_vec = double(1), log_tau_vec = double(1), nu = double(), 
    nID = double(2), cond_on_y = double(2), N = double(), k = double(), d = double(0), 
    M = double(0, default = 0) ) {
    
    # Setup
    NN <- 2*N + M
    num_NZ <- 3*N + k*N - (k*(k+1)/2) + (k+1)*M # Number of non-zero entries in U
    num_neigbs <- c(0, seq(from = 1, to = k-1, by = 1), array(k, M+N-k))
    Uvals <- array(0, num_NZ)
    rowInd <- array(0, num_NZ)
    colInd <- array(0, num_NZ)
    
    # Calculate the position of the diagonal elements    
    dgIdx_vec <- array(-1, N+M)
    for(l in 1:k){
      dgIdx_vec[l] <- 1 + sum(num_neigbs[1:l])
    }
    dgIdx_vec[(k+1):(N+M)] <- seq(from=(k*(k+1)/2)+1, to=num_NZ - 2*N, by = k+1)
    
    # First: the y_j
    Uvals[1] <- exp(log_sigma_vec[1])^2
    rowInd[1] <- 1
    colInd[1] <- 1
    for(i in 2:(N+M)){ # y_j
  
      if(i<=k)     nNei <- i-1      else      nNei <- k
      ind <- nID[i,1:nNei]
      ## these arrays must be extracted, before pssing to nsCorr() and nsCrosscorr() function:
      Xd1 <- array(dist1_3d[i, 1:nNei, (nNei + 1)], c(nNei, 1))     # Distances between location i and its neighbors
      Xd2 <- array(dist2_3d[i, 1:nNei, (nNei + 1)], c(nNei, 1))
      Xd12 <- array(dist12_3d[i, 1:nNei, (nNei + 1)], c(nNei, 1))
      S1 <- nimNumeric( value = Sigma11[i], length = 1)                           # Anisotropy parameters for location i
      S2 <- nimNumeric( value = Sigma22[i], length = 1)  
      S12 <- nimNumeric( value = Sigma12[i], length = 1)  
      
      d1 <- array(dist1_3d[i, 1:nNei, 1:nNei], c(nNei, nNei))          # Distances between the neighbors of location i
      d2 <- array(dist2_3d[i, 1:nNei, 1:nNei], c(nNei, nNei))
      d12 <- array(dist12_3d[i, 1:nNei, 1:nNei], c(nNei, nNei))
      xS1 <- Sigma11[ind]                        # Anisotropy parameters for the neighbors of location i
      xS2 <- Sigma22[ind]
      xS12 <- Sigma12[ind]
      
      # Cross-covariance between location and the conditioning set
      Crosscor <- nsCrosscorr(Xd1, Xd2, Xd12, S1, S2, S12, xS1, xS2, xS12, nu, d)
      if(length(ind) == 1) {
        sigmaMat_cond <- array(exp(log_sigma_vec[ind]), c(1,1))
      } else {
        sigmaMat_cond <- diag(exp(log_sigma_vec[ind]))
      }
      Crosscov <- array(exp(log_sigma_vec[i]), c(1,1)) %*% array(Crosscor, c(1,nNei)) %*% sigmaMat_cond # Formerly pt1
      
      # Covariance of conditioning set
      Cor_cond <- nsCorr(d1, d2, d12, xS1, xS2, xS12, nu, d)
      Cov_cond <- sigmaMat_cond %*% Cor_cond %*% sigmaMat_cond # Formerly pt2
      
      # Covariance of the process at the location
      Cov_loc <- exp(log_sigma_vec[i])^2 # Formerly pt3
      
      ####################################
      # b_i <- nimNumeric( value = (Crosscov %*% inverse(Cov_cond))[1,1:nNei], length = nNei)
      b_i <- nimNumeric( value = solve(Cov_cond, t(Crosscov))[1:nNei,1], length = nNei)
      r_i <- Cov_loc - inprod( b_i, nimNumeric(value = (Crosscov)[1,1:nNei], length = nNei) ) 
      ####################################

      # Store
      dgIdx <- dgIdx_vec[i]
      Uvals[dgIdx] <- 1/sqrt(r_i)
      if(i > N){
        rowInd[dgIdx] <- N + i 
        colInd[dgIdx] <- N + i 
      } else{
        rowInd[dgIdx] <- 2*(i-1)+1
        colInd[dgIdx] <- 2*(i-1)+1
      }
      
      for(j in 1:nNei){
        if(cond_on_y[i,j] == 1){ # condition on y_i
          Uvals[dgIdx + j] <- -b_i[j]/sqrt(r_i)
          if(i > N){ # Pred locations
            if(ind[j] > N){
              rowInd[dgIdx + j] <- N + ind[j]
            } else{
              rowInd[dgIdx + j] <- 2*(ind[j]-1)+1
            }
            colInd[dgIdx + j] <- N + i 
          } else{ # Obs locations
            rowInd[dgIdx + j] <- 2*(ind[j]-1)+1
            colInd[dgIdx + j] <- 2*(i-1)+1
          }
        } else{ # condition on z_i
          Uvals[dgIdx + j] <- -b_i[j]/sqrt(r_i)
          if(i > N){ # Pred locations
            stop("Error.")
            colInd[dgIdx + j] <- N + i
          } else{ # Obs locations
            rowInd[dgIdx + j] <- 2*ind[j]
            colInd[dgIdx + j] <- 2*(i-1)+1
          }
        }
      }
    } 
    
    # Next: the z_i
    Uvals[(num_NZ - (2*N) + 1):(num_NZ - N)] <- 1/exp(log_tau_vec[1:N])
    rowInd[(num_NZ - (2*N) + 1):(num_NZ - N)] <- seq(from = 2, to = (2*N), by = 2)
    colInd[(num_NZ - (2*N) + 1):(num_NZ - N)] <- seq(from = 2, to = (2*N), by = 2)
    Uvals[(num_NZ - N + 1):num_NZ] <- -1/exp(log_tau_vec[1:N])
    rowInd[(num_NZ - N + 1):num_NZ] <- seq(from = 1, to = (2*N), by = 2)
    colInd[(num_NZ - N + 1):num_NZ] <- seq(from = 2, to = (2*N), by = 2)
    
    # Combine
    U_ijx <- array(0, c(num_NZ, 3))
    U_ijx[,1] <- rowInd
    U_ijx[,2] <- colInd
    U_ijx[,3] <- Uvals
    
    returnType(double(2))
    return(U_ijx)
  }  ##, check = FALSE
)


R_sparse_tcrossprod <- function(i, j, x, subset = -1) {
  Asparse <- sparseMatrix(i = i, j = j, x = x)
  if(subset[1] < 0){ # No subset
    ans.dsCMatrix <- tcrossprod(Asparse)
  } else{
    ans.dsCMatrix <- tcrossprod(Asparse[subset,])
  }
  ans.dgTMatrix <- as(ans.dsCMatrix, 'dgTMatrix')
  i <- ans.dgTMatrix@i + 1
  j <- ans.dgTMatrix@j + 1
  x <- ans.dgTMatrix@x
  ijx <- cbind(i, j, x)
  return(ijx)
}

# ROxygen comments ----
#' nimble_sparse_tcrossprod
#' @param i TODO
#' @param j TODO
#' @param x TODO
#' @param subset TODO
#' @export
nimble_sparse_tcrossprod <- nimbleRcall(
  prototype = function(i = double(1), j = double(1), x = double(1), subset = double(1)) {},
  returnType = double(2),
  Rfun = 'R_sparse_tcrossprod'
)

R_sparse_crossprod <- function(i, j, x, z, n, subset = -1, transp = 1) {
  ## Mark: TODO.  What should be on the next line??  
  zSparse <- array(1:9, c(3,3))  ## sparseMatrix(i = 1:n, j = rep(1,n), x = as.numeric(z))
  if(transp == 1){ # use crossprod
    Asparse <- sparseMatrix(i = i, j = j, x = x)
    if(subset[1] < 0){ # No subset
      ans.dsCMatrix <- crossprod(Asparse, zSparse)
    } else{
      ans.dsCMatrix <- crossprod(Asparse[subset,], as.numeric(z))
    }
  } else{ # Use %*%
    Asparse <- sparseMatrix(i = j, j = i, x = x)
    if(subset[1] < 0){ # No subset
      ans.dsCMatrix <- crossprod(Asparse, zSparse)
    } else{
      ans.dsCMatrix <- crossprod(Asparse[,subset], as.numeric(z))
    }
  }
  return(ans.dsCMatrix@x)
}

# ROxygen comments ----
#' nimble_sparse_crossprod
#' @param i TODO
#' @param j TODO
#' @param x TODO
#' @param z TODO
#' @param n TODO
#' @param subset TODO
#' @param transp TODO
#' @export
nimble_sparse_crossprod <- nimbleRcall(
  prototype = function(i = double(1), j = double(1), x = double(1), z = double(1), n = double(), subset = double(1), transp = double()) {},
  returnType = double(1),
  Rfun = 'R_sparse_crossprod'
)

# ROxygen comments ----
#' R_sparse_chol
#' @param i TODO
#' @param j TODO
#' @param x TODO
#' @param n TODO
#' @export
R_sparse_chol <- function(i, j, x, n) {
  Asparse <- sparseMatrix(i = i, j = j, x = x)
  ans.dsCMatrix <- t(chol(Asparse[n:1,n:1]))
  ans.dgTMatrix <- as(ans.dsCMatrix, 'dgTMatrix')
  i <- ans.dgTMatrix@i + 1
  j <- ans.dgTMatrix@j + 1
  x <- ans.dgTMatrix@x
  ijx <- cbind(i, j, x)
  return(ijx)
}

# ROxygen comments ----
#' nimble_sparse_chol
#' @param i TODO
#' @param j TODO
#' @param x TODO
#' @param n TODO
#' @export
nimble_sparse_chol <- nimbleRcall(
  prototype = function(i = double(1), j = double(1), x = double(1), n = double()) {},
  returnType = double(2),
  Rfun = 'R_sparse_chol'
)

R_sparse_solve <- function(i, j, x, z) {
  # z3 <- solve(V_ord, rev(z2), system = "L")
  Asparse <- sparseMatrix(i = i, j = j, x = x)
  z_rev <- rev(z)
  ans.dsCMatrix <- solve(Asparse, z_rev, system = "L")
  return(ans.dsCMatrix@x)
}

# ROxygen comments ----
#' nimble_sparse_solve
#' @param i TODO
#' @param j TODO
#' @param x TODO
#' @param z TODO
#' @export
nimble_sparse_solve <- nimbleRcall(
  prototype = function(i = double(1), j = double(1), x = double(1), z = double(1)) {},
  returnType = double(1),
  Rfun = 'R_sparse_solve'
)


#==============================================================================
# Density function for the SGV approximation
#==============================================================================

# ROxygen comments ----
#' Function for the evaluating the SGV approximate density.
#'
#' \code{dmnorm_sgv} (and \code{rmnorm_sgv}) calculate the approximate SGV
#' likelihood for a fixed set of parameters (i.e., the U matrix). Finally,
#' the distributions must be registered within \code{nimble}.
#' 
#' @param x TODO
#' @param mean TODO
#' @param U TODO
#' @param N TODO
#' @param k TODO
#' @param log TODO
#' 
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @export
#' 
dmnorm_sgv <- nimbleFunction(
  run = function(x = double(1), mean = double(1), U = double(2),
                 N = double(), k = double(), log = double(0, default = 1)) {
    # Components
    zo_ord <- x
    z1 <- nimble_sparse_crossprod(
        i = U[,1], j = U[,2], x = U[,3], z = zo_ord - mean, n = N,
        subset = seq(from = 2, to = 2*N, by = 2), transp = 1)
    logdet_U <- -sum(log(U[U[,1] == U[,2],3]))
    z2 <- nimble_sparse_crossprod(
        i = U[,1], j = U[,2], x = U[,3], z = z1, n = N,
        subset = seq(from = 1, to = 2*N, by = 2), transp = 0)
    Amat <- nimble_sparse_tcrossprod(
        i = U[,1], j = U[,2], x = U[,3], 
        subset = seq(from = 1, to = 2*N, by = 2))
    Vmat_ord <- nimble_sparse_chol(i = Amat[,1], j = Amat[,2], x = Amat[,3], n = N)
    logdet_V <- sum(log(Vmat_ord[Vmat_ord[,1] == Vmat_ord[,2],3]))
    z3 <- nimble_sparse_solve(i = Vmat_ord[,1], j = Vmat_ord[,2], x = Vmat_ord[,3], z = z2)
    lp <- -(logdet_U + logdet_V + 0.5*sum(z1^2) - 0.5*sum(z3^2)) - 0.5*1.83787706649*N
    returnType(double())
    return(lp)
  }  ##, check = FALSE
)

rmnorm_sgv <- nimbleFunction(
  run = function(n = integer(), mean = double(1), U = double(2), N = double(), k = double()) {
    returnType(double(1))
    return(numeric(N))
  }
)

registerDistributions(list(
  dmnorm_sgv = list(
    BUGSdist = 'dmnorm_sgv(mean, U, N, k)',
    types = c('value = double(1)', 'mean = double(1)', 'U = double(2)', 'N = double()', 'k = double()'),
    mixedSizes = TRUE)
), verbose = FALSE)










