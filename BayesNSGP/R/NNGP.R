

#================================================
# Bayesian nonstationary Gaussian process 
# modeling in NIMBLE
# Mark Risser and Daniel Turek
# Lawrence Berkeley National Laboratory
# January, 2019
#================================================

#================================================
# Functions for the NNGP approximation
#================================================

## Script #3: nsgpNNGP.R (functions for the NNGP approximation)
## 
## - calculateAD_ns
## - calcQF
## - dmnorm_nngp (formerly dmnorm_nn2)
## - rmnorm_nngp (formerly rmnorm_nn2)


#==============================================================================
# Calculate the Gaussian quadratic form for the NNGP approximation
#==============================================================================




# ROxygen comments ----
#' Calculate the Gaussian quadratic form for the NNGP approximation
#' 
#' \code{calcQF} calculates the quadratic form in the multivariate Gaussian 
#' based on the NNGP approximation, for a specific parameter combination. The
#' quadratic form is \code{t(u)C^{-1}v}.
#' 
#' @param u Vector; left product.
#' @param v Vector; right product
#' @param AD N x (k+1) matrix; the first k columns are the 'A' matrix, and the
#' last column is the 'D' vector. Represents the Cholesky of \code{C^{-1}}.
#' @param nID N x k matrix of neighbor indices.
#' 
#' @return A list with two components: (1) an N x 2 array containing the 
#' same spatial coordinates, ordered by MMD, and (2) the same thing, but with 
#' any NA values removed.
#' 
#' @export
#' 
calcQF <- nimbleFunction(
  run = function(u = double(1), v = double(1), AD = double(2), nID = double(2)) {
    N <- dim(AD)[1]
    k <- dim(AD)[2] - 1
    qf <- u[1] * v[1] / AD[1,k+1]
    for(i in 2:N) {
      if(i<=k)     nNei <- i-1      else      nNei <- k
      qf <- qf + (u[i] - inprod( AD[i,1:nNei], u[nID[i,1:nNei]] )) *
        (v[i] - inprod( AD[i,1:nNei], v[nID[i,1:nNei]] )) / AD[i,k+1]
    }
    returnType(double())
    return(qf)
  }
)





#==============================================================================
# Calculate A and D matrices for the NNGP approximation
#==============================================================================

# ROxygen comments ----
#' Calculate A and D matrices for the NNGP approximation
#'
#' \code{calculateAD_ns} calculates A and D matrices (the Cholesky of the 
#' precision matrix) needed for the NNGP approximation.
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
#' @param nID N x k matrix of neighbor indices.
#' @param N Scalar; number of data measurements.
#' @param k Scalar; number of nearest neighbors.
#' @param nu Scalar; Matern smoothness parameter.
#' @param d Scalar; dimension of the spatial domain.
#' 
#' @return A N x (k+1) matrix; the first k columns are the 'A' matrix, and the
#' last column is the 'D' vector.
#'
#' @export
#' 
calculateAD_ns <- nimbleFunction(
  run = function(
    dist1_3d = double(3), dist2_3d = double(3), dist12_3d = double(3),
    Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1),
    log_sigma_vec = double(1), log_tau_vec = double(1), nID = double(2), N = double(), k = double(), nu = double(), d = double() ) {
    AD <- array(0, c(N,k+1))
    AD[1,k+1] <- exp(log_sigma_vec[1])^2 + exp(log_tau_vec[1])^2
    for(i in 2:N) {
      if(i<=k)     nNei <- i-1      else      nNei <- k
      ind <- c( nID[i,1:nNei], i )
      ## these arrays must be extracted, before pssing to nsCorr() function:
      d1 <- dist1_3d[i,1:(nNei+1),1:(nNei+1)]
      d2 <- dist2_3d[i,1:(nNei+1),1:(nNei+1)]
      d12 <- dist12_3d[i,1:(nNei+1),1:(nNei+1)]
      S1 <- Sigma11[ind];      S2 <- Sigma22[ind];      S12 <- Sigma12[ind]
      Cor <- nsCorr(d1, d2, d12, S1, S2, S12, nu, d)
      sigmaMat <- diag(exp(log_sigma_vec[ind]))
      Cov <- sigmaMat %*% Cor %*% sigmaMat
      C <- Cov + diag(exp(log_tau_vec[ind])^2)
      AD[i,1:nNei] <- solve( C[1:nNei,1:nNei], C[nNei+1,1:nNei] )
      AD[i,k+1] <- C[nNei+1,nNei+1] - inprod( C[nNei+1,1:nNei], AD[i,1:nNei] )
    } 
    returnType(double(2))
    return(AD)
  }, check = FALSE
)


#==============================================================================
# Density function for the NNGP approximation
#==============================================================================

# ROxygen comments ----
#' Function for the evaluating the NNGP approximate density.
#'
#' \code{dmnorm_nngp} (and \code{rmnorm_nngp}) calculate the approximate NNGP
#' likelihood for a fixed set of parameters (i.e., A and D matrices). Finally,
#' the distributions must be registered within \code{nimble}.
#' 
#' @param x N-vector of data.
#' @param mean N-vector with current values of the mean
#' @param AD N x (k+1) matrix; the first k columns are the 'A' matrix, and the
#' last column is the 'D' vector.
#' @param nID N x k matrix of neighbor indices.
#' @param N Scalar; number of data measurements.
#' @param k Scalar; number of nearest neighbors.
#' @param log Scalar; should the density be on the log scale (1) or not (0).
#' 
#' @return The NNGP approximate density.
#'
#' @export
#' 
dmnorm_nngp <- nimbleFunction(
  run = function(x = double(1), mean = double(1), AD = double(2), nID = double(2), N = double(), k = double(), log = double()) {
    xCentered <- x - mean
    qf <- calcQF(xCentered, xCentered, AD, nID)
    lp <- -0.5 * (1.83787706649*N + sum(log(AD[1:N,k+1])) + qf)      # log(2pi) = 1.8378770664
    returnType(double())
    return(lp)
  }, check = FALSE
)

# ROxygen comments ----
#' Function for the evaluating the NNGP approximate density.
#'
#' \code{dmnorm_nngp} (and \code{rmnorm_nngp}) calculate the approximate NNGP
#' likelihood for a fixed set of parameters (i.e., A and D matrices). Finally,
#' the distributions must be registered within \code{nimble}.
#' 
#' @param n N-vector of data.
#' @param mean N-vector with current values of the mean
#' @param AD N x (k+1) matrix; the first k columns are the 'A' matrix, and the
#' last column is the 'D' vector.
#' @param nID N x k matrix of neighbor indices.
#' @param N Scalar; number of data measurements.
#' @param k Scalar; number of nearest neighbors.
#' 
#' @return The NNGP approximate density.
#'
#' @export
#' 
rmnorm_nngp <- nimbleFunction(
  run = function(n = integer(), mean = double(1), AD = double(2), nID = double(2), N = double(), k = double()) {
    returnType(double(1))
    return(numeric(N))
  }
)

registerDistributions(list(
  dmnorm_nngp = list(
    BUGSdist = 'dmnorm_nngp(mean, AD, nID, N, k)',
    types = c('value = double(1)', 'mean = double(1)', 'AD = double(2)', 'nID = double(2)', 'N = double()', 'k = double()'),
    mixedSizes = TRUE)
), verbose = FALSE)







