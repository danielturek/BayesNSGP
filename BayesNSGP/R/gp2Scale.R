

#================================================
# Bayesian nonstationary Gaussian process 
# modeling in NIMBLE
# Mark Risser and Daniel Turek
# Lawrence Berkeley National Laboratory
# January, 2019
#================================================

#================================================
# Functions for the gp2Scale likelihood option
#================================================

## Functions needed
## 
## - sparseKernel (calculate the sparse kernel and determine nonzero entries)
## - calculateSparseC (calculate sparse covariance)
## - 
## - dmnorm_gp2Scale (dmnorm with sparse covariance) 
## - rmnorm_gp2Scale (rmnorm with sparse covariance)

#==============================================================================
# gp2Scale: calculate the sparse kernel and determine nonzero entries
#==============================================================================

# ROxygen comments ----
#' Calculate sparse kernel, core kernel, and determine nonzero entries
#'
#' \code{Cy_sm} calculates the normalized sparse kernel for a fixed
#' set of bump function hyperparameters and returns the nonzero entries. Note
#' that the matrix is calculated and returned in dense format.
#' 
#' @param dists N x N matrix of Euclidean distances
#' @param coords N x d matrix of coordinate/input locations
#' @param N Scalar; number of data measurements.
#' @param d Scalar; dimension of the spatial domain.
#' @param n1 Scalar; number of outer products.
#' @param n2 Scalar; number of bump functions in each outer product.
#' @param r0 Scalar; length-scale of sparse stationary kernel.
#' @param s0 Scalar; signal-variance of sparse stationary kernel.
#' @param cstat_opt Scalar; determines the compactly supported kernel. See Details.
#' @param normalize Logical; should C_sparse have 1's along the diagonal
#' @param bumpLocs Array of bump function locations (n2*d x n1)
#' @param rads Matrix of bump function radii (n1 x n2; denoted \eqn{r_{ij}})
#' @param ampls Matrix of bump function amplitudes (n1 x n2; denoted \eqn{a_{ij}})
#' @param shps Matrix of bump function shape parameters (n1 x n2; denoted \eqn{b_{ij}})
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
#' @param log_sigma_vec Vector of length N; log of the signal standard deviation.
#' @param lognuggetSD Vector of length N; log of the error standard deviation.
#' 
#' 
#' @return Returns a sparse matrix (N x 3) of the nonzero elements of the product 
#' between the core and sparse kernel.
#'
#' @export
#' 
Cy_sm <- nimbleFunction(  # Generate the sparse kernel in dense format
  run = function(
    dists = double(2), 
    coords = double(2), # coordinate / input locations
    N = double(0), # number of data points
    d = double(0), n1 = double(0), n2 = double(0), # dimension of input space; number of outer prods; number of bump fcns
    r0 = double(0), # scalar length-scale for compact stationary
    s0 = double(0), # Scalar; signal-variance of sparse stationary kernel.
    cstat_opt = double(0), # indicator for which compactly supported kernel is used
    normalize = double(0), # Logical; should C_sparse have 1's along the diagonal (1 = TRUE)
    bumpLocs = double(2), # array of bump function locations (n2*d x n1)
    rads = double(2), ampls = double(2), shps = double(2), # matrix of radii, amplitudes, shapes (n1 x n2)
    dist1_sq = double(2), dist2_sq = double(2), dist12 = double(2), 
    Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1), 
    nu = double(0), log_sigma_vec = double(1), lognuggetSD = double(1)
  ) {

    if(n1 == 1 | n2 == 1){
      bumpLocs <- matrix(bumpLocs, nrow=d*n2, ncol=n1)
      rads <- matrix(rads, nrow=n1, ncol=n2)
      ampls <- matrix(ampls, nrow=n1, ncol=n2)
      shps <- matrix(shps, nrow=n1, ncol=n2)
    }
    
    #=============== Sparse Kernel =================#
    # Calculate the bump function sums + cross products
    for(j in 1:n1){
      g_j <- rep(0,N)
      for(b in 1:n2){
        # Calculate coords x bump fcn distances
        tmp <- rep(0,N)
        for(l in 1:d){
          rwix <- (n2*(l-1) + b)
          tmp <- tmp + (coords[1:N,l] - rep(bumpLocs[rwix,j],N))^2
        }
        dist_jb <- sqrt(tmp)
        ind_nz <- which(dist_jb <= rads[j,b])
        g_j[ind_nz] <- g_j[ind_nz] + ampls[j,b] * exp(shps[j,b] - shps[j,b]/(1 - dist_jb[ind_nz]^2/rads[j,b]^2))
        # if(ampls[j,b] < 0){
        #   a_ij <- 0
        # } else{
        #   a_ij <- 1
        # }
        # g_j[ind_nz] <- g_j[ind_nz] + a_ij * exp(shps[j,b] - shps[j,b]/(1 - dist_jb[ind_nz]^2/rads[j,b]^2))
      }
      if(j == 1){
        crsPrd <- g_j %*% asRow(g_j)
      } else{
        crsPrd <- crsPrd + g_j %*% asRow(g_j)
      }
    }
    
    # Add on the sparse stationary; normalize and store as sparse matrix
    k_c <- matrix(0,N,N)
    for(i in 1:N){
      k_c[i,i] <- s0
      if(i < N){
        tmpDst <- dists[i,(i+1):N]
        ind_dist <- which(tmpDst <= r0)
        # Choices for compactly supported kernel
        if(cstat_opt == 1){ # TPF(a=1, v=a+0.5+(d-1)/2)
          a <- 1; v <- a + 0.5 + (d-1)/2
          tmpK <- (1-(tmpDst[ind_dist]/r0)^a)^v
        }
        if(cstat_opt == 2){ # Wendland degree 2
          tmpK <- (1-abs(tmpDst[ind_dist]/r0))^4*(4*abs(tmpDst[ind_dist]/r0) + 1)
        }
        if(cstat_opt == 3){ # Wendland degree 3
          tmpK <- (1/3)*(1-abs(tmpDst[ind_dist]/r0))^6*(35*abs(tmpDst[ind_dist]/r0)^2 + 18*abs(tmpDst[ind_dist]/r0) + 3)
        }
        if(cstat_opt == 4){ # Wendland degree 4
          tmpK <- (1-abs(tmpDst[ind_dist]/r0))^8*(32*abs(tmpDst[ind_dist]/r0)^3 + 25*abs(tmpDst[ind_dist]/r0)^2 + 8*abs(tmpDst[ind_dist]/r0) + 1)
        }
        if(cstat_opt == 5){ # M and R (2009)
          tmpK <- ((3.0*(tmpDst[ind_dist]/r0)^2*log((tmpDst[ind_dist]/r0)/(1+sqrt(1.0 - (tmpDst[ind_dist]/r0)^2)))) +
                     ((2.0*(tmpDst[ind_dist]/r0)^2+1.0)*sqrt(1.0-(tmpDst[ind_dist]/r0)^2))) #* (sqrt(2.0)/(3.0*sqrt(pi)))
        }
        k_c[i,ind_dist+i] <- s0*tmpK
        k_c[ind_dist+i,i] <- s0*tmpK
        
      }
    }
    sparseKernel_unnorm <- crsPrd + k_c
    sparse_diag_sqrt <- sqrt(diag(sparseKernel_unnorm))
    sparseKernel <- matrix(0,N,N)
    if(normalize == 1){
      diag(sparseKernel) <- 1
    } else{
      diag(sparseKernel) <- diag(sparseKernel_unnorm)
    }
    for(i in 1:(N-1)){
      # Sparse
      Ci <- sparseKernel_unnorm[i,(i+1):N]
      ind_nz_i <- which(Ci != 0)
      # n_nz_i <- length(ind_nz_i)
      if(normalize == 1){
        Ci_norm_nz <- Ci[ind_nz_i]/(sparse_diag_sqrt[ind_nz_i+i]*sparse_diag_sqrt[i])
      } else{
        Ci_norm_nz <- Ci[ind_nz_i]
      }
      sparseKernel[i,ind_nz_i+i] <- Ci_norm_nz
      sparseKernel[ind_nz_i+i,i] <- Ci_norm_nz
    }
    
    # Sparsify
    Cvals <- rep(0,N^2)
    rowInd <- rep(-1,N^2)
    colInd <- rep(-2,N^2)
    Cvals[1:N] <- diag(sparseKernel)
    rowInd[1:N] <- 1:N
    colInd[1:N] <- 1:N
    ctr <- N
    for(i in 1:(N-1)){
      Ci <- sparseKernel[i,(i+1):N]
      ind_nz_i <- which(Ci != 0)
      n_nz_i <- length(ind_nz_i)
      if(n_nz_i > 0){
        Ci_nz <- Ci[ind_nz_i]
        Cvals[(ctr+1):(ctr+2*n_nz_i)] <- c(Ci_nz,Ci_nz)
        rowInd[(ctr+1):(ctr+2*n_nz_i)] <- c(rep(i,n_nz_i),ind_nz_i+i)
        colInd[(ctr+1):(ctr+2*n_nz_i)] <- c(ind_nz_i+i,rep(i,n_nz_i))
        ctr <- ctr + 2*n_nz_i
      }
    }
    
    #=============== Core Kernel =================#
    Cvals[rowInd == colInd] <- (Cvals[rowInd == colInd] * (exp(log_sigma_vec)^2)) + exp(lognuggetSD)^2
    for(i in 1:(N-1)){
      # Ci <- Cvals[]
      ind_nz_i <- which(rowInd == i & colInd > i)
      ind_nz_i_v2 <- which(colInd == i & rowInd > i)
      if(length(ind_nz_i) > 0){
        colInd_rowI <- colInd[ind_nz_i]
        tmp11 <- Sigma11[colInd_rowI]
        tmp22 <- Sigma22[colInd_rowI]
        tmp12 <- Sigma12[colInd_rowI]
        cor_i <- nsCrosscorr( Xdist1_sq = matrix(dist1_sq[i,colInd_rowI], ncol = 1),
                              Xdist2_sq = matrix(dist2_sq[i,colInd_rowI], ncol = 1),
                              Xdist12 = matrix(dist12[i,colInd_rowI], ncol = 1),
                              Sigma11 = nimNumeric(Sigma11[i], length = 1), 
                              Sigma22 = nimNumeric(Sigma22[i], length = 1), 
                              Sigma12 = nimNumeric(Sigma12[i], length = 1),
                              PSigma11 = tmp11, PSigma22 = tmp22, PSigma12 = tmp12, 
                              nu = nu, d = d ) 
        cov_i <- nimNumeric(cor_i, length = length(colInd_rowI))*exp(log_sigma_vec[colInd_rowI])
        Cvals[ind_nz_i] <- Cvals[ind_nz_i] * cov_i * exp(log_sigma_vec[i])
        Cvals[ind_nz_i_v2] <- Cvals[ind_nz_i_v2] * cov_i * exp(log_sigma_vec[i])
      }
    }
    Cy_sm <- matrix(0,N^2,3) # cbind(rowInd, colInd, Cvals)
    Cy_sm[,1] <- rowInd
    Cy_sm[,2] <- colInd
    Cy_sm[,3] <- Cvals
    
    returnType(double(2))
    return(Cy_sm)
  }, check = FALSE, name = "Cy_sm"
)

Cy_sm_n11 <- nimbleFunction(  # Generate the sparse kernel in dense format
  run = function(
    dists = double(2), 
    coords = double(2), # coordinate / input locations
    N = double(0), # number of data points
    d = double(0), n1 = double(0), n2 = double(0), # dimension of input space; number of outer prods; number of bump fcns
    r0 = double(0), # scalar length-scale for compact stationary
    s0 = double(0), # Scalar; signal-variance of sparse stationary kernel.
    cstat_opt = double(0), # indicator for which compactly supported kernel is used
    normalize = double(0), # Logical; should C_sparse have 1's along the diagonal (1 = TRUE)
    bumpLocs = double(1), # array of bump function locations (n2*d x n1)
    rads = double(1), ampls = double(1), shps = double(1), # matrix of radii, amplitudes, shapes (n1 x n2)
    dist1_sq = double(2), dist2_sq = double(2), dist12 = double(2), 
    Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1), 
    nu = double(0), log_sigma_vec = double(1), lognuggetSD = double(1)
  ) {
    
    #=============== Sparse Kernel =================#
    # Calculate the bump function sums + cross products
    for(j in 1:n1){
      g_j <- rep(0,N)
      for(b in 1:n2){
        # Calculate coords x bump fcn distances
        tmp <- rep(0,N)
        for(l in 1:d){
          rwix <- (n2*(l-1) + b)
          tmp <- tmp + (coords[1:N,l] - rep(bumpLocs[rwix],N))^2
        }
        dist_jb <- sqrt(tmp)
        ind_nz <- which(dist_jb <= rads[b])
        g_j[ind_nz] <- g_j[ind_nz] + ampls[b] * exp(shps[b] - shps[b]/(1 - dist_jb[ind_nz]^2/rads[b]^2))
        # if(ampls[j,b] < 0){
        #   a_ij <- 0
        # } else{
        #   a_ij <- 1
        # }
        # g_j[ind_nz] <- g_j[ind_nz] + a_ij * exp(shps[j,b] - shps[j,b]/(1 - dist_jb[ind_nz]^2/rads[j,b]^2))
      }
      if(j == 1){
        crsPrd <- g_j %*% asRow(g_j)
      } else{
        crsPrd <- crsPrd + g_j %*% asRow(g_j)
      }
    }
    
    # Add on the sparse stationary; normalize and store as sparse matrix
    k_c <- matrix(0,N,N)
    for(i in 1:N){
      k_c[i,i] <- s0
      if(i < N){
        tmpDst <- dists[i,(i+1):N]
        ind_dist <- which(tmpDst <= r0)
        # Choices for compactly supported kernel
        if(cstat_opt == 1){ # TPF(a=1, v=a+0.5+(d-1)/2)
          a <- 1; v <- a + 0.5 + (d-1)/2
          tmpK <- (1-(tmpDst[ind_dist]/r0)^a)^v
        }
        if(cstat_opt == 2){ # Wendland degree 2
          tmpK <- (1-abs(tmpDst[ind_dist]/r0))^4*(4*abs(tmpDst[ind_dist]/r0) + 1)
        }
        if(cstat_opt == 3){ # Wendland degree 3
          tmpK <- (1/3)*(1-abs(tmpDst[ind_dist]/r0))^6*(35*abs(tmpDst[ind_dist]/r0)^2 + 18*abs(tmpDst[ind_dist]/r0) + 3)
        }
        if(cstat_opt == 4){ # Wendland degree 4
          tmpK <- (1-abs(tmpDst[ind_dist]/r0))^8*(32*abs(tmpDst[ind_dist]/r0)^3 + 25*abs(tmpDst[ind_dist]/r0)^2 + 8*abs(tmpDst[ind_dist]/r0) + 1)
        }
        if(cstat_opt == 5){ # M and R (2009)
          tmpK <- ((3.0*(tmpDst[ind_dist]/r0)^2*log((tmpDst[ind_dist]/r0)/(1+sqrt(1.0 - (tmpDst[ind_dist]/r0)^2)))) +
                     ((2.0*(tmpDst[ind_dist]/r0)^2+1.0)*sqrt(1.0-(tmpDst[ind_dist]/r0)^2))) #* (sqrt(2.0)/(3.0*sqrt(pi)))
        }
        k_c[i,ind_dist+i] <- s0*tmpK
        k_c[ind_dist+i,i] <- s0*tmpK
        
      }
    }
    sparseKernel_unnorm <- crsPrd + k_c
    sparse_diag_sqrt <- sqrt(diag(sparseKernel_unnorm))
    sparseKernel <- matrix(0,N,N)
    if(normalize == 1){
      diag(sparseKernel) <- 1
    } else{
      diag(sparseKernel) <- diag(sparseKernel_unnorm)
    }
    for(i in 1:(N-1)){
      # Sparse
      Ci <- sparseKernel_unnorm[i,(i+1):N]
      ind_nz_i <- which(Ci != 0)
      # n_nz_i <- length(ind_nz_i)
      if(normalize == 1){
        Ci_norm_nz <- Ci[ind_nz_i]/(sparse_diag_sqrt[ind_nz_i+i]*sparse_diag_sqrt[i])
      } else{
        Ci_norm_nz <- Ci[ind_nz_i]
      }
      sparseKernel[i,ind_nz_i+i] <- Ci_norm_nz
      sparseKernel[ind_nz_i+i,i] <- Ci_norm_nz
    }
    
    # Sparsify
    Cvals <- rep(0,N^2)
    rowInd <- rep(-1,N^2)
    colInd <- rep(-2,N^2)
    Cvals[1:N] <- diag(sparseKernel)
    rowInd[1:N] <- 1:N
    colInd[1:N] <- 1:N
    ctr <- N
    for(i in 1:(N-1)){
      Ci <- sparseKernel[i,(i+1):N]
      ind_nz_i <- which(Ci != 0)
      n_nz_i <- length(ind_nz_i)
      if(n_nz_i > 0){
        Ci_nz <- Ci[ind_nz_i]
        Cvals[(ctr+1):(ctr+2*n_nz_i)] <- c(Ci_nz,Ci_nz)
        rowInd[(ctr+1):(ctr+2*n_nz_i)] <- c(rep(i,n_nz_i),ind_nz_i+i)
        colInd[(ctr+1):(ctr+2*n_nz_i)] <- c(ind_nz_i+i,rep(i,n_nz_i))
        ctr <- ctr + 2*n_nz_i
      }
    }
    
    #=============== Core Kernel =================#
    Cvals[rowInd == colInd] <- (Cvals[rowInd == colInd] * (exp(log_sigma_vec)^2)) + exp(lognuggetSD)^2
    for(i in 1:(N-1)){
      # Ci <- Cvals[]
      ind_nz_i <- which(rowInd == i & colInd > i)
      ind_nz_i_v2 <- which(colInd == i & rowInd > i)
      if(length(ind_nz_i) > 0){
        colInd_rowI <- colInd[ind_nz_i]
        tmp11 <- Sigma11[colInd_rowI]
        tmp22 <- Sigma22[colInd_rowI]
        tmp12 <- Sigma12[colInd_rowI]
        cor_i <- nsCrosscorr( Xdist1_sq = matrix(dist1_sq[i,colInd_rowI], ncol = 1),
                              Xdist2_sq = matrix(dist2_sq[i,colInd_rowI], ncol = 1),
                              Xdist12 = matrix(dist12[i,colInd_rowI], ncol = 1),
                              Sigma11 = nimNumeric(Sigma11[i], length = 1), 
                              Sigma22 = nimNumeric(Sigma22[i], length = 1), 
                              Sigma12 = nimNumeric(Sigma12[i], length = 1),
                              PSigma11 = tmp11, PSigma22 = tmp22, PSigma12 = tmp12, 
                              nu = nu, d = d ) 
        cov_i <- nimNumeric(cor_i, length = length(colInd_rowI))*exp(log_sigma_vec[colInd_rowI])
        Cvals[ind_nz_i] <- Cvals[ind_nz_i] * cov_i * exp(log_sigma_vec[i])
        Cvals[ind_nz_i_v2] <- Cvals[ind_nz_i_v2] * cov_i * exp(log_sigma_vec[i])
      }
    }
    Cy_sm <- matrix(0,N^2,3) # cbind(rowInd, colInd, Cvals)
    Cy_sm[,1] <- rowInd
    Cy_sm[,2] <- colInd
    Cy_sm[,3] <- Cvals
    
    returnType(double(2))
    return(Cy_sm)
  }, check = FALSE, name = "Cy_sm_n11"
)


Csparse_dm <- nimbleFunction(  # Generate the sparse kernel in dense format
  run = function(
    dists = double(2), 
    coords = double(2), # coordinate / input locations
    N = double(0), # number of data points
    d = double(0), n1 = double(0), n2 = double(0), # dimension of input space; number of outer prods; number of bump fcns
    r0 = double(0), # scalar length-scale for compact stationary
    s0 = double(0), # Scalar; signal-variance of sparse stationary kernel.
    cstat_opt = double(0), # indicator for which compactly supported kernel is used
    normalize = double(0), # Logical; should C_sparse have 1's along the diagonal (1 = TRUE)
    bumpLocs = double(2), # array of bump function locations (n2*d x n1)
    rads = double(2), ampls = double(2), shps = double(2) # matrix of radii, amplitudes, shapes (n1 x n2)
    # dist1_sq = double(2), dist2_sq = double(2), dist12 = double(2), 
    # Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1), 
    # nu = double(0), log_sigma_vec = double(1), lognuggetSD = double(1)
  ) {
    
    #=============== Sparse Kernel =================#
    # Calculate the bump function sums + cross products
    for(j in 1:n1){
      g_j <- rep(0,N)
      for(b in 1:n2){
        # Calculate coords x bump fcn distances
        tmp <- rep(0,N)
        for(l in 1:d){
          rwix <- (n2*(l-1) + b)
          tmp <- tmp + (coords[1:N,l] - rep(bumpLocs[rwix,j],N))^2
        }
        dist_jb <- sqrt(tmp)
        ind_nz <- which(dist_jb <= rads[j,b])
        g_j[ind_nz] <- g_j[ind_nz] + ampls[j,b] * exp(shps[j,b] - shps[j,b]/(1 - dist_jb[ind_nz]^2/rads[j,b]^2))
        # if(ampls[j,b] < 0){
        #   a_ij <- 0
        # } else{
        #   a_ij <- 1
        # }
        # g_j[ind_nz] <- g_j[ind_nz] + a_ij * exp(shps[j,b] - shps[j,b]/(1 - dist_jb[ind_nz]^2/rads[j,b]^2))
      }
      if(j == 1){
        crsPrd <- g_j %*% asRow(g_j)
      } else{
        crsPrd <- crsPrd + g_j %*% asRow(g_j)
      }
    }
    
    # Add on the sparse stationary; normalize and store as sparse matrix
    k_c <- matrix(0,N,N)
    for(i in 1:N){
      k_c[i,i] <- s0
      if(i < N){
        tmpDst <- dists[i,(i+1):N]
        ind_dist <- which(tmpDst <= r0)
        # Choices for compactly supported kernel
        if(cstat_opt == 1){ # TPF(a=1, v=a+0.5+(d-1)/2)
          a <- 1; v <- a + 0.5 + (d-1)/2
          tmpK <- (1-(tmpDst[ind_dist]/r0)^a)^v
        }
        if(cstat_opt == 2){ # Wendland degree 2
          tmpK <- (1-abs(tmpDst[ind_dist]/r0))^4*(4*abs(tmpDst[ind_dist]/r0) + 1)
        }
        if(cstat_opt == 3){ # Wendland degree 3
          tmpK <- (1/3)*(1-abs(tmpDst[ind_dist]/r0))^6*(35*abs(tmpDst[ind_dist]/r0)^2 + 18*abs(tmpDst[ind_dist]/r0) + 3)
        }
        if(cstat_opt == 4){ # Wendland degree 4
          tmpK <- (1-abs(tmpDst[ind_dist]/r0))^8*(32*abs(tmpDst[ind_dist]/r0)^3 + 25*abs(tmpDst[ind_dist]/r0)^2 + 8*abs(tmpDst[ind_dist]/r0) + 1)
        }
        if(cstat_opt == 5){ # M and R (2009)
          tmpK <- ((3.0*(tmpDst[ind_dist]/r0)^2*log((tmpDst[ind_dist]/r0)/(1+sqrt(1.0 - (tmpDst[ind_dist]/r0)^2)))) +
                     ((2.0*(tmpDst[ind_dist]/r0)^2+1.0)*sqrt(1.0-(tmpDst[ind_dist]/r0)^2))) #* (sqrt(2.0)/(3.0*sqrt(pi)))
        }
        k_c[i,ind_dist+i] <- s0*tmpK
        k_c[ind_dist+i,i] <- s0*tmpK
        
      }
    }
    sparseKernel_unnorm <- crsPrd + k_c
    sparse_diag_sqrt <- sqrt(diag(sparseKernel_unnorm))
    sparseKernel <- matrix(0,N,N)
    if(normalize == 1){
      diag(sparseKernel) <- 1
    } else{
      diag(sparseKernel) <- diag(sparseKernel_unnorm)
    }
    for(i in 1:(N-1)){
      # Sparse
      Ci <- sparseKernel_unnorm[i,(i+1):N]
      ind_nz_i <- which(Ci != 0)
      # n_nz_i <- length(ind_nz_i)
      if(normalize == 1){
        Ci_norm_nz <- Ci[ind_nz_i]/(sparse_diag_sqrt[ind_nz_i+i]*sparse_diag_sqrt[i])
      } else{
        Ci_norm_nz <- Ci[ind_nz_i]
      }
      sparseKernel[i,ind_nz_i+i] <- Ci_norm_nz
      sparseKernel[ind_nz_i+i,i] <- Ci_norm_nz
    }

    returnType(double(2))
    return(sparseKernel)
  }, check = FALSE, name = "Csparse_dm"
)

Cy_dm <- nimbleFunction(  # Generate the sparse kernel in dense format
  run = function(
    dists = double(2), 
    coords = double(2), # coordinate / input locations
    N = double(0), # number of data points
    d = double(0), n1 = double(0), n2 = double(0), # dimension of input space; number of outer prods; number of bump fcns
    r0 = double(0), # scalar length-scale for compact stationary
    s0 = double(0), # Scalar; signal-variance of sparse stationary kernel.
    cstat_opt = double(0), # indicator for which compactly supported kernel is used
    normalize = double(0), # Logical; should C_sparse have 1's along the diagonal (1 = TRUE)
    bumpLocs = double(2), # array of bump function locations (n2*d x n1)
    rads = double(2), ampls = double(2), shps = double(2), # matrix of radii, amplitudes, shapes (n1 x n2)
    dist1_sq = double(2), dist2_sq = double(2), dist12 = double(2), 
    Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1), 
    nu = double(0), log_sigma_vec = double(1), lognuggetSD = double(1)
  ) {
    
    #=============== Sparse Kernel =================#
    # Calculate the bump function sums + cross products
    for(j in 1:n1){
      g_j <- rep(0,N)
      for(b in 1:n2){
        # Calculate coords x bump fcn distances
        tmp <- rep(0,N)
        for(l in 1:d){
          rwix <- (n2*(l-1) + b)
          tmp <- tmp + (coords[1:N,l] - rep(bumpLocs[rwix,j],N))^2
        }
        dist_jb <- sqrt(tmp)
        ind_nz <- which(dist_jb <= rads[j,b])
        g_j[ind_nz] <- g_j[ind_nz] + ampls[j,b] * exp(shps[j,b] - shps[j,b]/(1 - dist_jb[ind_nz]^2/rads[j,b]^2))
      }
      if(j == 1){
        crsPrd <- g_j %*% asRow(g_j)
      } else{
        crsPrd <- crsPrd + g_j %*% asRow(g_j)
      }
    }
    
    # Add on the sparse stationary; normalize and store as sparse matrix
    k_c <- matrix(0,N,N)
    for(i in 1:N){
      k_c[i,i] <- s0
      if(i < N){
        tmpDst <- dists[i,(i+1):N]
        ind_dist <- which(tmpDst <= r0)
        # Choices for compactly supported kernel
        if(cstat_opt == 1){ # TPF(a=1, v=a+0.5+(d-1)/2)
          a <- 1; v <- a + 0.5 + (d-1)/2
          tmpK <- (1-(tmpDst[ind_dist]/r0)^a)^v
        }
        if(cstat_opt == 2){ # Wendland degree 2
          tmpK <- (1-abs(tmpDst[ind_dist]/r0))^4*(4*abs(tmpDst[ind_dist]/r0) + 1)
        }
        if(cstat_opt == 3){ # Wendland degree 3
          tmpK <- (1/3)*(1-abs(tmpDst[ind_dist]/r0))^6*(35*abs(tmpDst[ind_dist]/r0)^2 + 18*abs(tmpDst[ind_dist]/r0) + 3)
        }
        if(cstat_opt == 4){ # Wendland degree 4
          tmpK <- (1-abs(tmpDst[ind_dist]/r0))^8*(32*abs(tmpDst[ind_dist]/r0)^3 + 25*abs(tmpDst[ind_dist]/r0)^2 + 8*abs(tmpDst[ind_dist]/r0) + 1)
        }
        if(cstat_opt == 5){ # M and R (2009)
          tmpK <- ((3.0*(tmpDst[ind_dist]/r0)^2*log((tmpDst[ind_dist]/r0)/(1+sqrt(1.0 - (tmpDst[ind_dist]/r0)^2)))) +
                     ((2.0*(tmpDst[ind_dist]/r0)^2+1.0)*sqrt(1.0-(tmpDst[ind_dist]/r0)^2))) #* (sqrt(2.0)/(3.0*sqrt(pi)))
        }
        k_c[i,ind_dist+i] <- s0*tmpK
        k_c[ind_dist+i,i] <- s0*tmpK
        
      }
    }
    sparseKernel_unnorm <- crsPrd + k_c
    sparse_diag_sqrt <- sqrt(diag(sparseKernel_unnorm))
    sparseKernel <- matrix(0,N,N)
    if(normalize == 1){
      diag(sparseKernel) <- 1
    } else{
      diag(sparseKernel) <- diag(sparseKernel_unnorm)
    }
    for(i in 1:(N-1)){
      # Sparse
      Ci <- sparseKernel_unnorm[i,(i+1):N]
      ind_nz_i <- which(Ci != 0)
      # n_nz_i <- length(ind_nz_i)
      if(normalize == 1){
        Ci_norm_nz <- Ci[ind_nz_i]/(sparse_diag_sqrt[ind_nz_i+i]*sparse_diag_sqrt[i])
      } else{
        Ci_norm_nz <- Ci[ind_nz_i]
      }
      sparseKernel[i,ind_nz_i+i] <- Ci_norm_nz
      sparseKernel[ind_nz_i+i,i] <- Ci_norm_nz
    }
    
    # # Sparsify
    # Cvals <- rep(0,N^2)
    # rowInd <- rep(-1,N^2)
    # colInd <- rep(-2,N^2)
    # Cvals[1:N] <- diag(sparseKernel)
    # rowInd[1:N] <- 1:N
    # colInd[1:N] <- 1:N
    # ctr <- N
    # for(i in 1:(N-1)){
    #   Ci <- sparseKernel[i,(i+1):N]
    #   ind_nz_i <- which(Ci != 0)
    #   n_nz_i <- length(ind_nz_i)
    #   if(n_nz_i > 0){
    #     Ci_nz <- Ci[ind_nz_i]
    #     Cvals[(ctr+1):(ctr+2*n_nz_i)] <- c(Ci_nz,Ci_nz)
    #     rowInd[(ctr+1):(ctr+2*n_nz_i)] <- c(rep(i,n_nz_i),ind_nz_i+i)
    #     colInd[(ctr+1):(ctr+2*n_nz_i)] <- c(ind_nz_i+i,rep(i,n_nz_i))
    #     ctr <- ctr + 2*n_nz_i
    #   }
    # }
    
    #=============== Core Kernel =================#
    Cy <- matrix(0,N,N)
    diag(Cy) <- (exp(log_sigma_vec)^2)*diag(sparseKernel) + exp(lognuggetSD)^2
    for(i in 1:(N-1)){
      # Ci <- Cvals[]
      # ind_nz_i <- which(rowInd == i & colInd > i)
      # ind_nz_i_v2 <- which(colInd == i & rowInd > i)
      Ci <- sparseKernel[i,(i+1):N]
      ind_nz_i <- which(Ci != 0)
      if(length(ind_nz_i) > 0){
        tmp11 <- Sigma11[ind_nz_i+i]
        tmp22 <- Sigma22[ind_nz_i+i]
        tmp12 <- Sigma12[ind_nz_i+i]
        cor_i <- nsCrosscorr( Xdist1_sq = matrix(dist1_sq[i,ind_nz_i+i], ncol = 1),
                              Xdist2_sq = matrix(dist2_sq[i,ind_nz_i+i], ncol = 1),
                              Xdist12 = matrix(dist12[i,ind_nz_i+i], ncol = 1),
                              Sigma11 = nimNumeric(Sigma11[i], length = 1), 
                              Sigma22 = nimNumeric(Sigma22[i], length = 1), 
                              Sigma12 = nimNumeric(Sigma12[i], length = 1),
                              PSigma11 = tmp11, PSigma22 = tmp22, PSigma12 = tmp12, 
                              nu = nu, d = d ) 
        cov_i <- nimNumeric(cor_i, length = length(ind_nz_i+i))*exp(log_sigma_vec[ind_nz_i+i])
        # if(any(is.na(Ci[ind_nz_i] * cov_i * exp(log_sigma_vec[i])))) stop(i)
        Cy[i,ind_nz_i+i] <- Ci[ind_nz_i] * cov_i * exp(log_sigma_vec[i])
        Cy[ind_nz_i+i,i] <- Ci[ind_nz_i] * cov_i * exp(log_sigma_vec[i])
      }
    }
    # Cy_sm <- matrix(0,N^2,3) # cbind(rowInd, colInd, Cvals)
    # Cy_sm[,1] <- rowInd
    # Cy_sm[,2] <- colInd 
    # Cy_sm[,3] <- Cvals
    
    returnType(double(2))
    return(Cy)
  }, check = FALSE, name = "Cy_dm"
)

# ROxygen comments ----
#' Calculate sparse kernel, core kernel, and determine nonzero entries
#'
#' \code{Cy_sm} calculates the normalized sparse kernel for a fixed
#' set of bump function hyperparameters and returns the nonzero entries. Note
#' that the matrix is calculated and returned in dense format.
#' 
#' @param dists N x N matrix of Euclidean distances
#' @param coords N x d matrix of coordinate/input locations
#' @param N Scalar; number of data measurements.
#' @param d Scalar; dimension of the spatial domain.
#' @param n1 Scalar; number of outer products.
#' @param n2 Scalar; number of bump functions in each outer product.
#' @param r0 Scalar; length-scale of sparse stationary kernel.
#' @param cstat_opt Scalar; determines the compactly supported kernel. See Details.
#' @param bumpLocs Array of bump function locations (n2*d x n1)
#' @param rads Matrix of bump function radii (n1 x n2; denoted \eqn{r_{ij}})
#' @param ampls Matrix of bump function amplitudes (n1 x n2; denoted \eqn{a_{ij}})
#' @param shps Matrix of bump function shape parameters (n1 x n2; denoted \eqn{b_{ij}})
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
#' @param log_sigma_vec Vector of length N; log of the signal standard deviation.
#' @param lognuggetSD Vector of length N; log of the error standard deviation.
#' 
#' 
#' @return Returns a sparse matrix (N x 3) of the nonzero elements of the product 
#' between the core and sparse kernel.
#'
#' @export
#' 
crossCy_sm <- nimbleFunction(  # Generate the sparse kernel in dense format
  run = function(
    Xdists = double(2), 
    coords = double(2), # coordinate / input locations
    Pcoords = double(2), # coordinate / input locations
    d = double(0), n1 = double(0), n2 = double(0), # dimension of input space; number of outer prods; number of bump fcns
    r0 = double(0), # scalar length-scale for compact stationary
    s0 = double(0), # Scalar; signal-variance of sparse stationary kernel.
    cstat_opt = double(0), # indicator for which compactly supported kernel is used
    normalize = double(0), # Logical; should C_sparse have 1's along the diagonal (1 = TRUE)
    bumpLocs = double(2), # array of bump function locations (n2*d x n1)
    rads = double(2), ampls = double(2), shps = double(2), # matrix of radii, amplitudes, shapes (n1 x n2)
    Xdist1_sq = double(2), Xdist2_sq = double(2), Xdist12 = double(2), 
    Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1), 
    PSigma11 = double(1), PSigma22 = double(1), PSigma12 = double(1), 
    nu = double(0), log_sigma_vec = double(1), Plog_sigma_vec = double(1)
  ) {
    
    N <- nrow(coords)
    M <- nrow(Pcoords)
    #=============== Sparse Kernel =================#
    # Calculate the bump function sums + cross products
    for(j in 1:n1){
      g_j <- rep(0,N)
      Pg_j <- rep(0,M)
      for(b in 1:n2){
        # Calculate coords x bump fcn distances
        tmp <- rep(0,N)
        Ptmp <- rep(0,M)
        for(l in 1:d){
          rwix <- (n2*(l-1) + b)
          tmp <- tmp + (coords[1:N,l] - rep(bumpLocs[rwix,j],N))^2
          Ptmp <- Ptmp + (Pcoords[1:M,l] - rep(bumpLocs[rwix,j],M))^2
        }
        dist_jb <- sqrt(tmp)
        ind_nz <- which(dist_jb <= rads[j,b])
        g_j[ind_nz] <- g_j[ind_nz] + ampls[j,b] * exp(shps[j,b] - shps[j,b]/(1 - dist_jb[ind_nz]^2/rads[j,b]^2))
        
        Pdist_jb <- sqrt(Ptmp)
        Pind_nz <- which(Pdist_jb <= rads[j,b])
        Pg_j[Pind_nz] <- Pg_j[Pind_nz] + ampls[j,b] * exp(shps[j,b] - shps[j,b]/(1 - Pdist_jb[Pind_nz]^2/rads[j,b]^2))
      }
      if(j == 1){
        crsPrd <- Pg_j %*% asRow(g_j)
        bf_diag <- g_j^2
        Pbf_diag <- Pg_j^2
      } else{
        crsPrd <- crsPrd + Pg_j %*% asRow(g_j)
        bf_diag <- bf_diag + g_j^2
        Pbf_diag <- Pbf_diag + Pg_j^2
      }
    }
    bf_diag <- bf_diag + s0
    Pbf_diag <- Pbf_diag + s0
    
    # Add on the sparse stationary; normalize and store as sparse matrix
    k_c <- matrix(0,M,N)
    for(i in 1:N){
      tmpDst <- Xdists[,i]
      ind_dist <- which(tmpDst <= r0)
      # Choices for compactly supported kernel
      if(cstat_opt == 1){ # TPF(a=1, v=a+0.5+(d-1)/2)
        a <- 1; v <- a + 0.5 + (d-1)/2
        tmpK <- (1-(tmpDst[ind_dist]/r0)^a)^v
      }
      if(cstat_opt == 2){ # Wendland degree 2
        tmpK <- (1-abs(tmpDst[ind_dist]/r0))^4*(4*abs(tmpDst[ind_dist]/r0) + 1)
      }
      if(cstat_opt == 3){ # Wendland degree 3
        tmpK <- (1/3)*(1-abs(tmpDst[ind_dist]/r0))^6*(35*abs(tmpDst[ind_dist]/r0)^2 + 18*abs(tmpDst[ind_dist]/r0) + 3)
      }
      if(cstat_opt == 4){ # Wendland degree 4
        tmpK <- (1-abs(tmpDst[ind_dist]/r0))^8*(32*abs(tmpDst[ind_dist]/r0)^3 + 25*abs(tmpDst[ind_dist]/r0)^2 + 8*abs(tmpDst[ind_dist]/r0) + 1)
      }
      if(cstat_opt == 5){ # M and R (2009)
        tmpK <- ((3.0*(tmpDst[ind_dist]/r0)^2*log((tmpDst[ind_dist]/r0)/(1+sqrt(1.0 - (tmpDst[ind_dist]/r0)^2)))) +
                   ((2.0*(tmpDst[ind_dist]/r0)^2+1.0)*sqrt(1.0-(tmpDst[ind_dist]/r0)^2))) #* (sqrt(2.0)/(3.0*sqrt(pi)))
      }
      k_c[ind_dist,i] <- tmpK*s0
    }
    sparseKernel_unnorm <- crsPrd + k_c
    # sparseKernel <- matrix(0,M,N)
    # for(i in 1:N){
    #   # Sparse
    #   Ci <- sparseKernel_unnorm[,i]
    #   ind_nz_i <- which(Ci != 0)
    #   # n_nz_i <- length(ind_nz_i)
    #   Ci_norm_nz <- Ci[ind_nz_i]/sqrt(bf_diag[i]*Pbf_diag[ind_nz_i])
    #   sparseKernel[ind_nz_i,i] <- Ci_norm_nz
    # }
    if(normalize == 1){
      sparseKernel <- diag(1/sqrt(Pbf_diag)) %*% sparseKernel_unnorm %*% diag(1/sqrt(bf_diag))
    }  else{
      sparseKernel <- sparseKernel_unnorm
    }
    
    # Sparsify
    Cvals <- rep(0,N*M)
    rowInd <- rep(-1,N*M)
    colInd <- rep(-2,N*M)
    ctr <- 0
    for(i in 1:N){
      Ci <- sparseKernel[,i]
      ind_nz_i <- which(Ci != 0)
      n_nz_i <- length(ind_nz_i)
      if(n_nz_i > 0){
        Ci_nz <- Ci[ind_nz_i]
        Cvals[(ctr+1):(ctr+n_nz_i)] <- Ci_nz
        rowInd[(ctr+1):(ctr+n_nz_i)] <- ind_nz_i # c(rep(i,n_nz_i),ind_nz_i+i)
        colInd[(ctr+1):(ctr+n_nz_i)] <- rep(i,n_nz_i) # c(ind_nz_i+i,rep(i,n_nz_i))
        ctr <- ctr + n_nz_i
      }
    }
    
    #=============== Core Kernel =================#
    # Cvals[rowInd == colInd] <- Cvals[rowInd == colInd] * exp(log_sigma_vec)^2 + exp(lognuggetSD)^2
    for(i in 1:N){
      ind_nz_i <- which(colInd == i)
      if(length(ind_nz_i) > 0){
        rowInd_colI <- rowInd[ind_nz_i]
        tmp11 <- PSigma11[rowInd_colI]
        tmp22 <- PSigma22[rowInd_colI]
        tmp12 <- PSigma12[rowInd_colI]
        cor_i <- nsCrosscorr( Xdist1_sq = matrix(Xdist1_sq[rowInd_colI,i], ncol = 1),
                              Xdist2_sq = matrix(Xdist2_sq[rowInd_colI,i], ncol = 1),
                              Xdist12 = matrix(Xdist12[rowInd_colI,i], ncol = 1),
                              Sigma11 = nimNumeric(Sigma11[i], length = 1), 
                              Sigma22 = nimNumeric(Sigma22[i], length = 1), 
                              Sigma12 = nimNumeric(Sigma12[i], length = 1),
                              PSigma11 = tmp11, PSigma22 = tmp22, PSigma12 = tmp12, 
                              nu = nu, d = d ) 
        cov_i <- nimNumeric(cor_i, length = length(rowInd_colI))*exp(Plog_sigma_vec[rowInd_colI])
        Cvals[ind_nz_i] <- Cvals[ind_nz_i] * cov_i * exp(log_sigma_vec[i])
      }
    }
    crossCy_sm <- matrix(0,N*M,3) # cbind(rowInd, colInd, Cvals)
    crossCy_sm[,1] <- rowInd
    crossCy_sm[,2] <- colInd
    crossCy_sm[,3] <- Cvals
    
    returnType(double(2))
    return(crossCy_sm)
  }, check = FALSE, name = "crossCy_sm"
)

crossCy_dm <- function(
    Xdists = double(2), 
    coords = double(2), # coordinate / input locations
    Pcoords = double(2), # coordinate / input locations
    d = double(0), n1 = double(0), n2 = double(0), # dimension of input space; number of outer prods; number of bump fcns
    r0 = double(0), # scalar length-scale for compact stationary
    s0 = double(0), # Scalar; signal-variance of sparse stationary kernel.
    cstat_opt = double(0), # indicator for which compactly supported kernel is used
    normalize = double(0), # Logical; should C_sparse have 1's along the diagonal (1 = TRUE)
    bumpLocs = double(2), # array of bump function locations (n2*d x n1)
    rads = double(2), ampls = double(2), shps = double(2), # matrix of radii, amplitudes, shapes (n1 x n2)
    Xdist1_sq = double(2), Xdist2_sq = double(2), Xdist12 = double(2), 
    Sigma11 = double(1), Sigma22 = double(1), Sigma12 = double(1), 
    PSigma11 = double(1), PSigma22 = double(1), PSigma12 = double(1), 
    nu = double(0), log_sigma_vec = double(1), Plog_sigma_vec = double(1)
) {
  
  N <- nrow(coords)
  M <- nrow(Pcoords)
  #=============== Sparse Kernel =================#
  # Calculate the bump function sums + cross products
  for(j in 1:n1){
    g_j <- rep(0,N)
    Pg_j <- rep(0,M)
    for(b in 1:n2){
      # Calculate coords x bump fcn distances
      tmp <- rep(0,N)
      Ptmp <- rep(0,M)
      for(l in 1:d){
        rwix <- (n2*(l-1) + b)
        tmp <- tmp + (coords[1:N,l] - rep(bumpLocs[rwix,j],N))^2
        Ptmp <- Ptmp + (Pcoords[1:M,l] - rep(bumpLocs[rwix,j],M))^2
      }
      dist_jb <- sqrt(tmp)
      ind_nz <- which(dist_jb <= rads[j,b])
      g_j[ind_nz] <- g_j[ind_nz] + ampls[j,b] * exp(shps[j,b] - shps[j,b]/(1 - dist_jb[ind_nz]^2/rads[j,b]^2))
      
      Pdist_jb <- sqrt(Ptmp)
      Pind_nz <- which(Pdist_jb <= rads[j,b])
      Pg_j[Pind_nz] <- Pg_j[Pind_nz] + ampls[j,b] * exp(shps[j,b] - shps[j,b]/(1 - Pdist_jb[Pind_nz]^2/rads[j,b]^2))
    }
    if(j == 1){
      crsPrd <- Pg_j %*% asRow(g_j)
      bf_diag <- g_j^2
      Pbf_diag <- Pg_j^2
    } else{
      crsPrd <- crsPrd + Pg_j %*% asRow(g_j)
      bf_diag <- bf_diag + g_j^2
      Pbf_diag <- Pbf_diag + Pg_j^2
    }
  }
  bf_diag <- bf_diag + s0
  Pbf_diag <- Pbf_diag + s0
  
  # Add on the sparse stationary; normalize and store as sparse matrix
  k_c <- matrix(0,M,N)
  for(i in 1:N){
    tmpDst <- Xdists[,i]
    ind_dist <- which(tmpDst <= r0)
    # Choices for compactly supported kernel
    if(cstat_opt == 1){ # TPF(a=1, v=a+0.5+(d-1)/2)
      a <- 1; v <- a + 0.5 + (d-1)/2
      tmpK <- (1-(tmpDst[ind_dist]/r0)^a)^v
    }
    if(cstat_opt == 2){ # Wendland degree 2
      tmpK <- (1-abs(tmpDst[ind_dist]/r0))^4*(4*abs(tmpDst[ind_dist]/r0) + 1)
    }
    if(cstat_opt == 3){ # Wendland degree 3
      tmpK <- (1/3)*(1-abs(tmpDst[ind_dist]/r0))^6*(35*abs(tmpDst[ind_dist]/r0)^2 + 18*abs(tmpDst[ind_dist]/r0) + 3)
    }
    if(cstat_opt == 4){ # Wendland degree 4
      tmpK <- (1-abs(tmpDst[ind_dist]/r0))^8*(32*abs(tmpDst[ind_dist]/r0)^3 + 25*abs(tmpDst[ind_dist]/r0)^2 + 8*abs(tmpDst[ind_dist]/r0) + 1)
    }
    if(cstat_opt == 5){ # M and R (2009)
      tmpK <- ((3.0*(tmpDst[ind_dist]/r0)^2*log((tmpDst[ind_dist]/r0)/(1+sqrt(1.0 - (tmpDst[ind_dist]/r0)^2)))) +
                 ((2.0*(tmpDst[ind_dist]/r0)^2+1.0)*sqrt(1.0-(tmpDst[ind_dist]/r0)^2))) #* (sqrt(2.0)/(3.0*sqrt(pi)))
    }
    k_c[ind_dist,i] <- s0*tmpK
  }
  sparseKernel_unnorm <- crsPrd + k_c
  if(normalize == 1){
    sparseKernel <- diag(1/sqrt(Pbf_diag)) %*% sparseKernel_unnorm %*% diag(1/sqrt(bf_diag))
  }  else{
    sparseKernel <- sparseKernel_unnorm
  }
  # sparseKernel <- matrix(0,M,N)
  # for(i in 1:N){
  #   # Sparse
  #   Ci <- sparseKernel_unnorm[,i]
  #   ind_nz_i <- which(Ci != 0)
  #   Ci_norm_nz <- Ci[ind_nz_i]/sqrt(bf_diag[i]*Pbf_diag[ind_nz_i])
  #   sparseKernel[ind_nz_i,i] <- Ci_norm_nz
  # }
  
  #=============== Core Kernel =================#
  crossCy_dm <- matrix(0,M,N)
  for(i in 1:N){
    Ci <- sparseKernel[,i]
    ind_nz_i <- which(Ci != 0)
    if(length(ind_nz_i) > 0){
      rowInd_colI <- ind_nz_i
      tmp11 <- PSigma11[rowInd_colI]
      tmp22 <- PSigma22[rowInd_colI]
      tmp12 <- PSigma12[rowInd_colI]
      cor_i <- nsCrosscorr( Xdist1_sq = matrix(Xdist1_sq[rowInd_colI,i], ncol = 1),
                            Xdist2_sq = matrix(Xdist2_sq[rowInd_colI,i], ncol = 1),
                            Xdist12 = matrix(Xdist12[rowInd_colI,i], ncol = 1),
                            Sigma11 = nimNumeric(Sigma11[i], length = 1), 
                            Sigma22 = nimNumeric(Sigma22[i], length = 1), 
                            Sigma12 = nimNumeric(Sigma12[i], length = 1),
                            PSigma11 = tmp11, PSigma22 = tmp22, PSigma12 = tmp12, 
                            nu = nu, d = d ) 
      cov_i <- nimNumeric(cor_i, length = length(rowInd_colI))*exp(Plog_sigma_vec[rowInd_colI])
      crossCy_dm[ind_nz_i,i] <- Ci[ind_nz_i] * cov_i * exp(log_sigma_vec[i])
    }
  }
  return(crossCy_dm)
}

#==============================================================================
# Density function for the gp2Scale kernel + sparse chol/solve
#==============================================================================

calc_num_nz <- nimbleFunction(  # Generate the sparse kernel in dense format
  run = function( Csparse = double(2) ) {
    Nnz <- sum(Csparse[,1] > 0)
    returnType(double())
    return(Nnz)
  }, name = "calc_num_nz"
)


# ROxygen comments ----
#' R_sparse_chol
#' @param i Vector of row indices.
#' @param j Vector of column indices.
#' @param x Vector of values in the matrix.
#' @export
R_sparse_cholesky <- function(i, j, x) {
  # Nnz <- sum(x != 0)
  # Asparse <- sparseMatrix(i = i[1:Nnz], j = j[1:Nnz], x = x[1:Nnz])
  Asparse <- sparseMatrix(i = i, j = j, x = x)
  ans.dsCMatrix <- chol(Asparse)
  ans.dgTMatrix <- as(ans.dsCMatrix, 'dgTMatrix')
  i <- ans.dgTMatrix@i + 1
  j <- ans.dgTMatrix@j + 1
  x <- ans.dgTMatrix@x
  ijx <- cbind(i, j, x)
  return(ijx)
}

# ROxygen comments ----
#' nimble_sparse_chol
#' @param i Vector of row indices.
#' @param j Vector of column indices.
#' @param x Vector of values in the matrix.
#' @export
nimble_sparse_cholesky <- nimbleRcall(
  prototype = function(i = double(1), j = double(1), x = double(1)) {},
  returnType = double(2),
  Rfun = 'R_sparse_cholesky'
)


# ROxygen comments ----
#' nimble_sparse_crossprod
#' @param i Vector of row indices.
#' @param j Vector of column indices.
#' @param x Vector of values in the matrix.
#' @param z Vector to calculate the cross-product with.
#' @param transp Optional indicator of using the transpose
#' @export
R_sparse_solveMat <- function(i, j, x, z, transp = 1) {
  if(transp == 0){ # No transpose
    Asparse <- sparseMatrix(i = i, j = j, x = x)
  } else{ # Transpose
    Asparse <- sparseMatrix(i = j, j = i, x = x)
  }
  ans.dsCMatrix <- solve(Asparse, as.numeric(z))
  return(ans.dsCMatrix@x)
}

# ROxygen comments ----
#' nimble_sparse_crossprod
#' @param i Vector of row indices.
#' @param j Vector of column indices.
#' @param x Vector of values in the matrix.
#' @param z Vector to calculate the cross-product with.
#' @param transp Optional indicator of using the transpose
#' @export
nimble_sparse_solveMat <- nimbleRcall(
  prototype = function(i = double(1), j = double(1), x = double(1), z = double(1), transp = double()) {},
  returnType = double(1),
  Rfun = 'R_sparse_solveMat'
)

# ROxygen comments ----
#' Function for the evaluating the Gaussian likelihood with gp2Scale sparse covariance.
#'
#' \code{dmnorm_gp2Scale} (and \code{rmnorm_gp2Scale}) calculate the usual Gaussian
#' likelihood for a fixed set of parameters (but with sparse matrices). Finally,
#' the distributions must be registered within \code{nimble}.
#' 
#' @param x Vector of measurements
#' @param mean Vector of mean valiues
#' @param Cov Matrix of size N x N; sparse kernel
#' @param N Number of measurements in x
#' @param Nnz Number of measurements in x
#' @param log Logical; should the density be evaluated on the log scale.
#' 
#' @return Returns the Gaussian likelihood using the gp2Scale sparse covariance.
#'
#' @export
#' 
dmnorm_gp2Scale <- nimbleFunction(
  run = function(x = double(1), mean = double(1), Cov = double(2), N = double(), Nnz = double(), log = double(0, default = 1)) {
    
    # Sparse cholesky
    Usp <- nimble_sparse_cholesky(i = Cov[1:Nnz,1], j = Cov[1:Nnz,2], x = Cov[1:Nnz,3])
    # Log determinant
    logdet_C <- sum(log(Usp[Usp[,1] == Usp[,2],3]))
    # Quadratic form
    Usp_t_inv_Z <- nimble_sparse_solveMat(i = Usp[,1], j = Usp[,2], x = Usp[,3], z = x - mean, transp = 1)
    qf <- sum(Usp_t_inv_Z^2)
    # Combine
    lp <- -logdet_C - 0.5*qf - 0.5*1.83787706649*N
    
    returnType(double())
    return(lp)
  }, check = FALSE, name = "dmnorm_gp2Scale"
)



# ROxygen comments ----
#' Function for the evaluating the SGV approximate density.
#'
#' \code{dmnorm_gp2Scale} (and \code{rmnorm_gp2Scale}) calculate the usual Gaussian
#' likelihood for a fixed set of parameters (but with sparse matrices). Finally,
#' the distributions must be registered within \code{nimble}.
#' 
#' @param x Vector of measurements
#' @param mean Vector of mean valiues
#' @param Csparse Matrix of size N x N; sparse kernel
#' @param Ccore Matrix of size N x N; core kernel
#' @param lognuggetSD Vector of log nugget (measurement error) standard deviation. 
#' @param N Number of measurements in x
#' @param log Logical; should the density be evaluated on the log scale.
#' 
#' @return Not applicable.
#'
#' @export
rmnorm_gp2Scale <- nimbleFunction(
  run = function(n = integer(), mean = double(1), Cov = double(2), N = double(), Nnz = double()) {
    returnType(double(1))
    return(numeric(N))
  }, name = "rmnorm_gp2Scale"
)

