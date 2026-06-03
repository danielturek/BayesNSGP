
#================================================================================================
# Functions for computing the joint density of the Nearest-Neighbor Gaussian Process (NNGP)
# Author: Fabian R. Ketwaroo
# Affiliation: Swiss Ornithological Institute
#================================================================================================


#================================================
# Exponential Covariance function
#================================================


#' Exponential covariance function for Gaussian processes
#'
#' Computes the exponential covariance matrix.
#'
#' @param edists Numeric matrix (\eqn{M \times M}) of pairwise Euclidean distances between spatial locations.
#' @param rho Positive numeric scalar. Length-scale parameter \eqn{\rho} controlling the decay of correlation with distance.
#' @param sigma2 Positive numeric scalar. Signal variance parameter \eqn{\sigma^2} controlling the overall variance.
#'
#' @return A numeric matrix (\eqn{M \times M}) representing the covariance between spatial locations.
#'
#' @details
#' The covariance function is defined as:
#' \deqn{C_{i,j} = \sigma^2 \exp(-d_{i,j} / \rho)}
#' where \eqn{d_{ij}} is the distance between locations \eqn{i} and \eqn{j}.
#'
#' Larger values of \eqn{\rho} imply slower decay of correlation with distance
#' (i.e., stronger spatial autocorrelation). Larger values of \eqn{\sigma^2}
#' increase the marginal variance of the process.
#'
#' A small nugget term (1e-6) is added to the diagonal for numerical stability.
#'
#' @examples
#' coords <- matrix(runif(10), ncol = 2)
#' d <- as.matrix(dist(coords))
#' expcov(edists = d, rho = 0.5, sigma2 = 1)
#'
#' @author Fabian Ketwaroo
#'
#' @export
expcov <- nimbleFunction(     
    run = function(edists = double(2), rho = double(0), sigma2 = double(0)) {
        returnType(double(2))
        m <- dim(edists)[1]
        n <- dim(edists)[2]
        
        result <- sigma2*exp(-edists/rho)
        
        # Only add nugget if the matrix is square (m x m)
        if (m == n) {
            nugget <- diag(m) * 1e-06 # add nugget effect for numerical stability
            return(result + nugget)
        }
        

        return(result)
    })


#================================================
# Compute Covariance matrix 
#================================================

#' Construct local covariance matrix for NNGP conditional distributions
#'
#' Constructs the local covariance matrix associated with a focal location
#' and its neighbor set under the Nearest-Neighbor Gaussian Process (NNGP).
#' 
#' @param fdist Numeric matrix (\eqn{k \times k}) containing pairwise
#'   Euclidean distances among the neighbors of the focal location.
#' 
#' @param ndist Numeric vector (length \eqn{k+1}) containing distances
#'   between the focal location and its neighbors, followed by the
#'   self-distance (0) of the focal location.
#'
#' @param rho Positive numeric scalar. Length-scale parameter \eqn{\rho} controlling
#'   the decay of spatial correlation.
#'
#' @param sigma2 Positive numeric scalar. Marginal variance parameter \eqn{\sigma^2}.
#'
#' @return A numeric matrix of dimension \eqn{(k+1) \times (k+1)} representing the
#'   local covariance matrix for the focal location and its neighbors.
#'
#' @details
#' This function constructs the local covariance matrix used in the NNGP
#' approximation for conditional Gaussian distributions. The covariance
#' structure is defined using an exponential covariance function.
#'
#' The structure of the matrix is:
#' \itemize{
#'   \item Top-left block: covariance among neighbors
#'   \item Last row/column: covariance between focal location and neighbors
#'   \item Bottom-right: marginal variance of the focal location
#' }
#'
#' A small nugget term (1e-6) is added to the diagonal for numerical stability.
#'
#' @author Fabian Ketwaroo
#'
#' @export
computeC <- nimbleFunction(
    run = function(fdist = double(2), ndist = double(1), rho = double(0), sigma2 = double(0) ){
        returnType(double(2))
        
        k1 <- length(ndist)-1
        C <- matrix(0,k1+1, k1+1 ) 
        C[1:k1, 1:k1] <- expcov(nimMatrix(fdist,k1,k1), rho, sigma2) # C value between neighbors of i 
        C[k1+1, 1:(k1+1)] <- sigma2*exp(-ndist/rho) # C value between i and its neighbors
        C[ k1+1, k1+1 ] <- C[ k1+1, k1+1 ] + 1e-06 # add nugget for computational stability
        return(C)
        
    }
)


#=========================================================================
# Compute NNGP regression coefficients (A) and conditional variances (D)
#=========================================================================

#' Compute NNGP regression coefficients and conditional variances
#'
#' Computes the Nearest-Neighbor Gaussian Process (NNGP)
#' factorization components for each location:
#' the regression coefficients (\eqn{A}) and conditional variances (\eqn{D})
#' following the algorithms described by Finley et al. (2019)
#'
#' @param edist Numeric matrix (\eqn{M \times M}). Pairwise Euclidean distances between spatial locations.
#'
#' @param nid.dist Numeric matrix (\eqn{M \times (k+1)}). Distances from each location to its
#'   ordered neighbor set, including self-distance in the last column.
#'
#' @param neighbors.id Integer matrix (\eqn{M \times k}). Indices of neighbor locations
#'   for each spatial location, based on the chosen ordering.
#'
#' @param rho Positive numeric scalar. Length-scale parameter \eqn{\rho}.
#'
#' @param sigma2 Positive numeric scalar. Marginal variance \eqn{\sigma^2}.
#'
#' @param k Positive integer. Maximum number of neighbors.
#'
#' @return A numeric matrix (\eqn{M \times (k+1)}) where:
#' \itemize{
#'   \item First \eqn{k} columns: regression coefficients (\eqn{A} matrix)
#'   \item Last column: conditional variances (\eqn{D} diagonal entries)
#' }
#'
#' @details
#' This function computes the parameters of the NNGP factorization
#' by expressing the joint distribution as a product of conditional Gaussian densities:
#'
#' \deqn{
#' p(\mathbf{x}) = \prod_{i=1}^M p(x_i \mid x_{N(i)})
#' }
#'
#' where \eqn{N(i)} denotes the set of neighbors of location \eqn{i}.
#'
#' For each location \eqn{i}, the conditional distribution is:
#'
#' \deqn{
#' x_i \mid x_{N(i)} \sim \mathcal{N}\left( A_i x_{N(i)}, \; D_i \right)
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{A_i} is a vector of regression coefficients
#'   \item \eqn{D_i} is the conditional variance
#' }
#'
#' These quantities are obtained from the covariance structure:
#'
#' \itemize{
#'   \item \eqn{A_i = C_{i,N(i)} C_{N(i),N(i)}^{-1}}
#'   \item \eqn{D_i = C_{i,i} - C_{i,N(i)} C_{N(i),N(i)}^{-1} C_{N(i),i}}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{C_{N(i),N(i)}} is the covariance matrix among neighbors
#'   \item \eqn{C_{i,N(i)}} is the covariance between location \eqn{i} and its neighbors
#' }
#'
#' The covariance structure is defined using an exponential covariance function.
#'
#' These parameters are used to efficiently evaluate the NNGP likelihood
#' and simulate from the process.
#' 
#' @references
#' Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E., & Banerjee, S. (2019). 
#' Efficient algorithms for Bayesian nearest neighbor Gaussian processes. 
#' *Journal of Computational and Graphical Statistics*, 28(2), 401–414. 
#' \doi{10.1080/10618600.2018.1537924}
#' 
#' @examples
#' coords <- matrix(runif(40), ncol = 2)
#' Nk = 5
#' res <- computeNeighbors(coords, k = Nk)
#' AD <- computeAD(edist = res$edist_sorted,
#'                 nid.dist = res$neighbors_dist,
#'                 neighbors.id = res$neighbor_idx,
#'                 rho = 0.1,
#'                 sigma2 = 0.3,
#'                 k = Nk)
#'
#' @author Fabian Ketwaroo
#' 
#' @export
computeAD <- nimbleFunction(
    run = function( edist = double(2), nid.dist = double(2), neighbors.id = double(2), rho = double(0), sigma2 = double(0), k = double(0)   ){
        returnType(double(2))
        
        # calculate A and D
        N <- dim(edist)[1]
        AD  <- matrix(0, N, k+1)
        AD[1, k+1] = sigma2 
        
        for (i in 1:(N-1)) {
            
            if(i < k) {
                
                
                k1 <- i
                
                # Select neighbors of i
                n.id <- neighbors.id[i+1, 1:k1 ]
                bdist <- edist[n.id, n.id] # distance between neighbors of i
                ndist <- nid.dist[i+1, 1:(k1+1)] # distance between i and its neighbors
                C <- computeC(bdist, ndist , rho, sigma2 )
                
                AD[i+1,  1:k1 ] = solve(C[1:k1, 1:k1], C[ i+1, 1:k1])
                AD[i+1, k+1] =  C[ i+1, k1+1] - inprod(C[i+1, 1:k1], AD[i+1, 1:k1])
                
                
            } else{
                
                k1 <- k
                
                # Select neighbors of i
                n.id <- neighbors.id[i+1, 1:k1 ]
                
                bdist <- edist[n.id, n.id] # distance between neighbors of i
                ndist <- nid.dist[i+1, 1:(k1+1)] # distance between i and its neighbors
                C <- computeC(bdist, ndist , rho, sigma2 )
                
                
                AD[i+1,  1:k1 ] = solve(C[1:k1, 1:k1], C[ k1+1, 1:k1])
                AD[i+1, k+1] =  C[ k1+1, k1+1] - inprod(C[k1+1, 1:k1], AD[i+1, 1:k1])
                
            } 
        }
        
        return(AD)
        
    }
)



#================================================
# Compute quadratic coefficient 
#================================================

#' Compute quadratic form for NNGP log-likelihood
#'
#' Evaluates the quadratic form appearing in the Gaussian log-likelihood
#' under the Nearest-Neighbor Gaussian Process (NNGP).
#'
#' @param u Numeric vector (length \eqn{M}). Typically centered data (\eqn{x - \mu}).
#'
#' @param v Numeric vector (length \eqn{M}). Typically identical to u in symmetric cases.
#'
#' @param AD Numeric matrix (\eqn{M \times (k+1)}). Output from \code{computeAD}, where:
#' \itemize{
#'   \item Columns \eqn{1:k} contain regression coefficients \eqn{A_i}
#'   \item Column \eqn{k+1} contains conditional variances \eqn{D_i}
#' }
#'
#' @param neighbors.id Integer matrix (\eqn{M \times k}). Neighbor indices for each location.
#'
#' @return A scalar value representing the quadratic form:
#' 
#'  \deqn{
#' \sum_{i=1}^M \frac{\left(u_i - A_i u_{N(i)}\right)
#' \left(v_i - A_i v_{N(i)}\right)}{D_i}
#' }
#'
#'
#' @details
#' This function computes the quadratic form associated with the NNGP precision
#' matrix without explicitly constructing the full covariance or precision matrix.
#' 
#' The expression is obtained by decomposing the joint Gaussian density into
#' a product of conditional densities:
#'
#' \deqn{
#' p(\mathbf{x}) = \prod_{i=1}^M p(x_i \mid x_{N(i)})
#' }
#' 
#' where \eqn{N(i)} denotes the set of neighbors of location \eqn{i}.
#'
#' Each term corresponds to a squared, standardized conditional residual:
#'
#' \deqn{
#' r_i = u_i - A_i u_{N(i)}
#' }
#'
#' and the quadratic form is:
#'
#' \deqn{
#' \sum_i \frac{r_i^2}{D_i}
#' }
#'
#' This provides an efficient way to evaluate the Gaussian log-likelihood under the NNGP approximation.
#'
#' @seealso computeAD
#'
#' @author Fabian Ketwaroo
#' 
#' @export
computeQF <- nimbleFunction(
    run = function(u = double(1), v = double(1), AD = double(2), neighbors.id = double(2) ){
        returnType(double(0))
        
        N <- dim(AD)[1]
        k <- dim(AD)[2]-1
        
        qf = u[1]*v[1]/AD[1,k+1]
        
        for(i in 2:N) { 
            if(i <= k) k1 <- i-1 else k1 <- k
            qf = qf + (u[i] - inprod(AD[i,1:k1], u[neighbors.id[i,1:k1]]))*(v[i] - inprod(AD[i, 1:k1 ], v[neighbors.id[i,1:k1] ]))/AD[i,k+1] 
        }
        
        return(qf)
        
    }
)


#================================================
# Efficiently simulate from NNGP
#================================================

#' Simulate from a Nearest-Neighbor Gaussian Process (NNGP)
#'
#' Efficiently simulates a realization from a Gaussian process using the
#' Nearest-Neighbor Gaussian Process (NNGP) approximation as
#' described by Datta (2022).
#' 
#' @param n Integer. Number of samples to generate. Currently only supports \code{n = 1}.
#'
#' @param mu Numeric vector (length \eqn{M}). Mean vector of the process.
#'
#' @param AD Numeric matrix (\eqn{M \times (k+1)}). Output from \code{computeAD}, where:
#' \itemize{
#'   \item Columns \eqn{1:k} contain regression coefficients (\eqn{A} matrix)
#'   \item Column \eqn{k+1} contains conditional variances (\eqn{D})
#' }
#'
#' @param neighbors.id Integer matrix (\eqn{M \times k}). Neighbor indices for each location,
#'   defining the conditioning set \eqn{N(i)}. Output from \code{computeNeighbors}.
#'
#' @return A numeric vector (length \eqn{M}) representing a simulated realization
#'   from the NNGP.
#'
#' @details
#' The simulation follows the recursive form derived from the sparse Cholesky 
#' factorization of the precision matrix (\eqn{\Sigma^{-1}}). Let \eqn{z = x - \mu} 
#' be the zero-mean spatial residuals. The NNGP simulates these residuals sequentially:
#'
#' \deqn{z_i = \sum_{j \in N(i)} A_{ij} z_j + \epsilon_i}
#'
#' where:
#' \itemize{
#'   \item \eqn{N(i)} is the neighbor set of location \eqn{i}, such that \eqn{j < i}.
#'   \item \eqn{A_{ij}} are the kriging weights computed from the local covariance.
#'   \item \eqn{\epsilon_i \sim \mathcal{N}(0, D_i)} are independent innovations.
#' }
#'
#' This sequential approach exploits the sparsity of the Cholesky factor \eqn{L}, 
#' where \eqn{\Sigma^{-1} \approx (I-A)^\top D^{-1} (I-A)} and allows for 
#' \eqn{O(Mk^3)} simulation, which is significantly faster than the standard 
#' \eqn{O(M^3)} Cholesky decomposition for large \eqn{M}.
#'
#' @references 
#' Datta, A. (2022). Nearest-neighbor sparse Cholesky matrices in spatial statistics. 
#' *Wiley Interdisciplinary Reviews: Computational Statistics*, 14(5), e1574.
#' \doi{10.1002/wics.1574}
#'
#' @examples
#' M <- 2000      # number of spatial locations
#' coords <- matrix(runif(2*M), ncol = 2)
#' Nk <- 15       # number of neighbors
#' res <- computeNeighbors(coords, k = Nk)
#' AD <- computeAD(edist = res$edist_sorted,
#'                 nid.dist = res$neighbors_dist,
#'                 neighbors.id = res$neighbor_idx,
#'                 rho = 0.1,
#'                 sigma2 = 0.3,
#'                 k = Nk)
#' w <- rmnorm_NN_GP(n = 1, mu = rep(0, M),
#'                  AD = AD[1:M, 1:(Nk+1)],
#'                  neighbors.id = res$neighbor_idx[1:M, 1:Nk])
#'
#' @author Fabian Ketwaroo
#' 
#' @export
rmnorm_NN_GP <- nimbleFunction(
    run = function(n = integer(0),  mu = double(1), AD = double(2), neighbors.id = double(2)) {
        returnType(double(1))
        
        M = dim(AD)[1]
        k = dim(AD)[2]-1
        A = AD[1:M, 1:k] # with elements b
        D = AD[1:M, k+1] # with elements f 
        
        # Auxiliary variable
        v = rnorm(M)
        z <- numeric(M) 
        
        z[1] <- v[1]*(D[1]^0.5)
        
        for (d in 2:k) {
            z[d] <- v[d]*(D[d]^0.5) + sum(A[d, (d-1):1]*z[(d-1):1] )
        }
        
        for (d in (k+1):M) {
            z[d] <- v[d]*(D[d]^0.5) + sum(A[d, ]*z[ neighbors.id[d,] ] )
        }
        
        return(z + mu)
        
    }
)




#================================================
# Compute NNGP joint density 
#================================================

#' Compute Nearest-Neighbor Gaussian Process (NNGP) log-density
#'
#' Evaluates the log-density of a multivariate normal distribution
#' under the Nearest-Neighbor Gaussian Process (NNGP) approximation
#' as defined by Data et al 2016.
#'
#' @param x Numeric vector (length \eqn{M}). Observed spatial process values.
#'
#' @param mu Numeric vector (length \eqn{M}). Mean vector.
#'
#' @param AD Numeric matrix (\eqn{M \times (k+1)}). Output from \code{computeAD}, where:
#' \itemize{
#'   \item Columns \eqn{1:k} contain regression coefficients \eqn{A_i}
#'   \item Column \eqn{k+1} contains conditional variances \eqn{D_i}
#' }
#'
#' @param neighbors.id Integer matrix (\eqn{M \times k}).  Neighbor indices defining
#'   the conditioning sets \eqn{N(i)}. Output from \code{computeNeighbors}
#'
#' @param log Logical/Integer. If \code{TRUE} (or 1), returns the log-density; 
#'   otherwise returns the density.
#'
#' @return A scalar value corresponding to the (log-)density.
#'
#' @details
#' The NNGP approximation factorizes the joint density as:
#'
#' \deqn{
#' p(\mathbf{x}) \approx \prod_{i=1}^M p(x_i \mid x_{N(i)})
#' }
#'
#' leading to the log-density:
#'
#' \deqn{
#' -\frac{1}{2} \left( N \log(2\pi) + \sum_i \log D_i + Q \right)
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{D_i} are conditional variances
#'   \item \eqn{Q} is the quadratic form computed by \code{computeQF}
#' }
#'
#' The quadratic form corresponds to:
#'
#' \deqn{
#' Q = \sum_{i=1}^M \frac{(x_i - \mu_i - A_i(x_{N(i)} - \mu_{N(i)}))^2}{D_i}
#' }
#'
#' This formulation avoids constructing the sparse covariance matrix,
#' enabling scalable likelihood evaluation for large spatial datasets.
#' 
#' @references 
#' Datta, A., Banerjee, S., Finley, A.O. and Gelfand, A.E., 2016. 
#' Hierarchical nearest-neighbor Gaussian process models for large geostatistical datasets. 
#' *Journal of the American Statistical Association*, 111(514), 800-812.
#' \doi{10.1080/01621459.2015.1044091}
#' 
#' @examples
#' M <- 2000      # number of spatial locations
#' coords <- matrix(runif(2*M), ncol = 2)
#' Nk <- 15       # number of neighbors
#' res <- computeNeighbors(coords, k = Nk)
#' AD <- computeAD(edist = res$edist_sorted,
#'                 nid.dist = res$neighbors_dist,
#'                 neighbors.id = res$neighbor_idx,
#'                 rho = 0.1,
#'                 sigma2 = 0.3,
#'                 k = Nk)
#' w <- rmnorm_NN_GP(n = 1, mu = rep(0, M),
#'                  AD = AD[1:M, 1:(Nk+1)],
#'                  neighbors.id = res$neighbor_idx[1:M, 1:Nk])
#' loglike <- dmnorm_NN_GP(w, mu = rep(0, M),
#'                        AD = AD[1:M, 1:(Nk+1)],
#'                        neighbors.id = res$neighbor_idx[1:M, 1:Nk],
#'                        log = TRUE)
#'
#' @author Fabian Ketwaroo
#' 
#' @export
dmnorm_NN_GP <- nimbleFunction(
    run = function( x= double(1), mu = double(1), AD = double(2), neighbors.id = double(2), log = integer(0, default =1) ){
        returnType(double(0))
        
        
        N = dim(AD)[1]
        M = dim(AD)[2]
        u = x- mu
        v = x - mu
        qf <- computeQF(u,v, AD, neighbors.id)
        
        loglik <- -0.5*(N*1.83787706649  +  sum(log(AD[1:N, M])) + qf )   # log(2pi) = 1.8378770664
        
        if(log)return(loglik)
        return(exp(loglik))
        
    }
)


registerDistributions(list(
    dmnorm_NN_GP = list(
        BUGSdist = 'dmnorm_NN_GP(mu, AD, neighbors.id)',
        types = c('value = double(1)', 'mu = double(1)', 'AD = double(2)', 'neighbors.id = double(2)'),
        mixedSizes = TRUE)
), verbose = FALSE)







#================================================
# Predictions from NNGP 
#================================================

#' Posterior predictive sampling for NNGP
#'
#' Generates a posterior predictive draw at a new location \eqn{x_0} using the
#' Nearest-Neighbor Gaussian Process (NNGP) approximation. This function 
#' automatically handles coordinate sorting, neighbor selection, and distance 
#' calculations internally following Algorithm 2 of Finley et al. (2019).
#'
#' @param x0 Numeric matrix (\eqn{1 \times 2}) or vector of length 2 representing 
#'   the coordinates of the new prediction location.
#' @param coords_sorted Numeric matrix (\eqn{M \times 2}) of the observed 
#'   training locations. Output from \code{computeNeighbors}.
#' @param rho Positive numeric scalar. Length-scale parameter \eqn{\rho}.
#' @param sigma2 Positive numeric scalar. Marginal variance parameter \eqn{\sigma^2}.
#' @param w.x0.all Numeric vector of length \eqn{M}. The realizations of the 
#'   spatial process (random effects) at all observed locations, ordered to 
#'   correspond with \code{coords_sorted}.
#' @param k Integer. The number of nearest neighbors to consider for the 
#'   prediction approximation.
#'
#' @return A numeric scalar representing a posterior predictive draw at location \eqn{x_0}.
#'
#' @details
#' The function merges the prediction point \eqn{x_0} into the training 
#' coordinates and identifies the \eqn{k} nearest neighbors from the set of 
#' points that precede it in the x-axis ordering. 
#' 
#' The predictive distribution is Gaussian:
#' \deqn{w(x_0) \mid \mathbf{w}_{N_0} \sim \mathcal{N}(m, v)}
#' where \eqn{m} is the **Kriging mean** and \eqn{v} is the **Kriging variance**
#' defined as:
#' \itemize{
#'   \item \eqn{m = c^\top C_{N_0}^{-1} \mathbf{w}_{N_0}}
#'   \item \eqn{v = C(x_0, x_0) - c^\top C_{N_0}^{-1} c}
#' }
#'
#' Here:
#' \itemize{
#'   \item \eqn{N_0} is the set of indices of the \eqn{k} nearest neighbors of \eqn{x_0} among the observed locations.
#'   \item \eqn{\mathbf{w}_{N_0}} is the vector of realizations of the spatial process at those neighbor locations (passed as \code{w.s.x0}).
#'   \item \eqn{c} is the \eqn{k \times 1} covariance vector between \eqn{x_0} and its neighbors in \eqn{N_0}.
#'   \item \eqn{C_{N_0}} is the \eqn{k \times k} covariance matrix among the neighbors in \eqn{N_0}.
#' }
#' 
#' The exponential covariance function is used: \eqn{C(d) = \sigma^2 \exp(-d / \rho)}.
#' A small nugget (\eqn{10^{-6}}) is added to the diagonal of \eqn{C_{N_0}} for 
#' numerical stability.
#' 
#' The internal subsetting of \code{w.x0.all} ensures that \eqn{\mathbf{w}_{N_0}} 
#' correctly represents the process realizations at the identified neighbor locations. 
#'
#' @references
#' Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E., & 
#' Banerjee, S. (2019). Efficient algorithms for Bayesian nearest neighbor 
#' Gaussian processes. *Journal of Computational and Graphical Statistics*, 
#' 28(2), 401–414. \doi{10.1080/10618600.2018.1537924}
#'
#' @examples
#' 
#' # Setup training data
#' M <- 50
#' coords <- matrix(runif(M * 2), ncol = 2)
#' w_all <- rnorm(M) # Mock spatial realizations
#' 
#' # Predict at a new location using k=5 neighbors
#' new_loc <- c(0.5, 0.5)
#' pred_draw <- NNGP.pred(x0 = new_loc, 
#'                        coords_sorted = coords, 
#'                        rho = 0.2, 
#'                        sigma2 = 1.0, 
#'                        w.x0.all = w_all,
#'                        k = 5)
#'
#' @author Fabian Ketwaroo
#' 
#' @export
NNGP.pred <- function( x0, coords_sorted, rho, sigma2, w.x0.all, k){
    
    Nk <- k
    M <- dim(coords_sorted)[1]
    
    # Add x0 to training locations
    coords.x0 = rbind(coords_sorted, x0) 
    
    # Sort coordinates by x axis 
    sorted_idx <- order(coords.x0[,1])
    x0.location = which(sorted_idx == M+1) 
    coords_sorted =coords.x0[sorted_idx, ]
    
    
    # Compute distance between all points
    coords.sf <- sf::st_as_sf(data.frame(x = coords_sorted[, 1], y = coords_sorted[, 2]), 
                              coords = c('x', 'y'))
    edist <- as.matrix(sf::st_distance(coords.sf))
    
    
    if(x0.location <= (Nk+1)){
        
        neighbors_idx.x0 <- as.numeric(1:(x0.location -  1)) # k nearest neighbors of x0 - N0
        neighbors_dist.x0 <-  edist[x0.location, neighbors_idx.x0] # dist between x0 and its neighbors
        neighbors_dist.N0 <- edist[neighbors_idx.x0, neighbors_idx.x0 ] # distance between neighbors of x0
        
    } else {
        
        curr_sorted_idx = order(edist[x0.location, ])
        curr_sorted_idx = curr_sorted_idx[curr_sorted_idx < x0.location] # take values to the left
        neighbors_idx.x0 = curr_sorted_idx[1:Nk]
        neighbors_dist.x0 <-  edist[x0.location, neighbors_idx.x0] # dist between x0 and its neighbors
        neighbors_dist.N0 <- edist[neighbors_idx.x0, neighbors_idx.x0 ] # distance between neighbors of x0
        
    }
    
    # Algorithm 2 of Finley et al. (2019).
    c = expcov(nimMatrix(neighbors_dist.x0), rho = rho, sigma2 = sigma2) # covariance between x0 and its k nearest neighbors
    if( any(x0.location== c(1,2) ) )  C_N0 = expcov(nimMatrix(neighbors_dist.N0) , rho = rho, sigma2 = sigma2) # covariance betwwen the k nearest neighbors
    else C_N0 = expcov(neighbors_dist.N0 , rho = rho, sigma2 = sigma2) + (diag( dim(neighbors_dist.N0 )[1] )*1e-06) # covariance betwwen the k nearest neighbors
    C_x0 = expcov(nimMatrix(0), rho = rho, sigma2 = sigma2) # covariance at prediction point x0
    
    w.x0 <- w.x0.all[neighbors_idx.x0]
    m = inprod(c, solve(C_N0 , w.x0) ) 
    v = C_x0 - inprod(c, solve(C_N0, c))
    a = rnorm(1, 0, 1)
    gp.pred = m + sqrt(v)*a
    return(gp.pred)
}



