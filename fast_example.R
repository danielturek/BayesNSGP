#================================================
# Fast example of BayesNSGP
# Mark Risser and Daniel Turek
# March, 2019
#================================================

library(nimble)
library(coda)
library(StatMatch)
source("BayesNSGP/R/core.R")
source("BayesNSGP/R/NNGP.R")
source("BayesNSGP/R/SGV.R")

# Simulated data ================================

# Fixed values
N <- 50
M <- 10^2
k <- 15
set.seed(8412)
locs <- matrix(runif(2*N), ncol = 2)
Xmat <- cbind(rep(1,N),locs[,2])
alpha_vec <- as.numeric(Xmat %*% c(-0.5, 1.5)) # Log process SD
delta_vec <- rep(log(sqrt(0.05)), N) # Log nugget SD
Sigma11_vec <- rep(0.4, N) # Kernel matrix element 1,1
Sigma22_vec <- rep(0.4, N) # Kernel matrix element 2,2
Sigma12_vec <- rep(0, N) # Kernel matrix element 1,2
mu_vec <- as.numeric(Xmat %*% c(0, 2)) # Mean
nu <- 0.5 # Smoothness

# Setup
dist_list <- nsDist(locs)
Cor_mat <- nsCorr( dist1_sq = dist_list$dist1_sq, dist2_sq = dist_list$dist2_sq, 
dist12 = dist_list$dist12, Sigma11 = Sigma11_vec, 
Sigma22 = Sigma22_vec, Sigma12 = Sigma12_vec, nu = nu )
Cov_mat <- diag(exp(alpha_vec)) %*% Cor_mat %*% diag(exp(alpha_vec))
D_mat <- diag(exp(delta_vec)^2) 
set.seed(110)
z <- as.numeric(mu_vec + t(chol(Cov_mat + D_mat)) %*% rnorm(N))
predCoords <- as.matrix(expand.grid(seq(0,1,length = sqrt(M)),seq(0,1,length = sqrt(M))))
Xmat_pred <- cbind(rep(1,M), predCoords[,2])

# Implementation ================================

# Constants setup
dist_fullGP <- nsDist(locs)
constants_fullGP <- list( dist1_sq = dist_fullGP$dist1_sq, dist2_sq = dist_fullGP$dist2_sq, 
                          dist12 = dist_fullGP$dist12, nu = 0.5, N = N, Sigma_HP1 = 2,
                          X_sigma = Xmat, p_sigma = ncol(Xmat),
                          X_mu = Xmat, p_mu = ncol(Xmat))
niter <- 500
strt.kp <- 401

# Defaults: tau_model = "constant", Sigma_model = "constant", mu_model = "constant"
Rmodel_fullGP <- nsgpModel( likelihood = "fullGP", sigma_model = "logLinReg", mu_model = "linReg",
                            constants = constants_fullGP, z = z )
compl_model_fullGP <- compileNimble( Rmodel_fullGP )
conf_model_fullGP <- configureMCMC( Rmodel_fullGP )
# IMPORTANT: allow the user to specify block samplers
conf_model_fullGP$removeSamplers( c("Sigma_coef1", "Sigma_coef2", "Sigma_coef3") )
conf_model_fullGP$addSampler( target = c("Sigma_coef1", "Sigma_coef2", "Sigma_coef3"), type = "RW_block" )
nim_mcmc_fullGP <- buildMCMC(conf_model_fullGP)
nim_Cmcmc_fullGP <- compileNimble(nim_mcmc_fullGP, project = Rmodel_fullGP)
nim_Cmcmc_fullGP$run(niter)

# Prediction
#### Note: these arguments need to be passed to the nsgpPredict() function. Using global assignment for now
X_sigma <<- Xmat
X_mu <<- Xmat
PX_sigma <<- Xmat_pred
PX_mu <<- Xmat_pred
postSamp_pred_fullGP <- nsgpPredict( Rmodel_fullGP, as.matrix(nim_Cmcmc_fullGP$mvSamples)[strt.kp:niter,], locs, predCoords )

