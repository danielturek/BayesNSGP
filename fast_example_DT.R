##================================================
## Fast example of BayesNSGP
## Mark Risser and Daniel Turek
## March, 2019
##================================================

library(nimble)
library(BayesNSGP)

## Simulated data ================================

## Fixed values
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

## Setup
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

#######################
## User's Workflow   ##
## (as I see it) -DT ##
#######################

## fullGP likelihood
Rmodel <- nsgpModel(likelihood = 'fullGP', sigma_model = 'logLinReg', mu_model = 'linReg', coords = locs, z = z, X_sigma = Xmat, X_mu = Xmat)

## fullGP likelihood, with approxGP for tau model
Rmodel <- nsgpModel(likelihood = 'fullGP', sigma_model = 'logLinReg', mu_model = 'linReg', tau_model = 'approxGP', coords = locs, z = z, X_sigma = Xmat, X_mu = Xmat, tau_knot_coords = locs)

## NNGP likelihood
Rmodel <- nsgpModel(likelihood = 'NNGP', sigma_model = 'logLinReg', mu_model = 'linReg', coords = locs, z = z, X_sigma = Xmat, X_mu = Xmat, k = 10)

## SGV likelihood
Rmodel <- nsgpModel(likelihood = 'SGV', sigma_model = 'logLinReg', mu_model = 'linReg', coords = locs, z = z, X_sigma = Xmat, X_mu = Xmat, k = 10)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
## optionally modify samplers here

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
samples <- runMCMC(Cmcmc, niter = 5000)

pred <- nsgpPredict(Rmodel, samples, predCoords, PX_sigma = Xmat_pred, PX_mu = Xmat_pred)




