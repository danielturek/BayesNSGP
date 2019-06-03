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
##
## Setup
##
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
Rmodel <- nsgpModel(likelihood = 'NNGP', sigma_model = 'logLinReg', mu_model = 'linReg', 
                    tau_model = 'approxGP', coords = locs, z = z, X_sigma = Xmat, X_mu = Xmat, 
                    tau_knot_coords = locs, k = 10)
Rmodel <- nsgpModel(likelihood = 'NNGP', tau_model = 'logLinReg', mu_model = 'linReg', 
                    sigma_model = 'approxGP', coords = locs, z = z, X_tau = Xmat, X_mu = Xmat, 
                    sigma_knot_coords = locs, k = 10, returnModelComponents = TRUE)
test <- nimbleModel(Rmodel$code, Rmodel$constants, Rmodel$data, Rmodel$inits)

## SGV likelihood
Rmodel <- nsgpModel(likelihood = 'SGV', sigma_model = 'logLinReg', mu_model = 'linReg', coords = locs, z = z, X_sigma = Xmat, X_mu = Xmat, k = 10)

Rmodel$calculate()
conf <- configureMCMC(Rmodel)
conf$printMonitors()
conf$printSamplers()
## optionally modify samplers here

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
samples <- runMCMC(Cmcmc, niter = 100)

pred <- nsgpPredict(Rmodel, samples, predCoords, PX_sigma = Xmat_pred, PX_mu = Xmat_pred)





########### MDR example
# Setup
N <- 100
set.seed(0)
coords <- matrix(runif(2*N), ncol = 2)
Xmat1 <- cbind(rep(1,N),locs[,1])
Xmat2 <- cbind(rep(1,N),locs[,2])
mu_vec <- as.numeric(Xmat2 %*% c(0, 2)) # Mean
alpha_vec <- as.numeric(Xmat1 %*% c(-0.5, 1.5)) # Log process SD
dist_list <- nsDist(coords)
Cor_mat <- nsCorr( dist1_sq = dist_list$dist1_sq, 
                   dist2_sq = dist_list$dist2_sq, dist12 = dist_list$dist12, 
                   Sigma11 = rep(0.4, N), Sigma22 = rep(0.4, N), 
                   Sigma12 = rep(0, N), nu = 0.5 )
Cov_mat <- diag(exp(alpha_vec)) %*% Cor_mat %*% diag(exp(alpha_vec))
D_mat <- diag(exp(rep(log(sqrt(0.05)), N))^2) 
# Draw data
set.seed(1)
z <- as.numeric(mu_vec + t(chol(Cov_mat + D_mat)) %*% rnorm(N))

constants <- list( X_sigma = Xmat1, X_Sigma = Xmat2, X_mu = Xmat1, tau_knot_coords = coords, k = 10 )

Rmodel <- nsgpModel( likelihood = "fullGP", 
                     sigma_model = "logLinReg", Sigma_model = "compRegIso",
                     mu_model = "linReg", tau_model = "approxGP", 
                     constants = constants, coords = coords, z = z )
conf <- configureMCMC(Rmodel)
conf$getSamplers()
conf$removeSamplers()
conf$addSampler(target = c("alpha[1]","alpha[2]"), type = "RW_block")
conf$addSampler(target = c("beta[1]","beta[2]"), type = "RW_block")
conf$addSampler(target = c("tauGP_mu","tauGP_phi","tauGP_sigma"), type = "RW_block")
conf$addSampler(target = c("Sigma_coef1[1]","Sigma_coef1[2]"), type = "AF_slice")
# First, subset the coordinates
groups <- list(which(coords[,1] < 0.5 & coords[,2] < 0.5))
groups[[2]] <- which(coords[,1] < 0.5 & coords[,2] >= 0.5)
groups[[3]] <- which(coords[,1] >= 0.5 & coords[,2] < 0.5)
groups[[4]] <- which(coords[,1] >= 0.5 & coords[,2] >= 0.5)
# Add block samplers
for(g in 1:4){
  conf$addSampler(target = paste0("w_tau[", groups[[g]], "]"), type = "RW_block")
}
conf$getSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
samples <- runMCMC(Cmcmc, niter = 5000, nburnin = 4000)


M <- 20^2
predCoords <- as.matrix(expand.grid(seq(0,1,length = sqrt(M)),seq(0,1,length = sqrt(M))))
Xmat1_pred <- cbind(rep(1,M), predCoords[,1])
Xmat2_pred <- cbind(rep(1,M), predCoords[,2])

pred_constants <- list( PX_sigma = Xmat1_pred, PX_Sigma = Xmat2_pred, PX_mu = Xmat1_pred )
pred <- nsgpPredict(model = Rmodel, samples = samples, 
                    coords.predict = predCoords, 
                    constants = pred_constants)




