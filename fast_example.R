##================================================
## Fast example of BayesNSGP
## Mark Risser and Daniel Turek
## March, 2019
##================================================

library(nimble)
library(BayesNSGP)
nimbleOptions(verbose = FALSE)

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
data <- as.numeric(mu_vec + t(chol(Cov_mat + D_mat)) %*% rnorm(N))
predCoords <- as.matrix(expand.grid(seq(0,1,length = sqrt(M)),seq(0,1,length = sqrt(M))))
Xmat_pred <- cbind(rep(1,M), predCoords[,2])

#######################
## User's Workflow   ##
## (as I see it) -DT ##
#######################

## fullGP likelihood
Rmodel <- nsgpModel(likelihood = 'fullGP', sigma_model = 'logLinReg', mu_model = 'linReg', coords = locs, data = data, X_sigma = Xmat, X_mu = Xmat)

## fullGP likelihood, with approxGP for tau model
Rmodel <- nsgpModel(likelihood = 'fullGP', sigma_model = 'logLinReg', mu_model = 'linReg', tau_model = 'approxGP', coords = locs, data = data, X_sigma = Xmat, X_mu = Xmat, tau_knot_coords = locs)

## NNGP likelihood
Rmodel <- nsgpModel(likelihood = 'NNGP', sigma_model = 'logLinReg', mu_model = 'linReg', coords = locs, data = data, X_sigma = Xmat, X_mu = Xmat, k = 10)
Rmodel <- nsgpModel(likelihood = 'NNGP', sigma_model = 'logLinReg', mu_model = 'linReg', 
                    tau_model = 'approxGP', coords = locs, data = data, X_sigma = Xmat, X_mu = Xmat, 
                    tau_knot_coords = locs, k = 10)
# Rmodel <- nsgpModel(likelihood = 'NNGP', tau_model = 'logLinReg', mu_model = 'linReg', 
#                     sigma_model = 'approxGP', coords = locs, data = data, X_tau = Xmat, X_mu = Xmat, 
#                     sigma_knot_coords = locs, k = 10, returnModelComponents = TRUE)
# test <- nimbleModel(Rmodel$code, Rmodel$constants, Rmodel$data, Rmodel$inits)

## SGV likelihood
Rmodel <- nsgpModel(likelihood = 'SGV', sigma_model = 'logLinReg', mu_model = 'linReg', coords = locs, data = data, X_sigma = Xmat, X_mu = Xmat, k = 10)

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





########### MDR example: used for description in Section 5
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
data <- as.numeric(mu_vec + t(chol(Cov_mat + D_mat)) %*% rnorm(N))
tau_knot_coords  <- as.matrix(expand.grid(seq(0,1,length = 10),seq(0,1,length = 10)))

constants <- list( X_sigma = Xmat1, X_Sigma = Xmat2, X_mu = Xmat1, tau_knot_coords = tau_knot_coords, k = 10 )

Rmodel <- nsgpModel( likelihood = "SGV", 
                     sigma_model = "logLinReg", Sigma_model = "compRegIso",
                     mu_model = "linReg", tau_model = "approxGP", 
                     constants = constants, coords = coords, data = data )
conf <- configureMCMC(Rmodel)
conf$printSamplers()
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
conf$printSamplers()
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

########### MDR example: comparing computational times
Nvec <- c(50,100,200,500,1000,2000,5000,10000)
tau_knot_coords  <- as.matrix(expand.grid(seq(0,1,length = 10),seq(0,1,length = 10)))
fullGPtime <- NNGPtime <- SGVtime <- rep(NA, length(Nvec))
for(n in 1:length(Nvec)){
  N <- Nvec[n]
  print(N)
  coords <- matrix(5*runif(2*N), ncol = 2)
  Xmat1 <- cbind(rep(1,N),coords[,1])
  Xmat2 <- cbind(rep(1,N),coords[,2])
  constants <- list( X_sigma = Xmat1, X_Sigma = Xmat2, X_mu = Xmat1, tau_knot_coords = tau_knot_coords, k = 10 )
  
  # fullGP
  if(N <= 1000){
    Rmodel_fullGP <- nsgpModel( likelihood = "fullGP", 
                                sigma_model = "logLinReg", Sigma_model = "compRegIso",
                                mu_model = "linReg", tau_model = "approxGP", 
                                constants = constants, coords = coords, data = rnorm(N) ) #, returnModelComponents = TRUE )
    fullGPtime[n] <- system.time(Rmodel_fullGP$calculate())[3]
    rm(Rmodel_fullGP)
  }
  
  # NNGP
  Rmodel_NNGP <- nsgpModel( likelihood = "NNGP", 
                            sigma_model = "logLinReg", Sigma_model = "compRegIso",
                            mu_model = "linReg", tau_model = "approxGP", 
                            constants = constants, coords = coords, data = rnorm(N) )
  NNGPtime[n] <- system.time(Rmodel_NNGP$calculate())[3]
  rm(Rmodel_NNGP)
  
  # SGV
  Rmodel_SGV <- nsgpModel( likelihood = "SGV", 
                           sigma_model = "logLinReg", Sigma_model = "compRegIso",
                           mu_model = "linReg", tau_model = "approxGP", 
                           constants = constants, coords = coords, data = rnorm(N) )
  SGVtime[n] <- system.time(Rmodel_SGV$calculate())[3]
  rm(Rmodel_SGV)
  
}

timeDF <- data.frame(
  N = Nvec,
  exactGP = fullGPtime,
  SGV = SGVtime,
  NNGP = NNGPtime
)
#       N exactGP   SGV  NNGP
# 1    50   0.034 0.070 0.060
# 2   100   0.045 0.081 0.061
# 3   200   0.085 0.123 0.078
# 4   500   0.918 0.242 0.130
# 5  1000   3.965 0.421 0.209
# 6  2000      NA 0.969 0.370
# 7  5000      NA 3.892 0.965
# 8 10000      NA 7.780 1.724

plot(Nvec, log(fullGPtime), type = "b", pch = "+")
lines(Nvec, log(NNGPtime), type = "b", pch = "+", col = 2)
lines(Nvec, log(SGVtime), type = "b", pch = "+", col = 4)

plotDF <- data.frame(
  N = rep(Nvec, 3),
  Time = c(fullGPtime, SGVtime, NNGPtime),
  Likelihood = rep(c("Exact GP", "SGV", "NNGP"), each = length(Nvec))
)
library(ggplot2)
ggplot(plotDF, aes(x = N, y = Time, color = Likelihood)) + geom_point() + geom_line() +
  scale_y_continuous(trans="log", breaks=10^(-2:3), limits = c(0.01,100), name = "Time (s)") +
  scale_x_continuous(breaks = c(50,1000,2000,5000,10000))




