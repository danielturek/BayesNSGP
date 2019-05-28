#================================================
# Application 3: NNGP
# Mark Risser
# Lawrence Berkeley National Laboratory
# May, 2019
#================================================

# Load packages
library(nimble)
library(coda)
library(StatMatch)
nimbleOptions(verbose = FALSE)
library(BayesNSGP)

# Load data =====================================
tasMax_trendDF <- read.csv("data/C20C_DJFtasMax_trend.csv")
latShift <- mean(tasMax_trendDF$latitude)
latScale <- sd(tasMax_trendDF$latitude)
tasMax_trendDF$Zlatitude <- (tasMax_trendDF$latitude - latShift)/latScale
# Convert lon/lat to x/y/z
xyz.crds <- matrix(NA,nrow(tasMax_trendDF),3)
# Transform degrees to radians
lat.radians <- tasMax_trendDF$latitude*(pi/180)
lon.radians <- tasMax_trendDF$longitude*(pi/180)
for(i in 1:nrow(xyz.crds)){
  xyz.crds[i,1] <- 6.371*cos(lat.radians[i])*cos(lon.radians[i]) # Earth diameter ~ 6371km
  xyz.crds[i,2] <- 6.371*cos(lat.radians[i])*sin(lon.radians[i])
  xyz.crds[i,3] <- 6.371*sin(lat.radians[i])
}
coords <- round(xyz.crds, 4)
z <- tasMax_trendDF$trendMax
Xmat <- unname(lm(trendMax ~ Zlatitude*ind_land, x = TRUE, 
                  data = tasMax_trendDF)$x)
N <- nrow(tasMax_trendDF)

# Constants for NNGP ============================
constants <- list( 
  nu = 0.5, k = 15, 
  X_tau = Xmat, tau_HP1 = 10, 
  X_sigma = Xmat, sigma_HP1 = 10,
  X_Sigma = Xmat, Sigma_HP1 = 10, 
  maxAnisoDist = 20, # maxDist = 22.07 km*1000
  X_mu = Xmat, mu_HP1 = 10 )

#================================================
# MCMC using the NNGP likelihood
#================================================
prt <- proc.time()
Rmodel <- nsgpModel(likelihood = "NNGP", constants = constants, 
                    coords = coords, z = z, tau_model = "logLinReg", 
                    sigma_model = "logLinReg", mu_model = "linReg", 
                    Sigma_model = "compRegIso")
conf <- configureMCMC(Rmodel)
conf$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
conf$addSampler(target = c("beta[1]","beta[2]","beta[3]","beta[4]"), type = "RW_block")
conf$removeSamplers(c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"))
conf$addSampler(target = c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"), type = "RW_block")
conf$removeSamplers(c("delta[1]","delta[2]","delta[3]","delta[4]"))
conf$addSampler(target = c("delta[1]","delta[2]","delta[3]","delta[4]"), type = "RW_block")
conf$removeSamplers(c("Sigma_coef1[1]","Sigma_coef1[2]","Sigma_coef1[3]","Sigma_coef1[4]"))
conf$addSampler(target = c("Sigma_coef1[1]","Sigma_coef1[2]","Sigma_coef1[3]","Sigma_coef1[4]"), type = "RW_block")
conf$getSamplers()
# Build
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
time_build <- proc.time() - prt
# Run
prt <- proc.time()
samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 0)
time_mcmc <- proc.time() - prt
save(samples, time_mcmc, time_build, file = "NNGP_samples_time.RData")

# par(ask=TRUE)
# for(h in 1:ncol(samples)) plot(samples[,h], type = "l", main = colnames(samples)[h])
# par(ask=FALSE)

#================================================
# MCMC using the NNGP likelihood
#================================================
prt <- proc.time()
Rmodel <- nsgpModel(likelihood = "SGV", constants = constants, 
                    coords = coords, z = z, tau_model = "logLinReg", 
                    sigma_model = "logLinReg", mu_model = "linReg", 
                    Sigma_model = "compRegIso")
conf <- configureMCMC(Rmodel)
conf$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
conf$addSampler(target = c("beta[1]","beta[2]","beta[3]","beta[4]"), type = "RW_block")
conf$removeSamplers(c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"))
conf$addSampler(target = c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"), type = "RW_block")
conf$removeSamplers(c("delta[1]","delta[2]","delta[3]","delta[4]"))
conf$addSampler(target = c("delta[1]","delta[2]","delta[3]","delta[4]"), type = "RW_block")
conf$removeSamplers(c("Sigma_coef1[1]","Sigma_coef1[2]","Sigma_coef1[3]","Sigma_coef1[4]"))
conf$addSampler(target = c("Sigma_coef1[1]","Sigma_coef1[2]","Sigma_coef1[3]","Sigma_coef1[4]"), type = "RW_block")
conf$getSamplers()
# Build
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
time_build <- proc.time() - prt
# Run
prt <- proc.time()
samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 0)
time_mcmc <- proc.time() - prt
save(samples, time_mcmc, time_build, file = "SGV_samples_time.RData")


