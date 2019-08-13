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
tasRV20_DF_all <- read.csv("data/C20C_DJFtasRV20_trend.csv")
tasRV20_DF_all <- tasRV20_DF_all[abs(tasRV20_DF_all$latitude) < 80,]
latShift <- mean(tasRV20_DF_all$latitude)
latScale <- sd(tasRV20_DF_all$latitude)
tasRV20_DF_all$Zlatitude <- (tasRV20_DF_all$latitude - latShift)/latScale
tasRV20_DF_all$rv20[tasRV20_DF_all$rv20 > 310] <- NA

# Regional lengthscale
region <- rep(NA, nrow(tasRV20_DF_all))
region[tasRV20_DF_all$latitude <= 0 & tasRV20_DF_all$longitude <= 0 & tasRV20_DF_all$longitude > -120] <- 1
region[tasRV20_DF_all$latitude <= 0 & (tasRV20_DF_all$longitude <= -120 | tasRV20_DF_all$longitude > 120)] <- 2
region[tasRV20_DF_all$latitude <= 0 & tasRV20_DF_all$longitude > 0 & tasRV20_DF_all$longitude <= 120] <- 3
region[tasRV20_DF_all$latitude > 0 & tasRV20_DF_all$longitude <= 0 & tasRV20_DF_all$longitude > -120] <- 4
region[tasRV20_DF_all$latitude > 0 & (tasRV20_DF_all$longitude <= -120 | tasRV20_DF_all$longitude > 120)] <- 5
region[tasRV20_DF_all$latitude >= 0 & tasRV20_DF_all$longitude > 0 & tasRV20_DF_all$longitude <= 120] <- 6
tasRV20_DF_all$region <- factor(region)

# Remove NA's for fitting
tasRV20_DF <- tasRV20_DF_all[!is.na(tasRV20_DF_all$rv20),]

# Subset (temporarily)
tasRV20_DF <- tasRV20_DF[sample(nrow(tasRV20_DF), 5000),]

# Standardize
z <- (tasRV20_DF$rv20 - 275)/20

# Design matrix
Xmat <- unname(lm(rv20 ~ region, x = TRUE, data = tasRV20_DF)$x)
N <- nrow(tasRV20_DF)

# Convert lon/lat to x/y/z
xyz.crds <- matrix(NA,nrow(tasRV20_DF),3)
# Transform degrees to radians
lat.radians <- tasRV20_DF$latitude*(pi/180)
lon.radians <- tasRV20_DF$longitude*(pi/180)
for(i in 1:nrow(xyz.crds)){
  xyz.crds[i,1] <- 6.371*cos(lat.radians[i])*cos(lon.radians[i]) # Earth diameter ~ 6371km
  xyz.crds[i,2] <- 6.371*cos(lat.radians[i])*sin(lon.radians[i])
  xyz.crds[i,3] <- 6.371*sin(lat.radians[i])
}
coords <- round(xyz.crds, 4)
# coords <- tasRV20_DF[,1:2]

# # Range
# library(geoR)
# range_vec <- rep(NA, 6)
# kappa_vec <- rep(NA, 6)
# sigma_vec <- rep(NA, 6)
# tau_vec <- rep(NA, 6)
# temp <- tempFitMat <- tempFitExp <- list()
# for(i in 1:6){
#   cat(i, " ")
#   crds <- coords[tasRV20_DF$region == i, ]
#   dat <- z[tasRV20_DF$region == i ]
#   temp[[i]] <- variog(coords = crds, data = dat, max.dist = 90)
#   tempFitMat[[i]] <- variofit(temp[[i]], cov.model = "matern", fix.kappa = TRUE, kappa = 5)
#   tempFitExp[[i]] <- variofit(temp[[i]], cov.model = "exponential")
# }
# par(ask = TRUE)
# for(i in 1:6){
#   plot(temp[[i]], main = i)
#   lines(tempFitExp[[i]], col = 2)
#   lines(tempFitMat[[i]], col = 4)
# }

# Constants for NNGP ============================
constants <- list( 
  nu = 5, k = 15, 
  X_sigma = Xmat, sigma_HP1 = 10,
  X_Sigma = Xmat, Sigma_HP1 = 10, 
  maxAnisoDist = 20 # maxDist = 22.07 km*1000
)

#================================================
# MCMC using the NNGP likelihood
#================================================
prt <- proc.time()
Rmodel <- nsgpModel(likelihood = "NNGP", constants = constants, 
                    coords = coords, z = z, tau_model = "constant", 
                    sigma_model = "logLinReg", mu_model = "constant", 
                    Sigma_model = "compRegIso")
conf <- configureMCMC(Rmodel)
conf$removeSamplers(c("alpha[1:6]","Sigma_coef1[1:6]"))
conf$addSampler(target = c("alpha[1]","alpha[2]","alpha[3]",
                           "Sigma_coef1[1]","Sigma_coef1[2]","Sigma_coef1[3]"), type = "RW_block")
conf$addSampler(target = c("alpha[4]","alpha[5]","alpha[6]",
                           "Sigma_coef1[4]","Sigma_coef1[5]","Sigma_coef1[6]"), type = "RW_block")
conf$getSamplers()
# Build
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
time_build <- proc.time() - prt
# Run
prt <- proc.time()
samples <- runMCMC(Cmcmc, niter = 7000, nburnin = 0)
time_mcmc <- proc.time() - prt
save(samples, time_mcmc, time_build, file = "app3_alt3.RData")

# par(ask=TRUE)
# for(h in 1:ncol(samples)) plot(samples[-(1:1000),h], type = "l", main = colnames(samples)[h])
# par(ask=FALSE)
