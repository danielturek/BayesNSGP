#================================================
# Application 2: fullGP and SGV
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
CONUSprecip <- read.csv("data/CONUS_WY2018.csv")

# Subset to western US
CONUSprecip <- CONUSprecip[CONUSprecip$longitude < -90,]
x_min <- min(CONUSprecip[,"longitude"])
x_max <- max(CONUSprecip[,"longitude"])
y_min <- min(CONUSprecip[,"latitude"])
y_max <- max(CONUSprecip[,"latitude"])

# Set up constants and design matrices
coords <- as.matrix(CONUSprecip[,c("longitude", "latitude")])
z <- CONUSprecip$logPR
elevShift <- mean(CONUSprecip$Xelevation)
elevScale <- sd(CONUSprecip$Xelevation)
lonShift <- mean(CONUSprecip$longitude)
lonScale <- sd(CONUSprecip$longitude)
CONUSprecip$Zelevation <- (CONUSprecip$Xelevation - elevShift)/elevScale
CONUSprecip$Zlongitude <- (CONUSprecip$longitude - lonShift)/lonScale
Xmat1 <- unname(lm(logPR ~ Zelevation, x = TRUE, data = CONUSprecip)$x)
Xmat2 <- unname(lm(logPR ~ Zelevation*Zlongitude, x = TRUE, data = CONUSprecip)$x)
N <- nrow(CONUSprecip)

# Set up the knot locations
knot_coords <- expand.grid(
  lon = seq(from = x_min - 0.25, to = x_max + 0.25, length = 9),
  lat = seq(from = y_min - 0.75, to = y_max, length = 6)
)
knot_coords <- as.matrix(knot_coords)
knot_coords <- knot_coords[-c(1:5,10),]

# plot(coords, asp = 1, pch = 16, cex = 0.25, ylim = c(25,50))
# US(add=T)
# points(knot_coords, col = 2, pch = "+")

# Load prediction grid: elevation and slope =====
CONUS_predDF <- read.csv("data/CONUS_predGrid.csv")
CONUS_predDF <- CONUS_predDF[CONUS_predDF$longitude < -90,]
CONUS_predDF$Zelevation <- (CONUS_predDF$Xelevation - elevShift)/elevScale
CONUS_predDF$Zlongitude <- (CONUS_predDF$longitude - lonShift)/lonScale
PXmat1 <- unname(lm(rnorm(nrow(CONUS_predDF)) ~ Zelevation, x = TRUE, data = CONUS_predDF)$x)
PXmat2 <- unname(lm(rnorm(nrow(CONUS_predDF)) ~ Zelevation*Zlongitude, x = TRUE, data = CONUS_predDF)$x)
predCoords <- CONUS_predDF[,c("longitude","latitude")]

# Constants for fullGP and SGV ==================
constants <- list( 
  nu = 0.5, k = 15, tau_HP1 = 10, sigma_knot_coords = knot_coords, sigma_HP1 = 10, 
  sigma_HP2 = 5, sigma_HP3 = 10, sigma_HP4 = 10, X_Sigma = Xmat1, Sigma_HP1 = 5, 
  maxAnisoDist = max(dist(coords)), X_mu = Xmat2, mu_HP1 = 10 )

#================================================
# MCMC using the fullGP likelihood
#================================================
Rmodel <- nsgpModel(likelihood = "fullGP", constants = constants, 
                    coords = coords, z = z, tau_model = "constant", 
                    sigma_model = "approxGP", mu_model = "linReg", 
                    Sigma_model = "compReg")
conf <- configureMCMC(Rmodel)
conf$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
conf$removeSamplers(c("Sigma_coef1[1]","Sigma_coef1[2]",
                      "Sigma_coef2[1]","Sigma_coef2[2]",
                      "Sigma_coef3[1]","Sigma_coef3[2]"))
conf$removeSamplers(c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"))
conf$addSampler(target = c("beta[1]","beta[2]","beta[3]","beta[4]"), type = "RW_block")
conf$addSampler(target = c("Sigma_coef1[1]","Sigma_coef1[2]","Sigma_coef2[1]",
                           "Sigma_coef2[2]","Sigma_coef3[1]","Sigma_coef3[2]"), type = "RW_block")
conf$addSampler(target = c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"), type = "RW_block")
conf$getSamplers()
conf$addMonitors("w_sigma")
# Build/run
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
prt <- proc.time()
samples_fullGP <- runMCMC(Cmcmc, niter = 40000, nburnin = 0)
time_fullGP <- proc.time() - prt
save(samples_fullGP, time_fullGP, file = "fullGP_samples_time.RData")

#================================================
# MCMC using the SGV likelihood
#================================================
Rmodel <- nsgpModel(likelihood = "SGV", constants = constants, 
                    coords = coords, z = z, tau_model = "constant", 
                    sigma_model = "approxGP", mu_model = "linReg", 
                    Sigma_model = "compReg")
conf <- configureMCMC(Rmodel)
conf$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
conf$removeSamplers(c("Sigma_coef1[1]","Sigma_coef1[2]",
                      "Sigma_coef2[1]","Sigma_coef2[2]",
                      "Sigma_coef3[1]","Sigma_coef3[2]"))
conf$removeSamplers(c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"))
conf$addSampler(target = c("beta[1]","beta[2]","beta[3]","beta[4]"), type = "RW_block")
conf$addSampler(target = c("Sigma_coef1[1]","Sigma_coef1[2]","Sigma_coef2[1]",
                           "Sigma_coef2[2]","Sigma_coef3[1]","Sigma_coef3[2]"), type = "RW_block")
conf$addSampler(target = c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"), type = "RW_block")
conf$getSamplers()
conf$addMonitors("w_sigma")
# Build/run
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
prt <- proc.time()
samples_SGV <- runMCMC(Cmcmc, niter = 40000, nburnin = 0)
time_SGV <- proc.time() - prt
save(samples_SGV, time_SGV, file = "SGV_samples_time.RData")













