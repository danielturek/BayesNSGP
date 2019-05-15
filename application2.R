library(nimble)
library(coda)
library(StatMatch)
nimbleOptions(verbose = FALSE)
library(BayesNSGP)


# Load data =====================================
CONUSprecip <- read.csv("data/CONUS_WY2018.csv")

# # Testing: use only a subsample
# set.seed(1)
# subsamp <- sort(sample(nrow(CONUSprecip), 750))
# CONUSprecip <- CONUSprecip[subsamp,]
Csub <- CONUSprecip[CONUSprecip$longitude < -98,]
CONUSprecip <- Csub
x_min <- min(CONUSprecip[,"longitude"])
x_max <- max(CONUSprecip[,"longitude"])
y_min <- min(CONUSprecip[,"latitude"])
y_max <- max(CONUSprecip[,"latitude"])

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
# knot_coords <- expand.grid(
#   lon = seq(from = x_min + 0.5*(x_max - x_min)/7, to = x_max - 0.5*(x_max - x_min)/7, length = 7),
#   lat = seq(from = y_min + 0.5*(y_max - y_min)/6, to = y_max - 0.5*(y_max - y_min)/6, length = 6)
# )
# knot_coords <- knot_coords[-c(1:5,8),]
knot_coords <- expand.grid(
  lon = seq(from = x_min, to = x_max, length = 7),
  lat = seq(from = y_min, to = y_max, length = 6)
)
knot_coords <- as.matrix(knot_coords)
knot_coords <- knot_coords[-c(1:5,8),]

# plot(coords, asp = 1, pch = 16, cex = 0.25)
# points(knot_coords, col = 2, pch = "+")

# Load prediction grid: elevation and slope =====
CONUS_predDF <- read.csv("data/CONUS_predGrid.csv")

# MCMC setup ====================================
# Distance matrices and constants
constants <- list( nu = 0.5, 
  tau_knot_coords = knot_coords, tau_HP1 = 10, tau_HP2 = 5, tau_HP3 = 10, tau_HP4 = 10,
  sigma_knot_coords = knot_coords, sigma_HP1 = 10, sigma_HP2 = 5, sigma_HP3 = 10, sigma_HP4 = 10,
  X_Sigma = Xmat1, Sigma_HP1 = 5, Sigma_HP5 = 50, X_mu = Xmat2, mu_HP1 = 10 )

# Setup
Rmodel <- nsgpModel(likelihood = "fullGP", constants = constants, 
                    coords = coords, z = z, tau_model = "approxGP", 
                    sigma_model = "approxGP", mu_model = "linReg", 
                    Sigma_model = "compReg")
conf <- configureMCMC(Rmodel)
conf$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
conf$removeSamplers(c("Sigma_coef1[1]","Sigma_coef1[2]",
                      "Sigma_coef2[1]","Sigma_coef2[2]",
                      "Sigma_coef3[1]","Sigma_coef3[2]"))
conf$removeSamplers(c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"))
conf$removeSamplers(c("tauGP_mu","tauGP_phi","tauGP_sigma"))
conf$addSampler(target = c("beta[1]","beta[2]","beta[3]",
                           "Sigma_coef1[1]","Sigma_coef1[2]",
                           "Sigma_coef2[1]","Sigma_coef2[2]"), type = "RW_block")
conf$addSampler(target = c("Sigma_coef3[1]","Sigma_coef3[2]"), type = "RW_block")
conf$addSampler(target = c("beta[4]"), type = "RW_block")
conf$addSampler(target = c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"), type = "RW_block")
conf$addSampler(target = c("tauGP_mu","tauGP_phi","tauGP_sigma"), type = "RW_block")
conf$getSamplers()
conf$addMonitors(c("w_sigma", "w_tau"))

# Build/run
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
prt <- proc.time()
samples_fullGP <- runMCMC(Cmcmc, niter = 40000, nburnin = 0)
time_fullGP <- proc.time() - prt

save(samples_fullGP, time_fullGP, file = "application2_fullGP_full_long.RData")
# 
# # MCMC setup ====================================
# # Distance matrices and constants
# constants <- list( nu = 0.5, k = 15,
#                    tau_knot_coords = knot_coords, tau_HP1 = 10, tau_HP2 = 5, tau_HP3 = 10, tau_HP4 = 10,
#                    sigma_knot_coords = knot_coords, sigma_HP1 = 10, sigma_HP2 = 5, sigma_HP3 = 10, sigma_HP4 = 10,
#                    X_Sigma = Xmat1, Sigma_HP1 = 5, Sigma_HP5 = 50, X_mu = Xmat2, mu_HP1 = 10 )
# 
# # Setup
# Rmodel <- nsgpModel(likelihood = "SGV", constants = constants, 
#                     coords = coords, z = z, tau_model = "approxGP", 
#                     sigma_model = "approxGP", mu_model = "linReg", 
#                     Sigma_model = "compReg")
# conf <- configureMCMC(Rmodel)
# conf$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
# conf$removeSamplers(c("Sigma_coef1[1]","Sigma_coef1[2]",
#                       "Sigma_coef2[1]","Sigma_coef2[2]",
#                       "Sigma_coef3[1]","Sigma_coef3[2]"))
# conf$removeSamplers(c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"))
# conf$removeSamplers(c("tauGP_mu","tauGP_phi","tauGP_sigma"))
# conf$addSampler(target = c("beta[1]","beta[2]","beta[3]",
#                            "Sigma_coef1[1]","Sigma_coef1[2]",
#                            "Sigma_coef2[1]","Sigma_coef2[2]"), type = "RW_block")
# conf$addSampler(target = c("Sigma_coef3[1]","Sigma_coef3[2]"), type = "RW_block")
# conf$addSampler(target = c("beta[4]"), type = "RW_block")
# conf$addSampler(target = c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"), type = "RW_block")
# conf$addSampler(target = c("tauGP_mu","tauGP_phi","tauGP_sigma"), type = "RW_block")
# conf$getSamplers()
# conf$addMonitors(c("w_sigma", "w_tau"))
# 
# # Build/run
# Rmcmc <- buildMCMC(conf)
# Cmodel <- compileNimble(Rmodel)
# Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
# prt <- proc.time()
# samples_SGV <- runMCMC(Cmcmc, niter = 100, nburnin = 0)
# time_SGV <- proc.time() - prt
# 
# save(samples_SGV, time_SGV, file = "application2_SGV_full.RData")
# 



# # Trace plots
# load("application2_fullGP_full.RData")
# par(ask=T)
# for(h in 1:16) plot(samples_fullGP[,h], type = "l", main = colnames(samples_fullGP)[h])
# 
# ggcorr(samples_fullGP)#[-(1:1000),])
# 
# ggcorr(samples_fullGP[-(1:1000),], palette = "RdBu")






