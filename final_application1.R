#================================================
# Application 1: 1981 CO precip, Paciorek/
#    Schervish vs. Risser/Calder analysis
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
COprecip <- read.csv("data/COprecip1981.csv")
coords <- as.matrix(COprecip[,c("Longitude", "Latitude")])
z <- COprecip$logPrecip
Xmat <- unname(lm(logPrecip ~ Zelevation*Zslope10, x = TRUE, data = COprecip)$x)
N <- nrow(COprecip)
elevMn <- mean(COprecip$Elevation)
elevScl <- sd(COprecip$Elevation)
slpMn <- mean(COprecip$Slope10)
slpScl <- sd(COprecip$Slope10)

# Set up the knot locations
x_min <- min(coords[,1])
x_max <- max(coords[,1])
y_min <- min(coords[,2])
y_max <- max(coords[,2])
N_knots <- 64
knot_coords <- expand.grid(
  lon = seq(from = x_min + 0.5*(x_max - x_min)/sqrt(N_knots), 
            to = x_max - 0.5*(x_max - x_min)/sqrt(N_knots), 
            length = sqrt(N_knots)),
  lat = seq(from = y_min + 0.5*(y_max - y_min)/sqrt(N_knots), 
            to = y_max - 0.5*(y_max - y_min)/sqrt(N_knots), 
            length = sqrt(N_knots)) )

# Load prediction grid: elevation and slope =====
COpredDF <- read.csv("data/COpredDF.csv")
COpredDF$Zelevation <- (COpredDF$elevation - elevMn)/elevScl
COpredDF$Zslope10 <- (COpredDF$slope - slpMn)/slpScl

# Subset
COpredDF_sub <- merge(
  x = expand.grid( longitude = unique(COpredDF$longitude)[seq(3,174,by = 4)], 
                   latitude = unique(COpredDF$latitude)[seq(5,108,by = 4)] ),
  y = COpredDF
)
predCoords <- COpredDF_sub[,c("longitude","latitude")]
Xmat_pred <- unname(lm(longitude ~ Zelevation*Zslope10, x = TRUE, data = COpredDF_sub)$x)

#================================================
# MCMC and prediction for Paciorek/Schervish 
#================================================
constants_PacScher <- list( nu = 2, Sigma_knot_coords = knot_coords,
                            Sigma_HP1 = c(10,10), Sigma_HP2 = rep(5,2), Sigma_HP3 = rep(3.85,2), 
                            Sigma_HP4 = c(10,20), maxAnisoDist = 16, mu_HP1 = 10 )

# Defaults: tau_model = "constant", sigma_model = "constant", mu_model = "constant"
Rmodel <- nsgpModel(likelihood = "fullGP", constants = constants_PacScher, 
                    coords = coords, z = z, Sigma_model = "npApproxGP")
conf <- configureMCMC(Rmodel)
conf$printSamplers()
# Edit latent process samplers
knot_groups <- matrix(c(1:4,9:12,17:20,25:28,
                        5:8,13:16,21:24,29:32,
                        32+c(1:4,9:12,17:20,25:28),
                        32+c(5:8,13:16,21:24,29:32)), ncol = 4)
# plot(coords, asp = 1.2)
# points(knot_coords[knot_groups[,1],], col = 2, pch = "+")
# points(knot_coords[knot_groups[,2],], col = 3, pch = "+")
# points(knot_coords[knot_groups[,3],], col = 4, pch = "+")
# points(knot_coords[knot_groups[,4],], col = 5, pch = "+")
conf$removeSamplers("w1_Sigma[1:64]")
conf$removeSamplers("w2_Sigma[1:64]")
conf$removeSamplers("w3_Sigma[1:64]")
for(h in 1:ncol(knot_groups)){
  conf$addSampler(target = c(paste0("w1_Sigma[",knot_groups[,h],"]")), type = "RW_block" )
  conf$addSampler(target = c(paste0("w2_Sigma[",knot_groups[,h],"]")), type = "RW_block" )
  conf$addSampler(target = c(paste0("w3_Sigma[",knot_groups[,h],"]")), type = "RW_block" )
}
# conf$getSamplers()
conf$addMonitors(c("Sigma11", "Sigma22", "Sigma12"))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
prt <- proc.time()
samples_PacScher <- runMCMC(Cmcmc, niter = 40000, nburnin = 20000)
time_PacScher <- proc.time() - prt

# Prediction
prt <- proc.time()
pred_PacScher <- nsgpPredict(model = Rmodel, 
                             samples = samples_PacScher[seq(5,20000,5),], 
                             coords.predict = predCoords)
time_PacScher_pred <- proc.time() - prt

save(samples_PacScher, time_PacScher, pred_PacScher, time_PacScher_pred, file = "app1_PacScher_mcmc_pred_time.RData")

#================================================
# MCMC and prediction for Risser/Calder 
#================================================
constants_RisserCalder <- list( nu = 0.5, X_mu = Xmat, mu_HP1 = 10,
                                X_sigma = Xmat, log_sigma_HP1 = 10,
                                X_Sigma = Xmat, Sigma_HP1 = c(10,10), 
                                Sigma_HP2 = c(2,2), maxAnisoDist = 16)

# Defaults: tau_model = "constant"
Rmodel <- nsgpModel(likelihood = "fullGP", constants = constants_RisserCalder, coords = coords, z = z, 
                    mu_model = "linReg", sigma_model = "logLinReg", Sigma_model = "covReg" )
conf <- configureMCMC(Rmodel)
conf$removeSamplers( c("psi11", "psi22", "rho") )
conf$addSampler( target = c("psi11", "psi22", "rho"), type = "RW_block" )
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
prt <- proc.time()
samples_RisserCalder <- runMCMC(Cmcmc, niter = 30000, nburnin = 10000)
time_RisserCalder <- proc.time() - prt

# Prediction
prt <- proc.time()
pred_RisserCalder <- nsgpPredict(model = Rmodel, samples = samples_RisserCalder[seq(5,20000,5),], 
                                 coords.predict = predCoords, 
                                 PX_sigma = Xmat_pred, PX_Sigma = Xmat_pred, PX_mu = Xmat_pred )
time_RisserCalder_pred <- proc.time() - prt

save(samples_RisserCalder, time_RisserCalder, pred_RisserCalder, time_RisserCalder_pred, file = "app1_RisserCalder_mcmc_pred_time.RData")


# library(mcmcse)
# summary(ess(samples_PacScher))
# plot(samples_PacScher[,"SigmaGP_mu[1]"], type = "l")
# plot(samples_PacScher[,"SigmaGP_mu[2]"], type = "l")
# plot(samples_PacScher[,"alpha"], type = "l")
# plot(samples_PacScher[,"beta"], type = "l")
# plot(samples_PacScher[,"delta"], type = "l")
# plot(samples_PacScher[,"w1_Sigma[21]"], type = "l")




