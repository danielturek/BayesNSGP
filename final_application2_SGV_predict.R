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
predCoords <- as.matrix(CONUS_predDF[,c("longitude","latitude")])

# Constants for fullGP and SGV ==================
constants <- list( 
  nu = 0.5, k = 15, tau_HP1 = 10, sigma_knot_coords = knot_coords, sigma_HP1 = 10, 
  sigma_HP2 = 5, sigma_HP3 = 10, sigma_HP4 = 10, X_Sigma = Xmat1, Sigma_HP1 = 5, 
  maxAnisoDist = max(dist(coords)), X_mu = Xmat2, mu_HP1 = 10 )

#================================================
# MCMC using the SGV likelihood
#================================================
Rmodel <- nsgpModel(likelihood = "SGV", constants = constants, 
                    coords = coords, z = z, tau_model = "constant", 
                    sigma_model = "approxGP", mu_model = "linReg", 
                    Sigma_model = "compReg")
# conf <- configureMCMC(Rmodel)
# conf$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
# conf$removeSamplers(c("Sigma_coef1[1]","Sigma_coef1[2]",
#                       "Sigma_coef2[1]","Sigma_coef2[2]",
#                       "Sigma_coef3[1]","Sigma_coef3[2]"))
# conf$removeSamplers(c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"))
# conf$addSampler(target = c("beta[1]","beta[2]","beta[3]","beta[4]"), type = "RW_block")
# conf$addSampler(target = c("Sigma_coef1[1]","Sigma_coef1[2]","Sigma_coef2[1]",
#                            "Sigma_coef2[2]","Sigma_coef3[1]","Sigma_coef3[2]"), type = "RW_block")
# conf$addSampler(target = c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"), type = "RW_block")
# # Edit latent process samplers
# knot_groups <- matrix(c(1:4,9:12,18:21,
#                         27:30,36:39,45:48,
#                         5:6,13:15,22:24,31:33,40:42,
#                         7:8,16:17,25:26,34:35,43:44), ncol = 4)
# conf$removeSamplers("w_sigma[1:48]")
# for(h in 1:ncol(knot_groups)){
#   conf$addSampler(target = c(paste0("w_sigma[",knot_groups[,h],"]")), type = "RW_block" )
# }
# conf$getSamplers()
# conf$addMonitors("w_sigma")
# # Build/run
# Rmcmc <- buildMCMC(conf)
# Cmodel <- compileNimble(Rmodel)
# Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
# prt <- proc.time()
# samples_SGV <- runMCMC(Cmcmc, niter = 40000, nburnin = 0)
# time_SGV <- proc.time() - prt
# save(samples_SGV, time_SGV, file = "app2_SGV_samples_time.RData")
load("app2_SGV_samples_time.RData")

# Prediction
prt <- proc.time()
# pred_SGV <- nsgpPredict(model = Rmodel, samples = samples_SGV[30001:30100,], 
pred_SGV <- nsgpPredict(model = Rmodel, samples = samples_SGV[seq(from=30002,to=40000,by=2),], 
                        coords.predict = predCoords, PX_Sigma = PXmat1, PX_mu = PXmat2)
time_SGV_pred <- proc.time() - prt
save(pred_SGV, time_SGV_pred, file = "app2_SGV_pred.RData")

# plot_conus( df = CONUS_predDF, color = colMeans(pred_SGV$pred), col.lim = c(-1.8,2.4),
#             col.type = "disc", # scale.trans = "log",
#             col.pal = viridis_pal(option = "viridis")(9)[9:1], plot.grid = TRUE,
#             brk_rd = 1, shp = 16, pt.size = 1.25, barheight = 20, n_brks = 10,
#             xlab = NULL, ylab = NULL, ttle.txt = "Posterior mean: mu(s)" )


