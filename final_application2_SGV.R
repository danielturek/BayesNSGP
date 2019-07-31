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
# Edit latent process samplers
knot_groups <- matrix(c(1:4,9:12,18:21,
                        27:30,36:39,45:48,
                        5:6,13:15,22:24,31:33,40:42,
                        7:8,16:17,25:26,34:35,43:44), ncol = 4)
# plot(coords, asp = 1.2)
# points(knot_coords[knot_groups[,1],], col = 2, pch = "+")
# points(knot_coords[knot_groups[,2],], col = 3, pch = "+")
# points(knot_coords[knot_groups[,3],], col = 4, pch = "+")
# points(knot_coords[knot_groups[,4],], col = 5, pch = "+")
conf$removeSamplers("w_sigma[1:48]")
for(h in 1:ncol(knot_groups)){
  conf$addSampler(target = c(paste0("w_sigma[",knot_groups[,h],"]")), type = "RW_block" )
}
conf$getSamplers()
conf$addMonitors("w_sigma")
# Build/run
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
prt <- proc.time()
samples_SGV <- runMCMC(Cmcmc, niter = 40000, nburnin = 0)
time_SGV <- proc.time() - prt
save(samples_SGV, time_SGV, file = "app2_SGV_samples_time.RData")



# Trace plots
par(ask=T)
for(h in 1:16) plot(samples_SGV[,h], type = "l", main = colnames(samples_SGV)[h])

# Calculate spatially-varying parameters ==========
postsamp <- samples_SGV[seq(from=20004,to=40000,by=4),]

constants <- Rmodel$isDataEnv$.BayesNSGP_constants_list

# p_tau <- constants$p_tau
# # tau_cross_dist_obs <- constants$tau_cross_dist
# tau_cross_dist_pred <- sqrt(nsCrossdist(Pcoords = predCoords, coords = constants$tau_knot_coords, isotropic = TRUE)$dist1_sq)
# tau_HP2 <- constants$tau_HP2
p_sigma <- constants$p_sigma
# sigma_cross_dist_obs <- constants$sigma_cross_dist
sigma_cross_dist_pred <- sqrt(nsCrossdist(Pcoords = predCoords, coords = constants$sigma_knot_coords, isotropic = TRUE)$dist1_sq)
sigma_HP2 <- constants$sigma_HP2
# X_Sigma <- constants$X_Sigma
PX_Sigma <- PXmat1
# X_mu <- constants$X_mu
PX_mu <- PXmat2
M <- nrow(CONUS_predDF)

subvec <- seq(50,5000,length=100)
Plog_sigma_vec <- PSigma11 <- PSigma12 <- PSigma22 <- Pmu <- matrix(NA, 100, M)
for(j in 1:100){
  cat(j, " ")
  samp_j <- postsamp[subvec[j],] 
  
  # w_tau_j <- as.numeric(samp_j[paste("w_tau[",1:p_tau,"]",sep = "")])
  
  # # Obs locations
  # Pmat_tau_obs_j <- matern_corr(tau_cross_dist_obs, samp_j["tauGP_phi"], tau_HP2)
  # log_tau_vec_j <- as.numeric(samp_j["tauGP_mu"]*rep(1,N) + samp_j["tauGP_sigma"] * Pmat_tau_obs_j %*% w_tau_j)
  
  # # Pred locations
  # Pmat_tau_pred_j <- matern_corr(tau_cross_dist_pred, samp_j["tauGP_phi"], tau_HP2)
  # Plog_tau_vec[j,] <- as.numeric(samp_j["tauGP_mu"]*rep(1,M) + samp_j["tauGP_sigma"] * Pmat_tau_pred_j %*% w_tau_j)
  
  # Required constants: p_sigma, sigma_cross_dist, sigma_cross_dist_pred, sigma_HP2
  w_sigma_j <- as.numeric(samp_j[paste("w_sigma[",1:p_sigma,"]",sep = "")])
  
  # # Obs locations
  # Pmat_sigma_obs_j <- matern_corr(sigma_cross_dist_obs, samp_j["sigmaGP_phi"], sigma_HP2)
  # log_sigma_vec_j <- as.numeric(samp_j["sigmaGP_mu"]*rep(1,N) + samp_j["sigmaGP_sigma"] * Pmat_sigma_obs_j %*% w_sigma_j)
  
  # Pred locations
  Pmat_sigma_pred_j <- matern_corr(sigma_cross_dist_pred, samp_j["sigmaGP_phi"], sigma_HP2)
  Plog_sigma_vec[j,] <- as.numeric(samp_j["sigmaGP_mu"]*rep(1,M) + samp_j["sigmaGP_sigma"] * Pmat_sigma_pred_j %*% w_sigma_j)
  
  # # Required constants: X_Sigma, PX_Sigma
  # eigen_comp1_j <- X_Sigma %*% samp_j[paste("Sigma_coef1[",1:ncol(X_Sigma),"]",sep = "")]
  # eigen_comp2_j <- X_Sigma %*% samp_j[paste("Sigma_coef2[",1:ncol(X_Sigma),"]",sep = "")]
  # eigen_comp3_j <- X_Sigma %*% samp_j[paste("Sigma_coef3[",1:ncol(X_Sigma),"]",sep = "")]
  # Sigma11_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 1))
  # Sigma12_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 3))
  # Sigma22_j <- as.numeric(inverseEigen(eigen_comp1_j, eigen_comp2_j, eigen_comp3_j, 2))
  
  Peigen_comp1_j <- PX_Sigma %*% samp_j[paste("Sigma_coef1[",1:ncol(PX_Sigma),"]",sep = "")]
  Peigen_comp2_j <- PX_Sigma %*% samp_j[paste("Sigma_coef2[",1:ncol(PX_Sigma),"]",sep = "")]
  Peigen_comp3_j <- PX_Sigma %*% samp_j[paste("Sigma_coef3[",1:ncol(PX_Sigma),"]",sep = "")]
  PSigma11[j,] <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 1))
  PSigma12[j,] <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 3)) 
  PSigma22[j,] <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 2))
  
  # mu <- as.numeric(X_mu %*% samp_j[paste("beta[",1:ncol(X_mu),"]",sep = "")])
  Pmu[j,] <- as.numeric(PX_mu %*% samp_j[paste("beta[",1:ncol(PX_mu),"]",sep = "")])
  
}

source("~/CASCADE_Projects/bayes_nsgp/paper_draft/plots_source_BayesNSGP.R")

plot_conus( df = CONUS_predDF, color = colMeans(Pmu), col.lim = c(-1, 1.4), 
            col.type = "disc", # scale.trans = "log", 
            col.pal = viridis_pal(option = "viridis")(9)[9:1], plot.grid = TRUE,
            brk_rd = 1, shp = 16, pt.size = 1.25, barheight = 20, n_brks = 10,
            xlab = NULL, ylab = NULL, ttle.txt = "Posterior mean: mu(s)" )
plot_conus( df = CONUS_predDF, color = colMeans(exp(Plog_sigma_vec)), col.lim = c(0.17, 1.9), 
            col.type = "disc", scale.trans = "log", 
            col.pal = viridis_pal(option = "viridis")(9)[9:1], plot.grid = TRUE,
            brk_rd = 2, shp = 16, pt.size = 1.25, barheight = 20, n_brks = 10,
            xlab = NULL, ylab = NULL, ttle.txt = "Posterior mean: sigma(s)" )
plot_conus( df = CONUS_predDF, color = colMeans(PSigma11), col.lim = c(5, 40), 
            col.type = "disc", # scale.trans = "log", 
            col.pal = viridis_pal(option = "viridis")(9)[9:1], plot.grid = TRUE,
            brk_rd = 1, shp = 16, pt.size = 1.25, barheight = 20, n_brks = 10,
            xlab = NULL, ylab = NULL, ttle.txt = "Posterior mean: Sigma11(s)" )
plot_conus( df = CONUS_predDF, color = colMeans(PSigma22), col.lim = c(5, 40), 
            col.type = "disc", # scale.trans = "log", 
            col.pal = viridis_pal(option = "viridis")(9)[9:1], plot.grid = TRUE,
            brk_rd = 1, shp = 16, pt.size = 1.25, barheight = 20, n_brks = 10,
            xlab = NULL, ylab = NULL, ttle.txt = "Posterior mean: Sigma22(s)" )
plot_conus( df = CONUS_predDF, color = colMeans(PSigma12), col.lim = c(-0.35, 0.35), 
            col.type = "disc", # scale.trans = "log", 
            col.pal = viridis_pal(option = "viridis")(9)[9:1], plot.grid = TRUE,
            brk_rd = 2, shp = 16, pt.size = 1.25, barheight = 20, n_brks = 10,
            xlab = NULL, ylab = NULL, ttle.txt = "Posterior mean: Sigma12(s)" )















