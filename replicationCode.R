#================================================
# Replication materials: BayesNSGP package for R
#================================================
# Mark Risser
# Lawrence Berkeley National Laboratory
# Daniel Turek
# Williams College
#================================================
# August 2019
#================================================

# Load required packages
library(nimble)
library(coda)
library(StatMatch)
nimbleOptions(verbose = FALSE)

# Install the BayesNSGP package from GitHub
devtools::install_github("danielturek/BayesNSGP", subdir = "BayesNSGP")
library(BayesNSGP)

#================================================
# Toy example: Section 5.1
#================================================

# Setup
N <- 100
set.seed(0)
coords <- matrix(runif(2*N), ncol = 2)
Xmat1 <- cbind(rep(1,N),coords[,1])
Xmat2 <- cbind(rep(1,N),coords[,2])
mu_vec <- as.numeric(Xmat2 %*% c(0, 2)) # Mean
alpha_vec <- as.numeric(Xmat1 %*% c(-0.5, 1.5)) # Log process SD
dist_list <- nsDist(coords)
Cor_mat <- nsCorr( dist1_sq = dist_list$dist1_sq,  dist2_sq = dist_list$dist2_sq, dist12 = dist_list$dist12, 
                   Sigma11 = rep(0.4, N), Sigma22 = rep(0.4, N), Sigma12 = rep(0, N), nu = 0.5 )
Cov_mat <- diag(exp(alpha_vec)) %*% Cor_mat %*% diag(exp(alpha_vec))
D_mat <- diag(exp(rep(log(sqrt(0.05)), N))^2) 
# Draw data
set.seed(1)
data <- as.numeric(mu_vec + t(chol(Cov_mat + D_mat)) %*% rnorm(N))
tau_knot_coords  <- as.matrix(expand.grid(seq(0,1,length = 10),seq(0,1,length = 10)))
# Define constants
constants <- list( X_sigma = Xmat1, X_Sigma = Xmat2, X_mu = Xmat1, tau_knot_coords = tau_knot_coords, k = 10 )
# Build the model
Rmodel <- nsgpModel( likelihood = "SGV", sigma_model = "logLinReg", Sigma_model = "compRegIso",
                     mu_model = "linReg", tau_model = "approxGP", constants = constants, coords = coords, data = data )
# Configure the model
conf <- configureMCMC(Rmodel) 
# Adjust the samplers
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
# Build/compile the model and the MCMC
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
# Run the MCMC
samples <- runMCMC(Cmcmc, niter = 5000, nburnin = 4000)
# Prediction: generate prediction locations and constants
M <- 20^2
predCoords <- as.matrix(expand.grid(seq(0,1,length = sqrt(M)),seq(0,1,length = sqrt(M))))
Xmat1_pred <- cbind(rep(1,M), predCoords[,1])
Xmat2_pred <- cbind(rep(1,M), predCoords[,2])
pred_constants <- list( PX_sigma = Xmat1_pred, PX_Sigma = Xmat2_pred, PX_mu = Xmat1_pred )
pred <- nsgpPredict(model = Rmodel, samples = samples, coords.predict = predCoords, constants = pred_constants)

#================================================
# Application 1: Colorado precipitation
#================================================
# Load data and setup
COprecip <- read.csv("data/COprecip1981.csv")
coords <- as.matrix(COprecip[,c("Longitude", "Latitude")])
data <- COprecip$logPrecip
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
# Load prediction grid: elevation and slope
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

# MCMC and prediction for Paciorek/Schervish 
constants_PacScher <- list( nu = 2, Sigma_knot_coords = knot_coords,
                            Sigma_HP1 = c(10,10), Sigma_HP2 = rep(5,2), Sigma_HP3 = rep(3.85,2), 
                            Sigma_HP4 = c(10,20), maxAnisoDist = 16, mu_HP1 = 10 )
# Defaults: tau_model = "constant", sigma_model = "constant", mu_model = "constant"
Rmodel <- nsgpModel(likelihood = "fullGP", constants = constants_PacScher, 
                    coords = coords, data = data, Sigma_model = "npApproxGP")
conf <- configureMCMC(Rmodel)
conf$printSamplers()
# Edit latent process samplers
knot_groups <- matrix(c(1:4,9:12,17:20,25:28,
                        5:8,13:16,21:24,29:32,
                        32+c(1:4,9:12,17:20,25:28),
                        32+c(5:8,13:16,21:24,29:32)), ncol = 4)
conf$removeSamplers("w1_Sigma[1:64]")
conf$removeSamplers("w2_Sigma[1:64]")
conf$removeSamplers("w3_Sigma[1:64]")
for(h in 1:ncol(knot_groups)){
  conf$addSampler(target = c(paste0("w1_Sigma[",knot_groups[,h],"]")), type = "RW_block" )
  conf$addSampler(target = c(paste0("w2_Sigma[",knot_groups[,h],"]")), type = "RW_block" )
  conf$addSampler(target = c(paste0("w3_Sigma[",knot_groups[,h],"]")), type = "RW_block" )
}
conf$addMonitors(c("Sigma11", "Sigma22", "Sigma12"))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
samples_PacScher <- runMCMC(Cmcmc, niter = 40000, nburnin = 20000)

# Prediction
pred_PacScher <- nsgpPredict(model = Rmodel, samples = samples_PacScher[seq(5,20000,5),], 
                             coords.predict = predCoords)

# MCMC and prediction for Risser/Calder 
constants_RisserCalder <- list( nu = 0.5, X_mu = Xmat, mu_HP1 = 10,
                                X_sigma = Xmat, log_sigma_HP1 = 10,
                                X_Sigma = Xmat, Sigma_HP1 = c(10,10), 
                                Sigma_HP2 = c(2,2), maxAnisoDist = 16)
# Defaults: tau_model = "constant"
Rmodel <- nsgpModel(likelihood = "fullGP", constants = constants_RisserCalder, coords = coords, data = data, 
                    mu_model = "linReg", sigma_model = "logLinReg", Sigma_model = "covReg" )
conf <- configureMCMC(Rmodel)
conf$removeSamplers( c("psi11", "psi22", "rho") )
conf$addSampler( target = c("psi11", "psi22", "rho"), type = "RW_block" )
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
samples_RisserCalder <- runMCMC(Cmcmc, niter = 30000, nburnin = 10000)
# Prediction
pred_RisserCalder <- nsgpPredict(model = Rmodel, samples = samples_RisserCalder[seq(5,20000,5),], 
                                 coords.predict = predCoords, 
                                 PX_sigma = Xmat_pred, PX_Sigma = Xmat_pred, PX_mu = Xmat_pred )

# PLOTS for Application 1

##### FIGURE 1 
pdf("COprecip_covariates.pdf", height = 5*1.1, width = 9*1.1)
grid.arrange(
  ggplot() +
    geom_polygon( data = map_data("state"), aes(x = long, y = lat, group = group),
                  color = "gray60", fill = "white", size = 0.2) +
    coord_fixed(ratio = 1.2) + theme_bw() + xlab("Longitude") + ylab("Latitude") + 
    geom_rect(aes(xmin = -110, xmax = -101, ymin = 36, ymax = 42), color = "red", fill = NA) +
    ggtitle("(a) Study area"),
  plot_colorado(COprecip, COprecip$logPrecip, pt.size = 1.25, shp = 16, brk_rd = 1,
                ttle.txt = "(b) Annual precipitation, 1981 (log mm)",
                col.type = "disc", n_brks = 10,
                col.pal = viridis_pal(option = "viridis")(9)[9:1],
                col.lim = c(5.2, 7.7) ) + theme_bw() + xlab(NULL) + ylab(NULL) +
    geom_point(data = data.frame( x = knot_coords[,1], y = knot_coords[,2], 
                                  Group = factor(c(rep(rep(1:2,each = 4),4), rep(rep(3:4,each = 4),4)))),
               mapping = aes(x = x, y = y, shape = Group), size = 1.65, alpha = 0.7) +
    scale_shape_manual(values = c(0,1,2,5)) + guides(shape = FALSE),
  plot_colorado(COpredDF, COpredDF$elevation/1000, brk_rd = 1, plot.grid = TRUE,
                col.type = "disc", n_brks = 10,
                ttle.txt = "(c) Elevation (km)", col.pal = terrain.colors(9)) +
    geom_point(data = COprecip, aes(x = longitude, y = latitude), size = 0.1) +
    theme_bw() + xlab(NULL) + ylab(NULL),
  plot_colorado(COpredDF, COpredDF$slope/1000, brk_rd = 1, plot.grid = TRUE, col.lim = c(-2.2,2.2),
                col.type = "disc", n_brks = 10, ttle.txt = "(d) Slope (km/deg longitude, west to east)", 
                col.pal = brewer.pal(9, "RdBu")) + 
    geom_point(data = COprecip, aes(x = longitude, y = latitude), size = 0.1) +
    theme_bw() + xlab(NULL) + ylab(NULL),
  layout_matrix = matrix(c(1,1,1,1,NA,2,2,2,2,2,
                           1,1,1,1,NA,2,2,2,2,2,
                           3,3,3,3,3,4,4,4,4,4,
                           3,3,3,3,3,4,4,4,4,4), nrow = 4, byrow = TRUE)
)
dev.off()

# Calculate spatial maps of the anisotropy processes
# Setup
x_min <- min(coords[,1])
x_max <- max(coords[,1])
y_min <- min(coords[,2])
y_max <- max(coords[,2])
N_knots <- 64
knot_coords <- expand.grid(
  lon = seq(from = x_min + 0.5*(x_max - x_min)/sqrt(N_knots), to = x_max - 0.5*(x_max - x_min)/sqrt(N_knots), 
            length = sqrt(N_knots)),
  lat = seq(from = y_min + 0.5*(y_max - y_min)/sqrt(N_knots), to = y_max - 0.5*(y_max - y_min)/sqrt(N_knots), 
            length = sqrt(N_knots))
)
dist_coords <- as.matrix(dist(coords, upper = T, diag = T))
Sigma_knot_dist <- as.matrix(dist(knot_coords, upper = T, diag = T))
Sigma_cross_dist <- mahalanobis.dist(coords, knot_coords, diag(2))

# ...and constants
constants_PacScher <- list( 
  nu = 2, N = N,
  Sigma_HP1 = c(10,10), Sigma_HP2 = rep(5,2), Sigma_HP3 = rep(3.85,2), 
  Sigma_HP4 = c(10,20), Sigma_HP5 = 16,
  dist = dist_coords, Sigma_cross_dist = Sigma_cross_dist, Sigma_knot_dist = Sigma_knot_dist,
  p_Sigma = N_knots, mu_HP1 = 10 )

p_Sigma <- constants_PacScher$p_Sigma
Sigma_cross_dist_obs <- constants_PacScher$Sigma_cross_dist_obs
Sigma_cross_dist_pred <- mahalanobis.dist(predCoords, knot_coords, diag(2))
Sigma_HP2 <- constants_PacScher$Sigma_HP2
M <- nrow(predCoords)
PSigma11 <- PSigma12 <- PSigma22 <- matrix(NA, M, 4000)
j_seq <- seq(1,20000,5)
for(j in 1:length(j_seq)){
  samp_j <- samples_PacScher[j_seq[j],]
  w1_Sigma_j <- samp_j[paste("w1_Sigma[",1:p_Sigma,"]",sep = "")]
  w2_Sigma_j <- samp_j[paste("w2_Sigma[",1:p_Sigma,"]",sep = "")]
  w3_Sigma_j <- samp_j[paste("w3_Sigma[",1:p_Sigma,"]",sep = "")]
  Pmat12_Sigma_pred_j <- matern_corr(Sigma_cross_dist_pred, samp_j["SigmaGP_phi[1]"], Sigma_HP2[1])
  Pmat3_Sigma_pred_j <- matern_corr(Sigma_cross_dist_pred, samp_j["SigmaGP_phi[2]"], Sigma_HP2[2])
  Peigen_comp1_j <- samp_j["SigmaGP_mu[1]"]*rep(1,M) + samp_j["SigmaGP_sigma[1]"] * Pmat12_Sigma_pred_j %*% w1_Sigma_j
  Peigen_comp2_j <- samp_j["SigmaGP_mu[1]"]*rep(1,M) + samp_j["SigmaGP_sigma[1]"] * Pmat12_Sigma_pred_j %*% w2_Sigma_j
  Peigen_comp3_j <- samp_j["SigmaGP_mu[2]"]*rep(1,M) + samp_j["SigmaGP_sigma[2]"] * Pmat3_Sigma_pred_j %*% w3_Sigma_j
  PSigma11[,j] <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 1))
  PSigma12[,j] <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 3))
  PSigma22[,j] <- as.numeric(inverseEigen(Peigen_comp1_j, Peigen_comp2_j, Peigen_comp3_j, 2))
  if(j %% 100 == 0) cat(j," ")
}

#### FIGURE 2(a), (b), (c)
pdf("postMean_Sigma.pdf", height = 4, width = 13)
grid.arrange(
  plot_colorado(COpredDF_sub, rowMeans(PSigma11), pt.size = 1.25, shp = 16, brk_rd = 2,
                ttle.txt = "(a) Posterior mean, Sigma11", scale.trans = "log",
                n_brks = 6, plot.grid = TRUE,
                col.pal = viridis_pal(option = "viridis")(9)[9:1],
                col.lim = c(0.005, 2.1) ) + theme_bw() + xlab(NULL) + ylab(NULL) +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_colorbar( barheight = 1, barwidth = 15 )),
  plot_colorado(COpredDF_sub, rowMeans(PSigma22), pt.size = 1.25, shp = 16, brk_rd = 2,
                ttle.txt = "(b) Posterior mean, Sigma22", scale.trans = "log",
                n_brks = 6, plot.grid = TRUE,
                col.pal = viridis_pal(option = "viridis")(9)[9:1],
                col.lim = c(0.005, 2.1) ) + theme_bw() + xlab(NULL) + ylab(NULL) +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_colorbar( barheight = 1, barwidth = 15 )),
  plot_colorado(COpredDF_sub, rowMeans(PSigma12), pt.size = 1.25, shp = 16, brk_rd = 2,
                ttle.txt = "(c) Posterior mean, Sigma12", 
                n_brks = 5, plot.grid = TRUE,
                col.pal = brewer.pal(9,"RdBu"),
                col.lim = c(-0.06, 0.06) ) + theme_bw() + xlab(NULL) + ylab(NULL) +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_colorbar( barheight = 1, barwidth = 15 )),
  nrow = 1
)
dev.off()

# Closest to knot locs
x_min <- min(coords[,1])
x_max <- max(coords[,1])
y_min <- min(coords[,2])
y_max <- max(coords[,2])
N_knots <- 24
knot_coords <- expand.grid(
  lon = seq(from = x_min + 0.5 * (x_max - x_min)/6, 
            to = x_max - 0.5 * (x_max - x_min)/6, 
            length = 6),
  lat = seq(from = y_min + 0.5 * (y_max - y_min)/4, 
            to = y_max - 0.5 * (y_max - y_min)/4, 
            length = 4)
)
Sigma_cross_dist2 <- mahalanobis.dist(predCoords, knot_coords, diag(2))
coord_knot_loc <- apply(Sigma_cross_dist2, 2, which.min)
Sigma_use <- cbind(rowMeans(PSigma11), rowMeans(PSigma12), rowMeans(PSigma22))[coord_knot_loc,]
coords_use <- as.matrix(predCoords[coord_knot_loc,])
Sigma_array <- array(NA, dim=c(2,2,N_knots))
for(i in 1:N_knots){
  Sigma_array[1,1,i] <- Sigma_use[i,1]
  Sigma_array[1,2,i] <- Sigma_use[i,2]
  Sigma_array[2,1,i] <- Sigma_use[i,2]
  Sigma_array[2,2,i] <- Sigma_use[i,3]
}

plotEllip <- function(locs, kernel_array, maintxt = NULL){
  mc.locations <- locs
  mc.kernels <- kernel_array
  K <- dim(mc.locations)[1]
  
  plot(ellipse(mc.kernels[, , 1], centre = mc.locations[1,], level = 0.5),
       type = "l", col = 2, asp = 1,
       xlab = "", ylab = "", ylim = c(36.5, 41.5), xlim = c(-109, -102), main = maintxt)
  points(mc.locations[1,1], mc.locations[1,2], cex = 1, pch="+" ) 
  
  for (k in 2:K) {
    lines(ellipse(mc.kernels[, , k], centre = mc.locations[k,], level = 0.5), col = 2 )
    points(mc.locations[k,1], mc.locations[k,2], cex = 1, pch="+" ) 
  }
}

#### FIGURE 2(d)
pdf("postMean_Sigma_ellipse.pdf", height = 5, width = 7)
plotEllip(coords_use, Sigma_array, maintxt = "(d) Posterior mean ellipses")
US(add=T, col = "gray")
points(coords, pch = 16, cex = 0.3, col = "gray50")
dev.off()

# Environmetrics results
oldResults <- data.frame(
  Parameter = c("beta0", "beta1", "beta2", "beta3", "sigmasq0", "alpha1", "alpha2", "alpha3",
                "gamma11", "gamma12", "gamma13", "gamma14", "gamma21", "gamma22", "gamma23", "gamma24",
                "psi11", "psi22", "psi12", "tausq" ),
  PostMean = c(6.308, 0.477, 0.053, 0.074, 0.163, 0.147, 0.101, -0.124,
               -0.292, 1.770, -0.752, 0.571, -0.312, -0.869, 0.876, 0.134,
               0.602, 1.240, -0.153, 0.010),
  CI_lb = c(6.155, 0.384, 0.013, 0.022, 0.115, -0.127, -0.089, -0.391,
            -1.429, 0.420, -1.723, -0.394, -1.278, -1.981, -0.687, -1.053,
            0.262, 0.535, -0.737, 0.006),
  CI_ub = c(6.488, 0.575, 0.091, 0.124, 0.227, 0.410, 0.308, 0.126, 
            0.756, 3.310, 0.200, 1.543, 0.596, 0.063, 2.034, 1.952,
            1.215, 2.671, 0.354, 0.016)
)

# Results for the current analysis
# Condition on gamma12 being positive
sampUse <- samples_RisserCalder[samples_RisserCalder[,"gamma1[2]"] >= 0,]
sampUse[,"alpha[1]"] <- exp(sampUse[,"alpha[1]"])^2 # Convert alpha[1] to sigmasq0
sampUse[,"rho"] <- sampUse[,"rho"]*sqrt(sampUse[,"psi11"]*sampUse[,"psi22"]) # Convert rho to psi12

newOrd <- c(5:8,1:4,10:20,9)
newResults <- data.frame(
  Parameter = c("beta0", "beta1", "beta2", "beta3", "sigmasq0", "alpha1", "alpha2", "alpha3",
                "gamma11", "gamma12", "gamma13", "gamma14", "gamma21", "gamma22", "gamma23", "gamma24",
                "psi11", "psi22", "psi12", "tausq" ),
  PostMean = apply(sampUse, 2, mean)[newOrd],
  CI_lb = apply(sampUse, 2, function(x){quantile(x, prob = 0.025)})[newOrd],
  CI_ub = apply(sampUse, 2, function(x){quantile(x, prob = 0.975)})[newOrd]
)

allResults <- rbind(oldResults, newResults)
allResults$Type <- ordered(rep(c("Risser/Calder","BayesNSGP"), each = 20), levels = c("Risser/Calder","BayesNSGP"))
rownames(allResults) <- NULL

#### FIGURE 3
pdf("compareEnvironmetrics.pdf", height = 4, width = 7)
ggplot(allResults, aes(x = Parameter, y = PostMean, ymin = CI_lb, ymax = CI_ub, color = Type)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.5) +
  scale_color_manual(values = c("red","blue"), name = "") + ylab("Posterior mean with 95% CrI") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

# Predictions/standard errors
plotDF <- as.data.frame(rbind(predCoords, predCoords))
plotDF$postMean <- c(apply(pred_PacScher$pred, 2, mean), apply(pred_RisserCalder$pred, 2, mean))
plotDF$postSD <- c(apply(pred_PacScher$pred, 2, sd), apply(pred_RisserCalder$pred, 2, sd))
plotDF$Type <- factor(rep(c("Paciorek & Schervish", "Risser & Calder"), each = nrow(predCoords) ), levels = c("Paciorek & Schervish", "Risser & Calder"))

#### FIGURE 4
pdf("postPrediction.pdf", height = 5.5, width = 9)
grid.arrange(
  plot_colorado(plotDF, plotDF$postMean, pt.size = 1.25, shp = 16, brk_rd = 2,
                ttle.txt = "(a) Posterior mean (log mm)", 
                n_brks = 5, plot.grid = TRUE,
                col.pal = viridis_pal(option = "viridis")(9)[9:1],
                col.lim = c(5, 7.2) ) + theme_bw() + xlab(NULL) + ylab(NULL) + facet_grid( ~ Type),
  plot_colorado(plotDF, plotDF$postSD, pt.size = 1.25, shp = 16, brk_rd = 2,
                ttle.txt = "(b) Posterior standard deviation (log mm)", scale.trans = "log",
                n_brks = 5, plot.grid = TRUE,
                col.pal = viridis_pal(option = "magma")(9)[9:1],
                col.lim = c(0.05, 0.6) ) + theme_bw() + xlab(NULL) + ylab(NULL) + facet_grid( ~ Type),
  nrow = 2
)
dev.off()


#================================================
# Application 2: CONUS precipitation
#================================================
# Load data 
CONUSprecip <- read.csv("data/CONUS_WY2018.csv")
CONUSprecip_all <- CONUSprecip

# Subset to western US
CONUSprecip <- CONUSprecip[CONUSprecip$longitude < -90,]
x_min <- min(CONUSprecip[,"longitude"])
x_max <- max(CONUSprecip[,"longitude"])
y_min <- min(CONUSprecip[,"latitude"])
y_max <- max(CONUSprecip[,"latitude"])

# Set up constants and design matrices
coords <- as.matrix(CONUSprecip[,c("longitude", "latitude")])
data <- CONUSprecip$logPR
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

# Load prediction grid: elevation and slope
CONUS_predDF <- read.csv("data/CONUS_predGrid.csv")
CONUS_predDF_all <- CONUS_predDF
CONUS_predDF <- CONUS_predDF[CONUS_predDF$longitude < -90,]
CONUS_predDF$Zelevation <- (CONUS_predDF$Xelevation - elevShift)/elevScale
CONUS_predDF$Zlongitude <- (CONUS_predDF$longitude - lonShift)/lonScale
PXmat1 <- unname(lm(rnorm(nrow(CONUS_predDF)) ~ Zelevation, x = TRUE, data = CONUS_predDF)$x)
PXmat2 <- unname(lm(rnorm(nrow(CONUS_predDF)) ~ Zelevation*Zlongitude, x = TRUE, data = CONUS_predDF)$x)
predCoords <- CONUS_predDF[,c("longitude","latitude")]

# Constants for fullGP and SGV
constants <- list( 
  nu = 0.5, k = 15, tau_HP1 = 10, sigma_knot_coords = knot_coords, sigma_HP1 = 10, 
  sigma_HP2 = 5, sigma_HP3 = 10, sigma_HP4 = 10, X_Sigma = Xmat1, Sigma_HP1 = 5, 
  maxAnisoDist = max(dist(coords)), X_mu = Xmat2, mu_HP1 = 10 )

# MCMC using the fullGP likelihood
Rmodel <- nsgpModel(likelihood = "fullGP", constants = constants, coords = coords, data = data, tau_model = "constant", 
                    sigma_model = "approxGP", mu_model = "linReg", Sigma_model = "compReg")
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
knot_groups <- matrix(c(1:4,9:12,18:21,27:30,36:39,45:48,
                        5:6,13:15,22:24,31:33,40:42,
                        7:8,16:17,25:26,34:35,43:44), ncol = 4)
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
samples_fullGP <- runMCMC(Cmcmc, niter = 40000, nburnin = 0)
time_fullGP <- proc.time() - prt
# Prediction
pred_fullGP <- nsgpPredict(model = Rmodel, samples = samples_fullGP[seq(from=30002,to=40000,by=2),], 
                        coords.predict = predCoords, PX_Sigma = PXmat1, PX_mu = PXmat2)

# MCMC using the SGV likelihood
Rmodel <- nsgpModel(likelihood = "SGV", constants = constants, coords = coords, data = data, tau_model = "constant", 
                    sigma_model = "approxGP", mu_model = "linReg", Sigma_model = "compReg")
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
knot_groups <- matrix(c(1:4,9:12,18:21,27:30,36:39,45:48,
                        5:6,13:15,22:24,31:33,40:42,
                        7:8,16:17,25:26,34:35,43:44), ncol = 4)
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
samples_SGV <- runMCMC(Cmcmc, niter = 40000, nburnin = 0)

# Prediction
pred_SGV <- nsgpPredict(model = Rmodel, samples = samples_SGV[seq(from=30002,to=40000,by=2),], 
                        coords.predict = predCoords, PX_Sigma = PXmat1, PX_mu = PXmat2)

# PLOTS for Application 2

#### Figure 5
pdf("CONUSprecip_elev.pdf", height = 5*1.6, width = 9*1.8/2)
grid.arrange(
  plot_conus( df = CONUSprecip_all, color = CONUSprecip_all$PR, col.lim = c(0.35, 11), 
              col.type = "disc", scale.trans = "log", 
              col.pal = viridis_pal(option = "viridis")(9)[9:1], plot.grid = FALSE,
              brk_rd = 1, shp = 16, pt.size = 1.25, barheight = 20, n_brks = 10,
              xlab = NULL, ylab = NULL, ttle.txt = "(a) Precipitation rate, 2018 water year (mm/day), with study region and knots" ) + 
    geom_rect(aes(xmin = -125.25, xmax = -90, ymin = 25.5, ymax = 49.5), color = "red", fill = NA, size = 0.5) +
    geom_point(data = data.frame( x = knot_coords[,1], y = knot_coords[,2], 
                                  Group = factor(knot_group_DF$grp) ),
               mapping = aes(x = x, y = y, shape = Group), size = 1.65, alpha = 0.7) +
    scale_shape_manual(values = c(0,1,2,5)) + guides(shape = FALSE),
  plot_conus( df = CONUS_predDF_all, color = CONUS_predDF_all$Xelevation/1000, col.lim = c(0, 2.5), 
              col.type = "disc",
              col.pal = terrain.colors(11), plot.grid = TRUE,
              brk_rd = 1, shp = 16, pt.size = 1.25, barheight = 20, n_brks = 12,
              xlab = NULL, ylab = NULL, ttle.txt = "(b) Elevation (km), with study region" ) +
    geom_point(data = CONUSprecip, mapping = aes(x = longitude, y = latitude), size = 0.05, alpha = 0.5) + 
    geom_rect(aes(xmin = -125.25, xmax = -90, ymin = 25.5, ymax = 49.5), color = "red", fill = NA, size = 0.5),
  nrow = 2
  # nrow = 1
)
dev.off()

sgvResults <- data.frame(
  Type = "SGV",
  Parameter = c("Aniso: x-dir intercept", "Aniso: x-dir elevation coef", "Aniso: y-dir intercept", 
                "Aniso: y-dir elevation coef", "Aniso: rot intercept", "Aniso: rot elevation coef", 
                "Mean: intercept", "Mean: elev main effect", "Mean: long main effect", "Mean: elev/long int.",
                "Nugget SD", "Spatial SD: lat GP mean", "Spatial SD: lat GP range", "Spatial SD: lat GP std dev"),
  PostMean = unname(apply(sampUse[,1:14], 2, mean)),
  CI_lb = unname(apply(sampUse[,1:14], 2, function(x){quantile(x, prob = 0.025)})),
  CI_ub = unname(apply(sampUse[,1:14], 2, function(x){quantile(x, prob = 0.975)}))
)
sampUse <- samples_fullGP[seq(from=30002,to=40000,by=2),]
fullGPResults <- data.frame(
  Type = "Exact GP",
  Parameter = c("Aniso: x-dir intercept", "Aniso: x-dir elevation coef", "Aniso: y-dir intercept", 
                "Aniso: y-dir elevation coef", "Aniso: rot intercept", "Aniso: rot elevation coef", 
                "Mean: intercept", "Mean: elev main effect", "Mean: long main effect", "Mean: elev/long int.",
                "Nugget SD", "Spatial SD: lat GP mean", "Spatial SD: lat GP range", "Spatial SD: lat GP std dev"),
  PostMean = unname(apply(sampUse[,1:14], 2, mean)),
  CI_lb = unname(apply(sampUse[,1:14], 2, function(x){quantile(x, prob = 0.025)})),
  CI_ub = unname(apply(sampUse[,1:14], 2, function(x){quantile(x, prob = 0.975)}))
)
app2_postMean <- rbind(fullGPResults,sgvResults)

##### FIGURE 6
pdf("app2_posteriors.pdf", height = 4, width = 7)
ggplot(app2_postMean, aes(x = Parameter, y = PostMean, ymin = CI_lb, ymax = CI_ub, color = Type)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.5) +
  scale_color_manual(values = c("red","blue"), name = "") + ylab("Posterior mean with 95% CrI") +
  theme(axis.text.x = element_text(angle=315, hjust=0.05)) + coord_cartesian(ylim = c(-4,3.25))
dev.off()


CONUS_predDF <- CONUS_predDF[CONUS_predDF$longitude < -90,]
CONUS_predDF$Zelevation <- (CONUS_predDF$Xelevation - elevShift)/elevScale
CONUS_predDF$Zlongitude <- (CONUS_predDF$longitude - lonShift)/lonScale
PXmat1 <- unname(lm(rnorm(nrow(CONUS_predDF)) ~ Zelevation, x = TRUE, data = CONUS_predDF)$x)
PXmat2 <- unname(lm(rnorm(nrow(CONUS_predDF)) ~ Zelevation*Zlongitude, x = TRUE, data = CONUS_predDF)$x)
predCoords <- as.matrix(CONUS_predDF[,c("longitude","latitude")])

postprdDF <- data.frame(
  meanSGV = c(colMeans(pred_SGV$pred), colMeans(pred_fullGP)),
  sdSGV = c(apply(pred_SGV$pred, 2, sd), apply(pred_fullGP, 2, sd)),
  likelihood = factor(rep(c("SGV", "Exact GP"), each = nrow(CONUS_predDF) ))
)
postprdDF <- cbind(postprdDF, rbind(CONUS_predDF,CONUS_predDF))

#### FIGURE 7
pdf("app2_postpred.pdf", height = 5*1.6, width = 9)
grid.arrange(
  plot_conus( df = postprdDF, color = postprdDF$meanSGV, col.lim = c(-2.6, 2.4), 
              col.type = "cont", plot.grid = TRUE,
              col.pal = viridis_pal(option = "viridis")(9)[9:1],
              brk_rd = 1, barheight = 13, n_brks = 8,
              xlab = NULL, ylab = NULL, 
              ttle.txt = "(a) Posterior mean (log mm per day)" ) + theme_bw() +
    coord_cartesian(xlim = c(-125,-89)) + facet_grid( . ~ likelihood),
  plot_conus( df = postprdDF, color = postprdDF$sdSGV, col.lim = c(0.01, 0.8), 
              col.type = "cont", plot.grid = TRUE, scale.trans = "log",
              col.pal = viridis_pal(option = "magma")(9)[9:1],
              brk_rd = 2, barheight = 13, n_brks = 8,
              xlab = NULL, ylab = NULL, 
              ttle.txt = "(b) Posterior standard deviation (log mm per day)" ) + theme_bw() +
    coord_cartesian(xlim = c(-125,-89)) + facet_grid( . ~ likelihood),
  nrow = 2
  # nrow = 1
)
dev.off()

#================================================
# Application 3: DJF 20-year return values
#================================================
# Load data 
tasRV20_DF_all <- read.csv("data/C20C_DJFtasRV20_trend.csv")
tasRV20_DF_all <- tasRV20_DF_all[abs(tasRV20_DF_all$latitude) < 85,]

# Remove NA's for fitting
tasRV20_DF <- tasRV20_DF_all[!is.na(tasRV20_DF_all$rv20),]

# Design matrix
Xmat <- unname(lm(rv20 ~ as.factor(region), x = TRUE, data = tasRV20_DF)$x)

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
# Constants for NNGP 
constants <- list( 
  nu = 5, k = 15, mu_HP1 = 100, X_sigma = Xmat, sigma_HP1 = 10,
  X_Sigma = Xmat, Sigma_HP1 = 10, maxAnisoDist = 20 # maxDist = 22.07 km*1000
)
# MCMC using the NNGP likelihood
Rmodel <- nsgpModel(likelihood = "NNGP", constants = constants, coords = round(xyz.crds, 4), data = tasRV20_DF$rv20, 
                    tau_model = "constant", sigma_model = "logLinReg", mu_model = "constant", Sigma_model = "compRegIso")
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
samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 0)
# Prediction
Xmat_pred <- unname(lm(rnorm(nrow(tasRV20_DF_all)) ~ region, x = TRUE, data = tasRV20_DF_all)$x)
# Convert lon/lat to x/y/z
xyz.crds <- matrix(NA,nrow(tasRV20_DF_all),3)
# Transform degrees to radians
lat.radians <- tasRV20_DF_all$latitude*(pi/180)
lon.radians <- tasRV20_DF_all$longitude*(pi/180)
for(i in 1:nrow(xyz.crds)){
  xyz.crds[i,1] <- 6.371*cos(lat.radians[i])*cos(lon.radians[i]) # Earth diameter ~ 6371km
  xyz.crds[i,2] <- 6.371*cos(lat.radians[i])*sin(lon.radians[i])
  xyz.crds[i,3] <- 6.371*sin(lat.radians[i])
}
pred_NNGP <- nsgpPredict(model = Rmodel, samples = samples[seq(from=5005,to=10000,by=5),],
                         coords.predict = round(xyz.crds, 4), PX_sigma = Xmat_pred, PX_Sigma = Xmat_pred)

# PLOTS for Application 3
col.lim <- c(242,308)
n_brks = 7; brk_rd = 0
col.pal <- brewer.pal(9, "Spectral")[9:1]
col.nm <- "K"
brks = round(seq(from = col.lim[1], to = col.lim[2], length = n_brks), brk_rd)
color <- tasRV20_DF_all$rv20
color[!is.na(color) & color < col.lim[1]] <- col.lim[1]
color[!is.na(color) & color > col.lim[2]] <- col.lim[2]
sc <- scale_fill_gradientn( colors = col.pal, name = col.nm, limits = col.lim, 
                            breaks = seq(245,305,length = 7) )
gd <- guides(fill = guide_colorbar( barheight = 15 ))

#### FIGURE 8
pdf("~/CASCADE_Projects/bayes_nsgp/paper_draft/tasRV20.pdf", height = 4.5, width = 8.5)
ggplot(tasRV20_DF_all, aes(x = longitude, y = latitude, fill = color)) +
  geom_tile() + coord_fixed(ratio = 1) +
  geom_polygon( data = map_data("world"), aes(x = long, y = lat, group = group),
                color = "black", fill = NA, size = 0.2) + sc + gd + 
  xlab(NULL) + ylab(NULL) + theme_bw() +
  geom_hline(yintercept = 0, size = 0.25) + geom_vline(xintercept = c(-120, 0, 120), size = 0.25) +
  annotate(geom= "text", label = paste0("Region ", c(5, 4, 6, 5, 2, 1, 3, 2)), 
           x = c(-150, -60, 60, 150, -150, -60, 60, 150 ), 
           y = c(rep(90,4), rep(-90,4)), size = 4)
dev.off()

postSamp <- samples[3001:5000,]
f_UB <- function(x){quantile(x, 0.975)}
f_LB <- function(x){quantile(x, 0.025)}

Sigma_post <- data.frame(
  Region = paste0("Region ", 1:6),
  Type = "(a) Spatial length scale (km)",
  Mean = 1000*c(mean(exp(postSamp[,1])), mean(exp(postSamp[,1] + postSamp[,2])),
                mean(exp(postSamp[,1] + postSamp[,3])), mean(exp(postSamp[,1] + postSamp[,4])),
                mean(exp(postSamp[,1] + postSamp[,5])), mean(exp(postSamp[,1] + postSamp[,6]))),
  UB = 1000*c(f_UB(exp(postSamp[,1])), f_UB(exp(postSamp[,1] + postSamp[,2])),
              f_UB(exp(postSamp[,1] + postSamp[,3])), f_UB(exp(postSamp[,1] + postSamp[,4])),
              f_UB(exp(postSamp[,1] + postSamp[,5])), f_UB(exp(postSamp[,1] + postSamp[,6]))),
  LB = 1000*c(f_LB(exp(postSamp[,1])), f_LB(exp(postSamp[,1] + postSamp[,2])),
              f_LB(exp(postSamp[,1] + postSamp[,3])), f_LB(exp(postSamp[,1] + postSamp[,4])),
              f_LB(exp(postSamp[,1] + postSamp[,5])), f_LB(exp(postSamp[,1] + postSamp[,6]))) )
sigma_post <- data.frame(
  Region = paste0("Region ", 1:6),
  Type = "(b) Spatial standard deviation (Kelvin)",
  Mean = c(mean(exp(postSamp[,7])), mean(exp(postSamp[,7] + postSamp[,8])),
           mean(exp(postSamp[,7] + postSamp[,9])), mean(exp(postSamp[,7] + postSamp[,10])),
           mean(exp(postSamp[,7] + postSamp[,11])), mean(exp(postSamp[,7] + postSamp[,12]))),
  UB = c(f_UB(exp(postSamp[,7])), f_UB(exp(postSamp[,7] + postSamp[,8])),
         f_UB(exp(postSamp[,7] + postSamp[,9])), f_UB(exp(postSamp[,7] + postSamp[,10])),
         f_UB(exp(postSamp[,7] + postSamp[,11])), f_UB(exp(postSamp[,7] + postSamp[,12]))),
  LB = c(f_LB(exp(postSamp[,7])), f_LB(exp(postSamp[,7] + postSamp[,8])),
         f_LB(exp(postSamp[,7] + postSamp[,9])), f_LB(exp(postSamp[,7] + postSamp[,10])),
         f_LB(exp(postSamp[,7] + postSamp[,11])), f_LB(exp(postSamp[,7] + postSamp[,12]))) )

#### FIGURE 9
pdf("app3_posteriors.pdf", height = 5, width = 8.5)
grid.arrange(
  ggplot(Sigma_post, aes(x = Region, y = Mean, ymin = LB, ymax = UB)) +
    geom_pointrange(size = 0.5) + ylab("Posterior mean with 95% CrI") + xlab(NULL) + facet_grid(~Type) +
    theme(strip.text.x = element_text(size=15)),
  ggplot(sigma_post, aes(x = Region, y = Mean, ymin = LB, ymax = UB)) +
    geom_pointrange(size = 0.5) + ylab("Posterior mean with 95% CrI") + xlab(NULL) + facet_grid(~Type) +
    theme(strip.text.x = element_text(size=13)),
  ncol = 2
)
dev.off()

postMean <- colMeans(pred_NNGP$pred)
postMean[postMean < col.lim[1]] <- col.lim[1]
postMean[postMean > col.lim[2]] <- col.lim[2]
postSD <- apply(pred_NNGP$pred, 2, sd)
SDlim <- c(0.09,1.25)
postSD[postSD < SDlim[1]] <- SDlim[1]
postSD[postSD > SDlim[2]] <- SDlim[2]

g1 <-   ggplot(tasRV20_DF_all, aes(x = longitude, y = latitude, fill = postMean)) +
  geom_tile() + coord_fixed(ratio = 1) +
  geom_polygon( data = map_data("world"), aes(x = long, y = lat, group = group),
                color = "black", fill = NA, size = 0.2) + sc + gd + 
  xlab(NULL) + ylab(NULL) + theme_bw() + ggtitle("(a) Posterior predictive mean")
g2 <- ggplot(tasRV20_DF_all, aes(x = longitude, y = latitude, fill = postSD)) +
  geom_tile() + coord_fixed(ratio = 1) +
  geom_polygon( data = map_data("world"), aes(x = long, y = lat, group = group),
                color = "black", fill = NA, size = 0.2) + 
  scale_fill_gradientn( colors = viridis_pal(option = "magma")(9), name = "K", limits = SDlim, 
                        trans = "log", breaks = c(0.1, 0.2, 0.3, 0.6, 1, 1.7)) + 
  gd + xlab(NULL) + ylab(NULL) + theme_bw() + ggtitle("(b) Posterior predictive standard deviation")

#### FIGURE 10
pdf("~/CASCADE_Projects/bayes_nsgp/paper_draft/app3_preds.pdf", height = 7.75, width = 8)
grid.arrange(g1, g2, ncol = 1)
dev.off()

#================================================
# Comparing computational times
#================================================
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
    Rmodel_fullGP <- nsgpModel( likelihood = "fullGP", sigma_model = "logLinReg", Sigma_model = "compRegIso",
                                mu_model = "linReg", tau_model = "approxGP", 
                                constants = constants, coords = coords, data = rnorm(N) ) #, returnModelComponents = TRUE )
    fullGPtime[n] <- system.time(Rmodel_fullGP$calculate())[3]
    rm(Rmodel_fullGP)
  }
  
  # NNGP
  Rmodel_NNGP <- nsgpModel( likelihood = "NNGP", sigma_model = "logLinReg", Sigma_model = "compRegIso",
                            mu_model = "linReg", tau_model = "approxGP", 
                            constants = constants, coords = coords, data = rnorm(N) )
  NNGPtime[n] <- system.time(Rmodel_NNGP$calculate())[3]
  rm(Rmodel_NNGP)
  
  # SGV
  Rmodel_SGV <- nsgpModel( likelihood = "SGV", sigma_model = "logLinReg", Sigma_model = "compRegIso",
                           mu_model = "linReg", tau_model = "approxGP", 
                           constants = constants, coords = coords, data = rnorm(N) )
  SGVtime[n] <- system.time(Rmodel_SGV$calculate())[3]
  rm(Rmodel_SGV)
  
}
plotDF <- data.frame( N = rep(Nvec, 3), Time = c(fullGPtime, SGVtime, NNGPtime),
                      Likelihood = rep(c("Exact GP", "SGV", "NNGP"), each = length(Nvec)) )

# FIGURE 11
pdf("compTime.pdf", height = 3, width = 7)
ggplot(plotDF, aes(x = N, y = Time, color = Likelihood)) + geom_point() + geom_line() +
  scale_y_continuous(trans="log", breaks=10^(-2:3), name = "Time (s)") +
  scale_x_continuous(breaks = c(50,1000,2000,5000,10000), name = "Sample size (N)") + 
  scale_color_manual(values = brewer.pal(3,"Set1"))
dev.off()












