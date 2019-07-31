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
library(splines)

# Load data =====================================
tasRV20_DF_all <- read.csv("data/C20C_DJFtasRV20_trend.csv")
tasRV20_DF_all <- tasRV20_DF_all[abs(tasRV20_DF_all$latitude) < 90,]
latShift <- mean(tasRV20_DF_all$latitude)
latScale <- sd(tasRV20_DF_all$latitude)
tasRV20_DF_all$Zlatitude <- (tasRV20_DF_all$latitude - latShift)/latScale
tasRV20_DF_all$rv20[tasRV20_DF_all$rv20 > 310] <- NA
tasRV20_DF <- tasRV20_DF_all[!is.na(tasRV20_DF_all$rv20),]

# Subset (temporarily)
# tasRV20_DF <- tasRV20_DF[sample(nrow(tasRV20_DF), 5000),]

z <- tasRV20_DF$rv20

# # Determine mean covariates
# plot(tasRV20_DF$latitude[tasRV20_DF$ind_land == 0], tasRV20_DF$rv20[tasRV20_DF$ind_land == 0],
#      pch = ".", xlab = "lat", ylab = "rv20", col = 2)
# lines(tasRV20_DF$latitude[tasRV20_DF$ind_land == 0 & !is.na(tasRV20_DF$rv20)],
#       predict(lm(rv20 ~ bs(Zlatitude,df = 6), data = tasRV20_DF[tasRV20_DF$ind_land == 0,])),
#       col = 2, lwd = 2)
# points(tasRV20_DF$latitude[tasRV20_DF$ind_land == 1], tasRV20_DF$rv20[tasRV20_DF$ind_land == 1],
#      pch = ".", xlab = "lat", ylab = "rv20", col = 4)
# lines(tasRV20_DF$latitude[tasRV20_DF$ind_land == 1 & !is.na(tasRV20_DF$rv20)],
#       predict(lm(rv20 ~ bs(Zlatitude,df = 6), data = tasRV20_DF[tasRV20_DF$ind_land == 1,])),
#       col = 4, lwd = 2)
# plot(tasRV20_DF$longitude[tasRV20_DF$ind_land == 0], tasRV20_DF$rv20[tasRV20_DF$ind_land == 0],
#      pch = ".", xlab = "lat", ylab = "rv20", col = 2)
# lines(tasRV20_DF$longitude[tasRV20_DF$ind_land == 0 & !is.na(tasRV20_DF$rv20)],
#       predict(lm(rv20 ~ bs(longitude,df = 6), data = tasRV20_DF[tasRV20_DF$ind_land == 0,])),
#       col = 2, lwd = 2)
# points(tasRV20_DF$longitude[tasRV20_DF$ind_land == 1], tasRV20_DF$rv20[tasRV20_DF$ind_land == 1],
#        pch = ".", xlab = "lat", ylab = "rv20", col = 4)
# lines(tasRV20_DF$longitude[tasRV20_DF$ind_land == 1 & !is.na(tasRV20_DF$rv20)],
#       predict(lm(rv20 ~ bs(longitude,df = 6), data = tasRV20_DF[tasRV20_DF$ind_land == 1,])),
#       col = 4, lwd = 2)

# For the mean, use basis splines with 6 degrees of freedom.
# The variance/spatial extent also represented using basis splines,
# but with only 4 degrees of freedom (and no land effect).
Xmat1 <- unname(lm(rv20 ~ bs(Zlatitude,df = 4)*ind_land, x = TRUE, 
                   data = tasRV20_DF)$x)
Xmat2 <- unname(lm(rv20 ~ bs(Zlatitude,df = 4), x = TRUE, 
                   data = tasRV20_DF)$x)
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

# Constants for NNGP ============================
constants <- list( 
  nu = 0.5, k = 15, 
  X_sigma = Xmat2, sigma_HP1 = 10,
  X_Sigma = Xmat2, Sigma_HP1 = 10, 
  maxAnisoDist = 20, # maxDist = 22.07 km*1000
  X_mu = Xmat1, mu_HP1 = 10 )

#================================================
# MCMC using the NNGP likelihood
#================================================
prt <- proc.time()
Rmodel <- nsgpModel(likelihood = "NNGP", constants = constants, 
                    coords = coords, z = z, tau_model = "constant", 
                    sigma_model = "logLinReg", mu_model = "linReg", 
                    Sigma_model = "compRegIso")
conf <- configureMCMC(Rmodel)
conf$removeSamplers(c("beta[1:10]"))
conf$addSampler(target = c("beta[1:10]"), type = "RW_block")
conf$removeSamplers(c("alpha[1:5]"))
conf$addSampler(target = c("alpha[1:5]"), type = "RW_block")
conf$removeSamplers(c("Sigma_coef1[1:5]"))
conf$addSampler(target = c("Sigma_coef1[1:5]"), type = "RW_block")
conf$getSamplers()
# Build
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
time_build <- proc.time() - prt
# Run
prt <- proc.time()
samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 0)
time_mcmc <- proc.time() - prt
save(samples, time_mcmc, time_build, file = "app3_rv20_NNGP_samples_time.RData")

# par(ask=TRUE)
# for(h in 1:ncol(samples)) plot(samples[,h], type = "l", main = colnames(samples)[h])
# par(ask=FALSE)
# 
# 
# summary(lm(rv20 ~ Zlatitude*ind_land, x = TRUE, 
#    data = tasRV20_DF))#$coefficients
# colMeans(samples[,9:12])

# Prediction
Xmat1_pred <- unname(lm(rnorm(nrow(tasRV20_DF_all)) ~ bs(Zlatitude,df = 4)*ind_land, x = TRUE, data = tasRV20_DF_all)$x)
Xmat2_pred <- unname(lm(rnorm(nrow(tasRV20_DF_all)) ~ bs(Zlatitude,df = 4), x = TRUE, data = tasRV20_DF_all)$x)
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
predCoords <- round(xyz.crds, 4)

prt <- proc.time()
pred_NNGP <- nsgpPredict(model = Rmodel, samples = samples[seq(from=10002,to=20000,by=2),], 
                         coords.predict = predCoords,
                         PX_sigma = Xmat2_pred, PX_Sigma = Xmat2_pred, PX_mu = Xmat1_pred)
time_NNGP_pred <- proc.time() - prt
save(pred_NNGP, time_NNGP_pred, file = "app3_rv20_NNGP_pred.RData")

