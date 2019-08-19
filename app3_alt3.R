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
tasRV20_DF_all <- tasRV20_DF_all[abs(tasRV20_DF_all$latitude) < 85,]
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

# # Subset (temporarily)
# tasRV20_DF <- tasRV20_DF[sample(nrow(tasRV20_DF), 1000),]

# Standardize
 # - 275)/20

# Design matrix
Xmat <- unname(lm(rv20 ~ region, x = TRUE, data = tasRV20_DF)$x)
# N <- nrow(tasRV20_DF)
z <- tasRV20_DF$rv20

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


# Constants for NNGP ============================
constants <- list( 
  nu = 5, k = 15, mu_HP1 = 100,
  X_sigma = Xmat, sigma_HP1 = 10,
  X_Sigma = Xmat, Sigma_HP1 = 10, 
  maxAnisoDist = 20 # maxDist = 22.07 km*1000
)

#================================================
# MCMC using the NNGP likelihood
#================================================
prt <- proc.time()
Rmodel <- nsgpModel(likelihood = "NNGP", constants = constants, 
                    coords = round(xyz.crds, 4), z = tasRV20_DF$rv20, 
                    tau_model = "constant", 
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
samples <- runMCMC(Cmcmc, niter = 5000, nburnin = 0)
time_mcmc <- proc.time() - prt
save(samples, time_mcmc, time_build, file = "app3_alt3_full.RData")
# 
# tmp <- samples
# 
# par(ask=TRUE)
# for(h in 1:ncol(samples)) plot(samples[4001:5000,h], type = "l", main = colnames(samples)[h])
# par(ask=FALSE)

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
predCoords <- round(xyz.crds, 4)

PREDconstants <- list(PX_sigma = Xmat_pred, PX_Sigma = Xmat_pred)

prt <- proc.time()
pred_NNGP <- nsgpPredict(model = Rmodel, samples = samples[seq(from=3002,to=5000,by=2),],
                         coords.predict = predCoords, constants = PREDconstants)
time_NNGP_pred <- proc.time() - prt
save(pred_NNGP, time_NNGP_pred, file = "app3_alt3_pred.RData")

# library(fields)
# quilt.plot(tasRV20_DF_all[,1:2], colMeans(pred_NNGP$pred))

