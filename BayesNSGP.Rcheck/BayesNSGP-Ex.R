pkgname <- "BayesNSGP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BayesNSGP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("calcQF")
### * calcQF

flush(stderr()); flush(stdout())

### Name: calcQF
### Title: Calculate the Gaussian quadratic form for the NNGP approximation
### Aliases: calcQF

### ** Examples

# TODO




cleanEx()
nameEx("calculateAD_ns")
### * calculateAD_ns

flush(stderr()); flush(stdout())

### Name: calculateAD_ns
### Title: Calculate A and D matrices for the NNGP approximation
### Aliases: calculateAD_ns

### ** Examples

# TODO




cleanEx()
nameEx("calculateU_ns")
### * calculateU_ns

flush(stderr()); flush(stdout())

### Name: calculateU_ns
### Title: Calculate the (sparse) matrix U
### Aliases: calculateU_ns

### ** Examples

# TODO




cleanEx()
nameEx("conditionLatentObs")
### * conditionLatentObs

flush(stderr()); flush(stdout())

### Name: conditionLatentObs
### Title: Assign conditioning sets for the SGV approximation
### Aliases: conditionLatentObs

### ** Examples

# TODO




cleanEx()
nameEx("determineNeighbors")
### * determineNeighbors

flush(stderr()); flush(stdout())

### Name: determineNeighbors
### Title: Determine the k-nearest neighbors for each spatial coordinate.
### Aliases: determineNeighbors

### ** Examples

# TODO




cleanEx()
nameEx("dmnorm_nngp")
### * dmnorm_nngp

flush(stderr()); flush(stdout())

### Name: dmnorm_nngp
### Title: Function for the evaluating the NNGP approximate density.
### Aliases: dmnorm_nngp

### ** Examples

# TODO




cleanEx()
nameEx("dmnorm_sgv")
### * dmnorm_sgv

flush(stderr()); flush(stdout())

### Name: dmnorm_sgv
### Title: Function for the evaluating the SGV approximate density.
### Aliases: dmnorm_sgv

### ** Examples

# TODO




cleanEx()
nameEx("inverseEigen")
### * inverseEigen

flush(stderr()); flush(stdout())

### Name: inverseEigen
### Title: Calculate covariance elements based on eigendecomposition
###   components
### Aliases: inverseEigen

### ** Examples

# TODO




cleanEx()
nameEx("matern_corr")
### * matern_corr

flush(stderr()); flush(stdout())

### Name: matern_corr
### Title: Calculate a stationary Matern correlation matrix
### Aliases: matern_corr

### ** Examples

# TODO




cleanEx()
nameEx("nsCorr")
### * nsCorr

flush(stderr()); flush(stdout())

### Name: nsCorr
### Title: Calculate a nonstationary Matern correlation matrix
### Aliases: nsCorr

### ** Examples

# TODO




cleanEx()
nameEx("nsCrosscorr")
### * nsCrosscorr

flush(stderr()); flush(stdout())

### Name: nsCrosscorr
### Title: Calculate a nonstationary Matern cross-correlation matrix
### Aliases: nsCrosscorr

### ** Examples

# TODO




cleanEx()
nameEx("nsCrossdist")
### * nsCrossdist

flush(stderr()); flush(stdout())

### Name: nsCrossdist
### Title: Calculate coordinate-specific cross-distance matrices
### Aliases: nsCrossdist

### ** Examples

# TODO




cleanEx()
nameEx("nsCrossdist3d")
### * nsCrossdist3d

flush(stderr()); flush(stdout())

### Name: nsCrossdist3d
### Title: Calculate coordinate-specific distance matrices, only for
###   nearest neighbors and store in an array
### Aliases: nsCrossdist3d

### ** Examples

# TODO




cleanEx()
nameEx("nsDist")
### * nsDist

flush(stderr()); flush(stdout())

### Name: nsDist
### Title: Calculate coordinate-specific distance matrices
### Aliases: nsDist

### ** Examples

# TODO




cleanEx()
nameEx("nsDist3d")
### * nsDist3d

flush(stderr()); flush(stdout())

### Name: nsDist3d
### Title: Calculate coordinate-specific distance matrices, only for
###   nearest neighbors and store in an array
### Aliases: nsDist3d

### ** Examples

# TODO




cleanEx()
nameEx("nsgpModel")
### * nsgpModel

flush(stderr()); flush(stdout())

### Name: nsgpModel
### Title: NIMBLE code for a generic nonstationary GP model
### Aliases: nsgpModel

### ** Examples

# TODO




cleanEx()
nameEx("orderCoordinatesMMD")
### * orderCoordinatesMMD

flush(stderr()); flush(stdout())

### Name: orderCoordinatesMMD
### Title: Order coordinates according to a maximum-minimum distance
###   criterion.
### Aliases: orderCoordinatesMMD

### ** Examples

# TODO




cleanEx()
nameEx("sgvSetup")
### * sgvSetup

flush(stderr()); flush(stdout())

### Name: sgvSetup
### Title: One-time setup wrapper function for the SGV approximation
### Aliases: sgvSetup

### ** Examples

# TODO




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
