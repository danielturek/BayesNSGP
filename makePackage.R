

library(devtools)
library(roxygen2)

## use roxygen comments to generate man/*.Rd help files
setwd('~/github/BayesNSGP/BayesNSGP')
document()

## build the package tarball file
setwd('~/github/BayesNSGP')
system('R CMD BUILD BayesNSGP')

## use check() to "check" the package
setwd('~/github/BayesNSGP/BayesNSGP')
check('.')

## install the package from GitHub:
library(devtools)
install_github('danielturek/BayesNSGP', subdir = 'BayesNSGP')
library(BayesNSGP)



