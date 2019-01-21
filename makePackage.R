

library(devtools)
library(roxygen2)

## change this if you want:
baseDir <- '~/github/BayesNSGP/'

## use roxygen comments to generate man/*.Rd help files:
document(paste0(baseDir, 'BayesNSGP'))

## build the package tarball file:
system(paste0('R CMD BUILD ', baseDir, 'BayesNSGP'))

## check the package for passing CRAN tests and checks:
check(paste0(baseDir, 'BayesNSGP'))

## install the BayesNSGP package from GitHub:
library(devtools)
remove.packages('BayesNSGP')
install_github('danielturek/BayesNSGP', subdir = 'BayesNSGP')
library(BayesNSGP)



