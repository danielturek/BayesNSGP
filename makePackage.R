

## Mark:
## use this code below to build and re-install the package,
## when you make changes to the source files.

library(devtools)
library(roxygen2)
library(nimble, warn.conflicts = FALSE)
if('GitHub' %in% list.files('~/Documents')) {   ## this should work for either of us
    baseDir <- '~/Documents/Github/BayesNSGP/'  ## Risser
} else { baseDir <- '~/github/BayesNSGP/' }     ## Turek
if(!('makePackage.R' %in% list.files(baseDir))) stop('change baseDir directory')
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
for(file in tarFiles) system(paste0('rm ', file))
system(paste0('rm -f ', paste0(baseDir, 'BayesNSGP/NAMESPACE')))
document(paste0(baseDir, 'BayesNSGP'))
namespaceFilename <- paste0(baseDir, 'BayesNSGP/NAMESPACE')
namespace <- readLines(namespaceFilename)
if(length(namespace) >= 2) namespace <- namespace[2:length(namespace)]
namespace <- c('import(StatMatch)', namespace)
namespace <- c('import(Matrix)', namespace)
namespace <- c('import(FNN)', namespace)
namespace <- c('import(nimble)', namespace)
namespace <- c('import(methods)', namespace)
namespace <- c('importFrom("stats", "dist", "rnorm")', namespace)
writeLines(namespace, namespaceFilename)
system(paste0('R CMD BUILD ', baseDir, 'BayesNSGP'))
##check(paste0(baseDir, 'BayesNSGP'))
try(remove.packages('BayesNSGP'), silent = TRUE)
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
lastTarFile <- tarFiles[length(tarFiles)]
system(paste0('R CMD install ', lastTarFile))

## now quit R

library(BayesNSGP)

##
## stop here
##



##matern_corr
##calcQF
##calculateAD_ns
##dmnorm_nngp
##C_calcQF <- nimble::compileNimble(calcQF)



library(nimble)
nsCorrC <- compileNimble(nsCorr)
calculateAD_nsC <- compileNimble(calculateAD_ns)   ## FAILING



as.name(as.character('funName'))
as.name(as.character('package::funName'))

class(quote(package::fun(arg1, arg2)))



## install the BayesNSGP package from GitHub:
library(devtools)
remove.packages('BayesNSGP')
install_github('danielturek/BayesNSGP', subdir = 'BayesNSGP')
library(BayesNSGP)



