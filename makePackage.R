

library(devtools)
library(roxygen2)
library(nimble)

baseDir <- '~/github/BayesNSGP/'
if(!('makePackage.R' %in% list.files(baseDir))) stop('change baseDir directory')
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
for(file in tarFiles) system(paste0('rm ', file))

if(TRUE) {
    ## remake NAMESPACE
    system(paste0('rm -f ', paste0(baseDir, 'BayesNSGP/NAMESPACE')))
    document(paste0(baseDir, 'BayesNSGP'))
    namespaceFilename <- paste0(baseDir, 'BayesNSGP/NAMESPACE')
    namespace <- readLines(namespaceFilename)
    if(length(namespace) >= 2) namespace <- namespace[2:length(namespace)]
    namespace <- c('importFrom("sf", "st_as_sf", "st_distance")', namespace)
    namespace <- c('import(ggplot2)', namespace)
    namespace <- c('import(StatMatch)', namespace)
    namespace <- c('importFrom("Matrix", "sparseMatrix")', namespace)
    namespace <- c('import(FNN)', namespace)
    namespace <- c('import(nimble)', namespace)
    namespace <- c('import(methods)', namespace)
    namespace <- c('importFrom("stats", "dist", "rnorm")', namespace)
    writeLines(namespace, namespaceFilename)
}

devtools::build(paste0(baseDir, 'BayesNSGP'))

check(paste0(baseDir, 'BayesNSGP'))

try(remove.packages('BayesNSGP'), silent = TRUE)
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
(lastTarFile <- tarFiles[length(tarFiles)])
system(paste0('/usr/local/bin/R CMD install ', lastTarFile, ' --build-vignettes'))

devtools::install('.', build_vignettes = TRUE)

q('no')    ## quit R

1          ## restart R

library(BayesNSGP)

browseVignettes('BayesNSGP')

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



