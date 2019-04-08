
1

library(devtools); library(roxygen2)
baseDir <- '~/github/BayesNSGP/'
(tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE))
for(file in tarFiles) system(paste0('rm ', file))
system(paste0('rm -f ', paste0(baseDir, 'BayesNSGP/NAMESPACE')))
document(paste0(baseDir, 'BayesNSGP'))
namespaceFilename <- paste0(baseDir, 'BayesNSGP/NAMESPACE')
namespace <- readLines(namespaceFilename)
if(length(namespace) >= 2) namespace <- namespace[2:length(namespace)]
namespace <- c('import(methods)', namespace)
namespace <- c('import(stats)',   namespace)
writeLines(namespace, namespaceFilename)
system(paste0('R CMD BUILD ', baseDir, 'BayesNSGP'))
check(paste0(baseDir, 'BayesNSGP'))
remove.packages('BayesNSGP')
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
(lastTarFile <- tarFiles[length(tarFiles)])
system(paste0('R CMD install ', lastTarFile))
qqq
1
library(BayesNSGP)
matern_corr
calcQF
calculateAD_ns
dmnorm_nngp
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



