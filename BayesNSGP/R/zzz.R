

# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

    packageStartupMessage("Loading BayesNSGP. Registering the following distribution: dmnorm_NN_GP\n")

    suppressMessages({
        registerDistributions(list(
            dmnorm_NN_GP = list(
                BUGSdist = 'dmnorm_NN_GP(mu, AD, neighbors.id)',
                types = c('value = double(1)', 'mu = double(1)', 'AD = double(2)', 'neighbors.id = double(2)'),
                mixedSizes = TRUE)
        ), verbose = FALSE)
    })
}
