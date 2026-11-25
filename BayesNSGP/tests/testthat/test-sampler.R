test_that("sampler_RW_NN_GP reproduces the generic RW sampler within an MCMC run", {
  # The local-update sampler is an algebraic shortcut: it evaluates the same
  # Metropolis-Hastings ratio as NIMBLE's scalar RW sampler, but from the
  # node's own conditional density and those of its reverse neighbours instead
  # of the full NNGP density. Both samplers draw rnorm(1) for the proposal and
  # then call decide(), and both use the same scale adaptation. So with the same
  # seed, initial values and sampler order the two chains must coincide, up to
  # floating-point rounding from the different arithmetic path.
  skip_on_cran()
  skip_if_not_installed("nimble")
  
  set.seed(1)
  M <- 30
  k <- 5
  
  # --- neighbour structure ---------------------------------------------------
  coords <- matrix(runif(2 * M), ncol = 2)
  nb <- computeNeighbors(coords, k = k, proj = 4326)
  
  # sf returns 'units' objects (metres); strip them and rescale so distances
  # lie in (0, 1] and the range prior below is sensible.
  dmax_raw <- max(as.numeric(nb$edist_sorted))
  edist    <- matrix(as.numeric(nb$edist_sorted), M, M) / dmax_raw
  nid.dist <- matrix(as.numeric(nb$neighbors_dist),
                     nrow(nb$neighbors_dist), ncol(nb$neighbors_dist)) / dmax_raw
  n.id     <- nb$neighbor_idx
  
  # --- simulate a small Poisson dataset --------------------------------------
  rho0    <- 0.3
  sigma20 <- 1
  AD0   <- computeAD(edist, nid.dist, n.id, rho0, sigma20, k)
  omega0 <- rmnorm_NN_GP(n = 1, mu = rep(0, M), AD = AD0, neighbors.id = n.id)
  y     <- rpois(M, exp(1 + omega0))
  
  constants <- list(
    N = M, k = k, mu = rep(0, M), mu0 = 1,
    edist = edist, nid.dist = nid.dist, n.id = n.id,
    dmin = min(edist[edist != 0]), dmax = max(edist)
  )
  
  code <- nimble::nimbleCode({
    sigma2 ~ dinvgamma(3, 1)
    rho    ~ dunif(dmin / 3, dmax / 3)
    AD[1:N, 1:(k + 1)] <- computeAD(edist[1:N, 1:N], nid.dist[1:N, 1:(k + 1)],
                                    n.id[1:N, 1:k], rho, sigma2, k)
    omega[1:N] ~ dmnorm_NN_GP(mu[1:N], AD[1:N, 1:(k + 1)], n.id[1:N, 1:k])
    for (i in 1:N) {
      log(lambda[i]) <- mu0 + omega[i]
      y[i] ~ dpois(lambda[i])
    }
  })
  
  inits <- list(omega = omega0, rho = 0.2, sigma2 = 1)
  
  run_chain <- function(sampler_type) {
    model <- nimble::nimbleModel(code, constants = constants,
                                 data = list(y = y), inits = inits)
    conf <- nimble::configureMCMC(model, monitors = c("rho", "sigma2", "omega"))
    conf$removeSamplers("omega")
    for (i in 1:M) {
      if (sampler_type == "RW_NN_GP") {
        conf$addSampler(target = paste0("omega[", i, "]"), type = "RW_NN_GP",
                        control = list(AD = "AD", neighbors.id = n.id))
      } else {
        conf$addSampler(target = paste0("omega[", i, "]"), type = "RW")
      }
    }
    mcmc  <- nimble::buildMCMC(conf)
    cmod  <- nimble::compileNimble(model)
    cmcmc <- nimble::compileNimble(mcmc, project = model)
    nimble::runMCMC(cmcmc, niter = 400, nburnin = 0, nchains = 1,
                    inits = inits, setSeed = 123, progressBar = FALSE)
  }
  
  s_rw   <- run_chain("RW")
  s_nngp <- run_chain("RW_NN_GP")
  
  expect_equal(dim(s_nngp), dim(s_rw))
  expect_equal(colnames(s_nngp), colnames(s_rw))
  
  # The chain must actually be moving, otherwise the comparison is vacuous.
  om <- s_nngp[, grep("^omega\\[", colnames(s_nngp))]
  expect_true(mean(apply(om, 2, function(x) length(unique(x)) > 1)) > 0.5)
  expect_true(length(unique(s_nngp[, "rho"])) > 1)
  
  # Same chain up to rounding. (expect_identical() is too strict here: the two
  # samplers sum the same terms in a different order.)
  expect_equal(s_nngp, s_rw, tolerance = 1e-8)
})
