test_that("dmnorm_NN_GP matches the full multivariate normal density when k = M - 1", {
  # When every location conditions on every earlier location, the NNGP
  # factorisation is exact rather than approximate, so it must reproduce the
  # ordinary multivariate normal log-density for the same covariance.
  set.seed(1)
  M      <- 8
  k      <- M - 1
  rho    <- 0.6
  sigma2 <- 1.3
  
  coords <- matrix(runif(2 * M), ncol = 2)
  nb <- computeNeighbors(coords, k = k, proj = 4326)
  # sf::st_distance() returns a 'units' object (metres) for a geographic CRS;
  # strip it here so plain arithmetic like exp() works below. computeAD()
  # itself doesn't need this -- NIMBLE's own type coercion handles it.
  edist <- matrix(as.numeric(nb$edist_sorted), M, M)
  
  AD <- computeAD(
    edist        = nb$edist_sorted,
    nid.dist     = nb$neighbors_dist,
    neighbors.id = nb$neighbor_idx,
    rho          = rho,
    sigma2       = sigma2,
    k            = k
  )
  
  # Built from the same distance matrix computeNeighbors returned, so this is
  # an exact comparison rather than one that depends on a second, independent
  # distance calculation.
  Sigma <- sigma2 * exp(-edist / rho) + diag(1e-6, M)
  
  x  <- rnorm(M)
  mu <- rep(0, M)
  
  ll_nngp <- dmnorm_NN_GP(x, mu, AD, nb$neighbor_idx, log = TRUE)
  
  R <- chol(Sigma)
  quad   <- sum(backsolve(R, x - mu, transpose = TRUE)^2)
  ll_ref <- -0.5 * (M * log(2 * pi) + 2 * sum(log(diag(R))) + quad)
  
  expect_equal(ll_nngp, ll_ref, tolerance = 1e-6)
})

test_that("rmnorm_NN_GP simulates with the correct covariance when k = M - 1", {
  set.seed(2)
  M      <- 6
  k      <- M - 1
  rho    <- 0.5
  sigma2 <- 1
  
  coords <- matrix(runif(2 * M), ncol = 2)
  nb <- computeNeighbors(coords, k = k, proj = 4326)
  edist <- matrix(as.numeric(nb$edist_sorted), M, M)
  AD <- computeAD(nb$edist_sorted, nb$neighbors_dist, nb$neighbor_idx, rho, sigma2, k)
  
  mu <- rep(0, M)
  R  <- 4000
  draws <- t(replicate(R, rmnorm_NN_GP(1, mu, AD, nb$neighbor_idx)))
  
  Sigma_true <- sigma2 * exp(-edist / rho) + diag(1e-6, M)
  Sigma_hat  <- cov(draws)
  
  # Monte Carlo comparison over many draws; a generous absolute tolerance
  # keeps the test from being flaky while still catching a wrong covariance.
  expect_true(max(abs(Sigma_hat - Sigma_true)) < 0.15)
})
