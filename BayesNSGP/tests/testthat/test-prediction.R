test_that("NNGP.pred implements exact simple kriging when all training points are used", {
  # With k = M, every training location is used as a neighbor, so the
  # predictive mean and variance should match ordinary (non-approximate)
  # simple kriging under the same exponential covariance.
  set.seed(4)
  M      <- 6
  rho    <- 0.4
  sigma2 <- 1
  coords <- matrix(runif(2 * M), ncol = 2)
  w      <- rnorm(M)
  x0     <- c(0.5, 0.5)

  n_draws <- 3000
  draws <- replicate(n_draws, NNGP.pred(x0, coords, rho, sigma2, w, k = M))

  # Independent reference computation using ordinary Euclidean distance.
  # NNGP.pred computes distances via sf::st_distance with no CRS set, which
  # for planar coordinates reduces to the same Euclidean metric as dist().
  d0 <- sqrt(rowSums((coords - matrix(x0, M, 2, byrow = TRUE))^2))
  D  <- as.matrix(dist(coords))
  C0 <- sigma2 * exp(-d0 / rho)
  C  <- sigma2 * exp(-D / rho) + diag(1e-6, M)

  m_true <- as.numeric(C0 %*% solve(C, w))
  v_true <- sigma2 - as.numeric(C0 %*% solve(C, C0))

  expect_lt(abs(mean(draws) - m_true), 0.1)
  expect_lt(abs(var(draws) - v_true), 0.3)
})
