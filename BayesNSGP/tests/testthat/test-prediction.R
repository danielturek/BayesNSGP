test_that("NNGP.pred implements exact simple kriging when all training points are used", {
  # NNGP.pred documents `coords_sorted` as "output from computeNeighbors", i.e.
  # already sorted by x-coordinate, with `w.x0.all` reordered to match -- the
  # function indexes its neighbor set by position in that sorted order. This
  # test therefore sorts the training data before calling it.
  #
  # k = M and x0 placed to the right of every training point together force
  # NNGP.pred's first branch to use all M training points as neighbors, so
  # the predictive mean and variance should match ordinary (non-approximate)
  # simple kriging under the same exponential covariance.
  set.seed(4)
  M      <- 6
  rho    <- 0.4
  sigma2 <- 1
  coords <- matrix(runif(2 * M), ncol = 2)
  w      <- rnorm(M)

  ord <- order(coords[, 1])
  coords_sorted <- coords[ord, ]
  w_sorted      <- w[ord]

  # Guaranteed to sort last, so all M training points become its neighbors.
  x0 <- c(max(coords[, 1]) + 1, 0.5)

  n_draws <- 3000
  draws <- replicate(
    n_draws,
    NNGP.pred(x0, coords_sorted, rho, sigma2, w_sorted, k = M)
  )

  # Independent reference computation using ordinary Euclidean distance.
  d0 <- sqrt(rowSums((coords_sorted - matrix(x0, M, 2, byrow = TRUE))^2))
  D  <- as.matrix(dist(coords_sorted))
  C0 <- sigma2 * exp(-d0 / rho)
  C  <- sigma2 * exp(-D / rho) + diag(1e-6, M)

  m_true <- as.numeric(C0 %*% solve(C, w_sorted))
  v_true <- sigma2 - as.numeric(C0 %*% solve(C, C0))

  expect_lt(abs(mean(draws) - m_true), 0.1)
  expect_lt(abs(var(draws) - v_true), 0.3)
})
