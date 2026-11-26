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


# ---------------------------------------------------------------------------
# Reference implementation: the straightforward version of Algorithm 2 that
# builds the full (M + 1) x (M + 1) distance matrix, as NNGP.pred() did before
# it was rewritten to search neighbors without forming that matrix.
# `dfun` maps a coordinate matrix to its distance matrix.
nngp_pred_reference <- function(x0, coords_sorted, rho, sigma2, w, k, dfun) {
  M    <- nrow(coords_sorted)
  allc <- rbind(coords_sorted, x0)
  o    <- order(allc[, 1])
  loc  <- which(o == M + 1)
  cs   <- allc[o, ]
  D    <- dfun(cs)
  if (loc <= k + 1) {
    idx <- seq_len(loc - 1)
  } else {
    ord <- order(D[loc, ]); ord <- ord[ord < loc]
    idx <- ord[1:k]
  }
  c0 <- sigma2 * exp(-D[loc, idx] / rho)
  C  <- sigma2 * exp(-matrix(D[idx, idx], length(idx)) / rho) + diag(1e-6, length(idx))
  list(m = sum(c0 * solve(C, w[idx])), v = sigma2 - sum(c0 * solve(C, c0)))
}

test_that("NNGP.pred matches the full-distance-matrix reference (Euclidean)", {
  set.seed(11)
  M <- 40; k <- 5
  coords <- matrix(runif(2 * M), ncol = 2)
  ord <- order(coords[, 1]); coords <- coords[ord, ]
  w <- rnorm(M)
  dfun <- function(cs) as.matrix(dist(cs))
  
  # one location in the middle (k nearest of many) and one near the left edge
  # (fewer than k + 1 points precede it, so all of them are neighbors)
  for (x0 in list(c(0.5, 0.5), c(min(coords[, 1]) + 0.01, 0.3))) {
    ref <- nngp_pred_reference(x0, coords, 0.3, 1.2, w, k, dfun)
    set.seed(1); a <- rnorm(1)
    set.seed(1); got <- NNGP.pred(x0, coords, 0.3, 1.2, w, k)
    expect_equal(got, ref$m + sqrt(ref$v) * a, tolerance = 1e-6)
  }
})

test_that("NNGP.pred handles many posterior draws and equals repeated single calls", {
  set.seed(12)
  M <- 30; k <- 4; nd <- 25
  coords <- matrix(runif(2 * M), ncol = 2)
  coords <- coords[order(coords[, 1]), ]
  W      <- matrix(rnorm(nd * M), nd, M)
  rho    <- runif(nd, 0.1, 0.5)
  sigma2 <- runif(nd, 0.5, 2)
  x0     <- c(0.6, 0.4)
  
  set.seed(3)
  vec <- NNGP.pred(x0, coords, rho, sigma2, W, k)
  set.seed(3)
  one <- vapply(seq_len(nd), function(d) NNGP.pred(x0, coords, rho[d], sigma2[d], W[d, ], k), numeric(1))
  # Same random numbers are consumed in the same order only if drawn one at a
  # time, so compare through a fixed reference of the Gaussian noise instead.
  set.seed(3); a <- rnorm(nd)
  ref <- vapply(seq_len(nd), function(d) {
    r <- nngp_pred_reference(x0, coords, rho[d], sigma2[d], W[d, ], k, function(cs) as.matrix(dist(cs)))
    r$m + sqrt(r$v) * a[d]
  }, numeric(1))
  
  expect_length(vec, nd)
  expect_equal(vec, ref, tolerance = 1e-6)
  expect_length(one, nd)   # the single-draw interface still returns a scalar each time
})

test_that("NNGP.pred uses proj and dist.scale consistently with computeNeighbors", {
  set.seed(13)
  M <- 30; k <- 5
  coords <- cbind(runif(M, -80, -70), runif(M, 35, 45))     # lon / lat
  coords <- coords[order(coords[, 1]), ]
  w <- rnorm(M)
  x0 <- c(-75, 40)
  # distances in km on the ellipsoid / sphere, the way sf computes them for lon/lat
  dfun <- function(cs) {
    p <- sf::st_as_sf(data.frame(x = cs[, 1], y = cs[, 2]), coords = c("x", "y"), crs = 4326)
    matrix(as.numeric(sf::st_distance(p)), nrow(cs)) / 1000
  }
  rho <- 150   # km
  ref <- nngp_pred_reference(x0, coords, rho, 1, w, k, dfun)
  set.seed(5); a <- rnorm(1)
  set.seed(5)
  got <- NNGP.pred(x0, coords, rho, 1, w, k, proj = 4326, dist.scale = 1000)
  expect_equal(got, ref$m + sqrt(ref$v) * a, tolerance = 1e-6)
  
  # ignoring the CRS measures distance in degrees, a different scale for rho
  set.seed(5)
  deg <- NNGP.pred(x0, coords, rho, 1, w, k)
  expect_false(isTRUE(all.equal(got, deg, tolerance = 1e-3)))
})

test_that("NNGP.pred draws from N(0, sigma2) when no training point precedes x0", {
  set.seed(14)
  M <- 10
  coords <- matrix(runif(2 * M, 1, 2), ncol = 2)
  coords <- coords[order(coords[, 1]), ]
  nd <- 4000
  draws <- NNGP.pred(c(0, 1.5), coords, rho = rep(0.3, nd), sigma2 = rep(2, nd),
                     w.x0.all = matrix(rnorm(nd * M), nd, M), k = 3)
  expect_lt(abs(mean(draws)), 0.15)
  expect_lt(abs(var(draws) - 2), 0.2)
})