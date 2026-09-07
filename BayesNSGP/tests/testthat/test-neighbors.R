test_that("computeNeighbors returns a well-formed Vecchia-type ordering", {
  set.seed(3)
  M <- 30
  k <- 5
  coords <- matrix(runif(2 * M), ncol = 2)

  nb <- computeNeighbors(coords, k = k, proj = 4326)

  expect_equal(dim(nb$neighbor_idx),  c(M, k))
  expect_equal(dim(nb$coords_sorted), c(M, 2))
  expect_equal(dim(nb$edist_sorted),  c(M, M))

  # Locations are ordered by x-coordinate.
  expect_equal(nb$coords_sorted[, 1], sort(nb$coords_sorted[, 1]))

  # The first location precedes everything else, so it has no neighbors.
  expect_true(all(nb$neighbor_idx[1, ] == 0))

  # From the (k+2)-th location on, every neighbor slot should be filled, and
  # every neighbor index must refer to an earlier location in the sorted
  # order -- the defining property of a Vecchia-type ordering that makes the
  # NNGP a valid, triangular factorisation.
  for (i in (k + 2):M) {
    idx <- nb$neighbor_idx[i, ]
    expect_true(all(idx > 0))
    expect_true(all(idx < i))
    expect_length(unique(idx), k)
  }
})
