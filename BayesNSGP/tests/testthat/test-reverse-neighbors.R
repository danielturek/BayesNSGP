test_that("get_single_reverse_neighbors matches its own documented example", {
  # Taken directly from the function's @examples: node 2 conditions on node 1
  # (row 3), and node 3 conditions on node 1 (row 4). Both are reverse
  # neighbors of node 1.
  nn_matrix <- matrix(c(0, 0,
                         1, 0,
                         1, 2,
                         1, 3), nrow = 4, byrow = TRUE)

  expect_equal(
    get_single_reverse_neighbors(target_node = 1, neighbor_idx = nn_matrix),
    c(2, 3, 4)
  )

  # The last row cannot be a neighbor of anything later, so it has none.
  expect_equal(
    get_single_reverse_neighbors(target_node = 4, neighbor_idx = nn_matrix),
    integer(0)
  )
})

test_that("RWNNGP_setup resolves the correct local graph on a hand-worked example", {
  # A small, fully hand-verifiable neighbor structure (M = 4, k = 2):
  #   location 1: no neighbors
  #   location 2: neighbor 1
  #   location 3: neighbors 1, 2
  #   location 4: neighbors 2, 3
  neighbors.id <- matrix(c(0, 0,
                            1, 0,
                            1, 2,
                            2, 3), nrow = 4, byrow = TRUE)
  N.neighbors <- apply(neighbors.id, 1, function(x) sum(x != 0))  # c(0, 1, 2, 2)

  node_id <- 2
  Rneighbors.id <- get_single_reverse_neighbors(node_id, neighbors.id)  # c(3, 4)

  out <- RWNNGP_setup(
    node_id       = node_id,
    AD            = "AD",
    neighbors.id  = neighbors.id,
    Rneighbors.id = Rneighbors.id,
    N.neighbors   = N.neighbors,
    k             = 2
  )

  # Node 2's only forward neighbor is node 1, in column 1 of its row.
  expect_equal(out$Fneighbors.id, 1)
  expect_equal(out$AFnodes, "AD[2,1]")

  # Node 2 is neighbor 2 of node 3, and neighbor 1 of node 4.
  expect_equal(out$ARnodes, c("AD[3,2]", "AD[4,1]"))

  # The nodes whose conditional densities change when node 2 is updated.
  expect_equal(out$update_id, c(2, 3, 4))

  # Every A-coefficient feeding into the local log-likelihood recalculation.
  expect_equal(
    out$A.neighbors,
    c("AD[2,1]", "AD[3,1]", "AD[3,2]", "AD[4,1]", "AD[4,2]")
  )
})
