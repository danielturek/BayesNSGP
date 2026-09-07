
#================================================================================================
# Functions for computing and plotting the neighbors of the Nearest-Neighbor Gaussian Process (NNGP)
# Author: Fabian R. Ketwaroo
# Affiliation: Swiss Ornithological Institute
#================================================================================================


#================================================
# Compute NNGP neighbors
#================================================


#' Compute neighbor sets for a Nearest-Neighbor Gaussian Process (NNGP)
#'
#' Computes ordered neighbor sets based on Vecchia's approximation for a set of spatial coordinates.
#' Locations are first ordered by their x-coordinate, and neighbors are selected from previously
#' ordered locations.
#'
#' @param coords Numeric matrix (\eqn{M \times 2}) of spatial coordinates, where each row corresponds to a location.
#' @param k Positive integer. Number of neighbors for each location.
#' @param proj Coordinate reference system (CRS) used to compute distances. Can be an EPSG code
#'   (e.g., 4326) or a CRS object compatible with \code{sf}.
#'
#' @return A list with components:
#' \describe{
#'   \item{coords_sorted}{Matrix of coordinates sorted by x-coordinate}
#'   \item{edist_sorted}{Matrix (\eqn{M \times M}) of pairwise distances between sorted coordinates}
#'   \item{neighbors}{Binary matrix (\eqn{M \times M}) indicating neighbor relationships}
#'   \item{neighbor_idx}{Matrix (\eqn{M \times k}) of neighbor indices for each location}
#'   \item{neighbors_dist}{Matrix (\eqn{M \times k}) of distances to neighbors}
#' }
#'
#' @details
#' Neighbor selection follows a Vecchia-type ordering. For each location \eqn{i},
#' neighbors are chosen from the set of previously ordered locations \eqn{1, \dots, i-1}.
#'
#' Distances are computed using \code{sf::st_distance} and depend on the CRS.
#'
#' The first location has no neighbors. For the first \eqn{k+1} locations, all previous
#' locations are used as neighbors. For subsequent locations, the \eqn{k} nearest
#' neighbors among previously ordered points are selected.
#'
#' @examples
#' coords <- matrix(runif(40), ncol = 2)
#' res <- computeNeighbors(coords, k = 5, proj = 4326)
#'
#' @seealso \code{\link[sf]{st_distance}}
#'
#' @author Fabian Ketwaroo
#' 
#' @export
computeNeighbors <- function(coords,k , proj){
    
    N <- dim(coords)[1]
    # Sort coordinates by x axis 
    sorted_idx <- order(coords[,1]) # in ascending order
    coords_sorted =coords[sorted_idx, ]
    
    
    coords_sorted.sf <- sf::st_as_sf(data.frame(x = coords_sorted[, 1], y = coords_sorted[, 2]),
                                     coords = c('x', 'y'), crs = proj)
    
    # Compute euclidean distance matrix 
    edist <- as.matrix(sf::st_distance(coords_sorted.sf))

    neighbors <- matrix(0,N,N) 
    neighbor_idx <- matrix(0, N, k)
    neighbors_dist <- matrix(0, N, k+1) #
    neighbors_dist[ 1, 1:k ] = edist[1, 1:k]
    
    # the 1st location has no neighbors and then k+1 location has k neighbors to the left
    for (i in  2:(k + 1)) {
        neighbor_idx[i, 1:(i - 1) ] = as.numeric(1:(i -  1))
        neighbors[i, neighbor_idx] = 1
        neighbors_dist[i,  1:(i - 1)  ] <- edist[i, neighbor_idx[i, 1:(i - 1) ]]
        neighbors_dist[ i, i:k] <- edist[i, i:k]
    }
    
    if ((k + 2) <= N) {
  for (i in (k + 2):N) {
    curr_sorted_idx <- order(edist[i, ])
    curr_sorted_idx <- curr_sorted_idx[curr_sorted_idx < i]
    neighbor_idx[i, ] <- curr_sorted_idx[1:min(i, k)]
    neighbors[i, neighbor_idx[i, ]] <- 1
    neighbors_dist[i, 1:k] <- edist[i, neighbor_idx[i, ]]
  }
}
    
    
    return(list(coords_sorted = coords_sorted, edist_sorted = edist, neighbors = neighbors, neighbor_idx = neighbor_idx, neighbors_dist = neighbors_dist))
    
}



#================================================
# Display NNGP neighbors based on ordering 
#================================================

#' Plot neighbor sets for a Nearest-Neighbor Gaussian Process (NNGP)
#'
#' Visualizes the neighbor set for a given location based on an ordered
#' Nearest-Neighbor Gaussian Process (NNGP) constructed using Vecchia's approximation.
#'
#' @param ind Integer index of the focal location (in the sorted order).
#' @param coords.sorted Numeric matrix (\eqn{M \times 2}) of spatial coordinates sorted
#'   according to the ordering used in the NNGP construction. Output from \code{computeNeighbors}.
#' @param neighbor_idx Integer matrix (\eqn{M \times k}) where each row contains the indices
#'   of the neighbors for a given location. Output from \code{computeNeighbors}.
#' @param proj Coordinate reference system (CRS) used for plotting. Can be an EPSG code
#'   or a CRS object compatible with \code{sf}.
#'
#' @return A \code{ggplot2} object showing:
#' \itemize{
#'   \item All spatial locations (black markers)
#'   \item The focal location (red point)
#'   \item Its neighbors (blue points)
#'   \item A vertical reference line at the focal location's x-coordinate
#' }
#'
#' @details
#' The function assumes that coordinates and neighbor indices are based on an
#' ordered representation of the data (e.g., sorted by x-coordinate or another
#' ordering used in Vecchia's approximation).
#'
#' The neighbor set for location \eqn{i} consists of previously ordered locations,
#' as defined by the NNGP construction.
#'
#' @examples
#' coords <- matrix(runif(40), ncol = 2)
#' res <- computeNeighbors(coords, k = 5, proj = 4326)
#' p <- displayNeighbors(
#'   ind = 10,
#'   coords.sorted = res$coords_sorted,
#'   neighbor_idx = res$neighbor_idx,
#'   proj = 4326
#' )
#' print(p)
#'
#' @seealso computeNeighbors
#'
#' @author Fabian Ketwaroo
#' 
#' @export
displayNeighbors <- function( ind, coords.sorted, neighbor_idx, proj ){
    
    coords <-coords.sorted
    
    # Coordinates for plotting
    coords.sf <- sf::st_as_sf(data.frame(x = coords[, 1], y = coords[, 2]),
                              coords = c('x', 'y'), crs = proj)

    
    idx <- neighbor_idx[ind, ]
    idx <- idx[idx > 0]
    
    nlocation <- sf::st_as_sf(
                         data.frame(
                             x = coords[idx, 1],
                             y = coords[idx, 2]
                         ),
                         coords = c("x", "y"),
                         crs = proj
                     )
    
    ##result <- ggplot2::ggplot() +
    ##    ggplot2::geom_sf(data = coords.sf, pch = 2) +
    ##    ggplot2::theme_bw(base_size = 14) + 
    ##    ggplot2::labs(x = "x", y = "y") +
    ##    ggplot2::geom_sf(data = nlocation, col = "blue") +
    ##    ggplot2::geom_point(
    ##                 ggplot2::aes(x = coords[ind, 1], y = coords[ind, 2]),
    ##                 color = "red",
    ##                 size = 3
    ##             ) + 
    ##    ggplot2::geom_vline(
    ##                 xintercept = coords[ind, 1],
    ##                 linetype = "dotted",
    ##                 linewidth = 1
    ##             )

    result <- ggplot() +
        geom_sf(data = coords.sf, pch = 2) +
        theme_bw(base_size = 14) + 
        labs(x = "x", y = "y") +
        geom_sf(data = nlocation, col = "blue") +
        geom_point(
            aes(x = coords[ind, 1], y = coords[ind, 2]),
            color = "red",
            size = 3
        ) + 
        geom_vline(
            xintercept = coords[ind, 1],
            linetype = "dotted",
            linewidth = 1
        )
    
    
    return(result)
    
}
