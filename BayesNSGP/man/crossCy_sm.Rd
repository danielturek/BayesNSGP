% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp2Scale.R
\name{crossCy_sm}
\alias{crossCy_sm}
\title{Calculate sparse kernel, core kernel, and determine nonzero entries}
\usage{
crossCy_sm(
  Xdists,
  coords,
  Pcoords,
  d,
  n1,
  n2,
  r0,
  s0,
  cstat_opt,
  normalize,
  bumpLocs,
  rads,
  ampls,
  shps,
  Xdist1_sq,
  Xdist2_sq,
  Xdist12,
  Sigma11,
  Sigma22,
  Sigma12,
  PSigma11,
  PSigma22,
  PSigma12,
  nu,
  log_sigma_vec,
  Plog_sigma_vec
)
}
\arguments{
\item{coords}{N x d matrix of coordinate/input locations}

\item{d}{Scalar; dimension of the spatial domain.}

\item{n1}{Scalar; number of outer products.}

\item{n2}{Scalar; number of bump functions in each outer product.}

\item{r0}{Scalar; length-scale of sparse stationary kernel.}

\item{cstat_opt}{Scalar; determines the compactly supported kernel. See Details.}

\item{bumpLocs}{Array of bump function locations (n2*d x n1)}

\item{rads}{Matrix of bump function radii (n1 x n2; denoted r_{ij})}

\item{ampls}{Matrix of bump function amplitudes (n1 x n2; denoted a_{ij})}

\item{shps}{Matrix of bump function shape parameters (n1 x n2; denoted b_{ij})}

\item{Sigma11}{Vector of length N; contains the 1-1 element of the 
anisotropy process for each station.}

\item{Sigma22}{Vector of length N; contains the 2-2 element of the 
anisotropy process for each station.}

\item{Sigma12}{Vector of length N; contains the 1-2 element of the 
anisotropy process for each station.}

\item{nu}{Scalar; Matern smoothness parameter. \code{nu = 0.5} corresponds 
to the Exponential correlation; \code{nu = Inf} corresponds to the Gaussian
correlation function.}

\item{log_sigma_vec}{Vector of length N; log of the signal standard deviation.}

\item{dists}{N x N matrix of Euclidean distances}

\item{N}{Scalar; number of data measurements.}

\item{dist1_sq}{N x N matrix; contains values of pairwise squared distances
in the x-coordinate.}

\item{dist2_sq}{N x N matrix; contains values of pairwise squared distances
in the y-coordinate.}

\item{dist12}{N x N matrix; contains values of pairwise signed cross-
distances between the x- and y-coordinates. The sign of each element is
important; see \code{nsDist} function for the details of this calculation.
in the x-coordinate.}

\item{lognuggetSD}{Vector of length N; log of the error standard deviation.}
}
\value{
Returns a sparse matrix (N x 3) of the nonzero elements of the product 
between the core and sparse kernel.
}
\description{
\code{Cy_sm} calculates the normalized sparse kernel for a fixed
set of bump function hyperparameters and returns the nonzero entries. Note
that the matrix is calculated and returned in dense format.
}
