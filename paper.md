---
title: 'Nearest-neighbor Gaussian processes for Bayesian hierarchical models in NIMBLE: Scalable densities, samplers and prediction in `BayesNSGP`'
tags:
  - R
  - NIMBLE
  - Bayesian inference
  - spatial statistics
  - Gaussian processes
  - Markov chain Monte Carlo
  - ecology
authors:
  - name: Fabian R. Ketwaroo
    orcid: 0009-0007-8540-5208
    corresponding: true
    affiliation: 1
  - name: Daniel Turek
    orcid: 0000-0002-1453-1908
    corresponding: false
    affiliation: 2
affiliations:
  - name: Swiss Ornithological Institute, Sempach, Switzerland
    index: 1
  - name: Lafayette College, Easton, PA, United States
    index: 2
date: 18 August 2026
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary
 
Ecologists, epidemiologists and environmental scientists routinely need
to account for the fact that observations taken close together in space
tend to be similar. Gaussian processes (GPs) are the standard
statistical tool for this, but they scale poorly for large datasets:
fitting a GP to $M$ locations requires operations that grow as $M^3$, so
even a few thousand sites can make a model impractical to fit. The
nearest-neighbor Gaussian process [NNGP, @datta2016hierarchical] removes
this bottleneck by letting each location depend on only a small number,
$k$, of nearby locations, which yields a valid probability model whose
cost grows linearly in $M$.
 
This paper describes new NNGP functionality added to the `BayesNSGP`
package [@risser2020bayesian] for the `NIMBLE` probabilistic programming
system [@devalpine2017programming]. The additions let a user build the
neighbor structure for a set of coordinates, simulate a spatial surface,
insert an NNGP-distributed spatial random effect directly into a
hierarchical model written in the BUGS language, sample it with a
purpose-built Markov chain Monte Carlo (MCMC) algorithm, and predict at
unobserved locations. The spatial process is formulated as a
user-defined probability distribution. It can therefore be combined with
any likelihood expressible in the BUGS language: count, occupancy,
capture–recapture, disease-prevalence and integrated models.
 
# Statement of need
 
Hierarchical spatial models in ecology and epidemiology are usually
non-Gaussian: counts, detections and survival outcomes are observed, and
a latent spatial field enters through a link function. In this setting
the spatial random effects cannot be integrated out analytically, so an
MCMC algorithm fitting the model must update the entire field at every
iteration. Two issues then need attention. First, evaluating the spatial
density is $\mathcal{O}(M^3)$ using a full GP model. Second, the
sampling algorithm matters as much as the GP density itself. The
`NIMBLE` package provides conditional autoregressive distributions for
areal data, but no GP approximation for point-referenced data, and its
automatic sampler assignment is necessarily generic: a multivariate node
receives a block Metropolis-Hastings random-walk sampler. This default is
not tailored to a high-dimensional spatial field, where joint proposals
are rarely accepted.
 
A fast density is therefore not sufficient on its own: an efficient
sampler is needed alongside it, and implementations that provide both
are generally tied to a specific model family. The components described
here target researchers who need a scalable spatial random effect as a
single component of a larger, custom hierarchical model, and who need
that component to be both statistically exact and computationally
efficient within the model they have written.
 
# State of the field
 
Several mature implementations of NNGP and related Vecchia-type
approximations exist. `spNNGP` [@finley2022spnngp] is the reference R
implementation: it is highly optimised, offers OpenMP parallelisation,
and handles Gaussian responses and binomial responses via a Pólya-Gamma
augmented Gibbs sampler. Its scope, however, is spatial regression:
users cannot embed its NNGP inside an arbitrary hierarchical model.
`GpGp` [@guinness2018permutation] and `BRISC` [@saha2018brisc] provide
fast maximum-likelihood and bootstrap inference for Vecchia
approximations, but are not Bayesian. `R-INLA` [@rue2009approximate;
@lindgren2011explicit] supports very large latent Gaussian models
through a stochastic partial differential equation representation, but
is restricted to the latent Gaussian model class and does not give the
user control over the sampling algorithm.
 
`BayesNSGP` [@risser2020bayesian] already contains an NNGP
approximation, available through `nsgpModel()` for (non)stationary GP
regression, where `dmnorm_nngp()` serves as the likelihood of the
observed response. `nsgpPredict()` supports this likelihood but operates
on a fitted `nsgpModel()` object. This is a coherent design for that
problem. When the observation model is non-Gaussian the process must
remain a latent field, which calls for a density that can be evaluated
on a latent node, a simulator, a suitable sampler and a prediction
algorithm that conditions on posterior draws; a user who writes such a
model by hand with the existing density receives `NIMBLE`'s default block
sampler. Rather than create a competing package, we contribute to
`BayesNSGP` a complementary, user-facing layer for the stationary case,
to be composed by hand inside BUGS code, together with the simulation,
sampling and prediction machinery that a latent-field formulation
requires. The new functions, `dmnorm_NN_GP()`, `rmnorm_NN_GP()`,
`RW_NN_GP` and `NNGP.pred()`, coexist with the package's existing
functions.
 
# Software design
 
Our contribution comprises three groups of functions and rests on three
design decisions.
 
**A compressed representation instead of sparse matrices.** The NNGP
implies a sparse Cholesky decomposition of the precision matrix
$\tilde{\mathbf{K}}^{-1} = (\mathbf{I}-\mathbf{A})^\top \mathbf{D}^{-1}(\mathbf{I}-\mathbf{A})$,
where $\mathbf{A}$ is strictly lower triangular with at most $k$
non-zero entries per row and $\mathbf{D}$ is diagonal. `NIMBLE`'s
domain-specific language has limited support for sparse matrix types, so
rather than storing $\mathbf{A}$ and $\mathbf{D}$ sparsely,
`computeAD()` packs their non-zero elements into a dense
$M \times (k+1)$ matrix: columns $1{:}k$ hold the neighbor weights, and
column $k+1$ holds the conditional variances. This keeps every operation
inside the subset of the language that `NIMBLE` compiles to C++ and,
because `computeAD()` is a `nimbleFunction`, lets the weights be a
deterministic node that is recomputed whenever the range or variance
parameter is updated.
 
**The process as a distribution, not a wrapper.**
`dmnorm_NN_GP(mu, AD, neighbors.id)` is registered as a user-defined
multivariate distribution for the latent spatial field
$\boldsymbol{\omega} = (\omega_1, \ldots, \omega_M)$, evaluating the
log-density as a sum of standardised conditional residuals via
`computeQF()`, never forming a covariance matrix. The deliberate
trade-off is that the user must construct and pass the neighbor
structures themselves — `computeNeighbors()` returns coordinates sorted
along the horizontal axis [@vecchia1988estimation], a distance matrix,
neighbor indices and neighbor distances. Sorting is exposed rather than
hidden precisely because the user must remap their response vector to
the sorted order, and automatic (under the hood) sorting introduces
potential for error. Distances are computed through `sf`
[@pebesma2018simple], so projected or geographic coordinate reference
systems are both supported. `rmnorm_NN_GP()` simulates the field
sequentially in $\mathcal{O}(Mk)$ operations [@datta2022nearest],
drawing each value from its conditional distribution given its
neighbors. Simulating a GP in the usual way costs $\mathcal{O}(M^3)$, so this is a
practical advantage in itself as well as a source of MCMC starting
values.
 
**A sampler that exploits the directed acyclic graph.** A generic scalar
random-walk sampler must reevaluate the full spatial density after each
individual scalar proposal, which is $\mathcal{O}(M)$.
`sampler_RW_NN_GP()` instead recognises that changing $\omega_i$ alters
only its own conditional density and those of its *reverse* neighbors —
the nodes that condition on it. `get_single_reverse_neighbors()` and
`RWNNGP_setup()` resolve this local graph once, in R, at MCMC
configuration time, and translate it into explicit node addresses; the
compiled sampler then performs an $\mathcal{O}(k^2)$ update, the cost of
which is independent of $M$. Proposal scales adapt following
@shaby2011exploring. The trade-off is that the sampler must be assigned
node by node with the neighbor matrix supplied at sampler assignment.
 
Prediction is handled as a post-processing step by `NNGP.pred()`, which
follows Algorithm 2 of @finley2019efficient: for each new location it
finds the $k$ nearest training neighbors under the same ordering rule,
and draws are generated using the resulting univariate conditional
distribution.
 
# Research impact statement
 
A [reproducible
vignette](https://danielturek.github.io/nngp/nngp_RW_sampler/demo.html)
fits a spatial Poisson generalised linear model with a latent NNGP field
($M = 1000$ locations, $N = 900$ for training, $P = 100$ for validation,
$k = 15$) and benchmarks three sampling strategies. `NIMBLE`'s default
block sampler left roughly 85% of the 900 spatial random effects with
$\widehat{R} > 1.2$. Scalar random-walk samplers converged at every node
but ran about 7.5 times slower, for a median 18-fold gain in effective
sample size per second. The `RW_NN_GP` sampler also converged at every
node while running only about 1.9 times slower than the block sampler,
with a median efficiency roughly 69 times that of the block sampler and
3.8 times that of the scalar sampler. Posterior means matched those of
the generic sampler (median relative bias $-0.093$ versus $-0.094$),
confirming that the local update targets the same posterior, and
out-of-sample predictive counts at the 100 held-out locations had a
median relative bias of $0.031$. Because the per-update cost depends on
$k$ and not $M$, we expect the advantage to grow with dataset size.
 
These components were developed for, and are already used in, the
spatially explicit integrated population model (sIPM) of
@ketwaroo2026sipm, which joins count and capture–recapture data from
many sampling locations to estimate spatio-temporal demographic rates
such as apparent survival and recruitment. There the NNGP models residual
spatial autocorrelation and generates spatial predictions in an analysis
of Gray Catbird (*Dumetella carolinensis*) data from the North American
Breeding Bird Survey and the Monitoring Avian Productivity and
Survivorship program along the eastern coast of North America. This is
the case the design targets: a latent spatial field embedded in a large,
custom, non-Gaussian hierarchical model that the packages discussed above
cannot express, sampled node by node across many thousands of MCMC
iterations. The sIPM first used these components as standalone R
functions, later generalised and released as part of BayesNSGP. The case
study has since been re-run with the package, where `RW_NN_GP` gave the
same posterior means of survival and recruitment as the scalar
random-walk sampler and ran 2.2 times faster (analysis repository:
<https://github.com/Fabian-Ketwaroo/Spaital_Integrated_Population_models>).
 
In accordance with JOSS policy on co-publication, we note that
@ketwaroo2026sipm is a related methodological manuscript currently
available as a preprint; the present paper describes the software, not
those research results.
 
# AI usage disclosure
 
Generative AI was used in preparing this manuscript. Claude (Anthropic)
was used to produce an initial draft of the manuscript text from the
existing package vignette and the Roxygen documentation of the
contributed functions, and for revisions. All statistical and
algorithmic design decisions, the implementation, the benchmarks
reported above and the final text are the work of the human authors, who
reviewed, edited and validated all AI-assisted output and are
responsible for its accuracy.
 
# Acknowledgements
 
We thank the developers of `NIMBLE` for an extensible design that made
this contribution possible, Jaume Badia for testing the code, and
Michael Schaub for reading the manuscript. This work was supported by
the Swiss National Science Foundation (grant 215689). The funder had no
involvement in the design, implementation or reporting of the software.
