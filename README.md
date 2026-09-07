# BayesNSGP

R package for Bayesian analysis of non-stationary Gaussian processes using NIMBLE.

The package provides fully Bayesian, nonstationary Gaussian process (GP) modelling, in which
covariance parameters are allowed to vary over space. Inference is carried out by Markov chain
Monte Carlo (MCMC) via the [nimble](https://r-nimble.org) package. Two approximate likelihoods
are available for large spatial datasets: the nearest-neighbor Gaussian process (NNGP) and the
sparse general Vecchia (SGV) approximation.

## New in version 0.3.0

The package now provides NNGP components for use directly inside user-written NIMBLE model
code, where the spatial process is retained as a latent field rather than integrated into the
likelihood of the observed data:

* **Neighbor-structure construction** — order coordinates and build the conditioning sets,
  with a plotting function for checking them
* **A latent-field density and matching simulation function** — an NNGP distribution that can
  be placed on a latent node in BUGS model code, and a simulator that draws a realisation of
  the field in O(Mk) operations
* **A purpose-built Metropolis-Hastings sampler** — updates the field one node at a time by
  exploiting the NNGP graph, so the cost of each update is independent of the number of
  locations
* **Posterior prediction at unobserved locations** — conditioning on posterior draws of the
  latent field

Together these support non-Gaussian hierarchical models — counts, occupancy,
capture–recapture and integrated models — in which the latent field must be sampled
explicitly.

## Installation

### Dependencies

`BayesNSGP` requires **R >= 3.4.0** and the following packages, all available from CRAN and
installed automatically by the commands below:

| Type | Packages |
| --- | --- |
| Depends | `nimble` |
| Imports | `FNN`, `Matrix`, `methods`, `StatMatch`, `sf`, `ggplot2` |

`nimble` compiles model code to C++, so a working C++ toolchain is required. See the
[NIMBLE installation guide](https://r-nimble.org/download) for platform-specific instructions
(Rtools on Windows, Xcode command line tools on macOS, `g++` on Linux).

### From CRAN (recommended)

The CRAN release is version 0.3.0, which includes the NNGP components described above:

```r
install.packages("BayesNSGP")
library(BayesNSGP)
```

### From GitHub (development version)

The package source lives in the `BayesNSGP/` subdirectory of this repository, so `subdir` must
be given:

```r
install.packages("remotes")            # if not already installed
remotes::install_github("danielturek/BayesNSGP", subdir = "BayesNSGP")
```

## Core functionality

Every function below has a help page, accessible with `?` in R (for example
`?dmnorm_NN_GP`) or `help(package = "BayesNSGP")` for the full index.

### Nonstationary GP regression

| Function | Purpose |
| --- | --- |
| `nsgpModel()` | Build a (non)stationary GP regression model, with `likelihood` set to `"fullGP"`, `"NNGP"` or `"SGV"` |
| `nsgpPredict()` | Posterior prediction at unobserved locations from a fitted `nsgpModel()` object |

### NNGP components for latent spatial fields

| Function | Purpose |
| --- | --- |
| `computeNeighbors()` | Order coordinates and find the `k` nearest neighbors of each location; returns sorted coordinates, distances and neighbor indices |
| `computeAD()` | Compute the sparse Cholesky factor of the NNGP precision matrix in compressed form, for use as a deterministic node in model code |
| `dmnorm_NN_GP()` | NNGP density for a latent spatial field, usable as a distribution in BUGS model code |
| `rmnorm_NN_GP()` | Simulate a realisation of the field in O(Mk) operations, for synthetic data or MCMC initial values |
| `sampler_RW_NN_GP()` | Adaptive Metropolis-Hastings sampler that updates the field one node at a time, exploiting the NNGP graph; assigned with `addSampler(type = "RW_NN_GP")` |
| `NNGP.pred()` | Posterior prediction at unobserved locations, conditioning on posterior draws of the latent field |

Supporting functions: `computeC()`, `computeQF()`, `expcov()`, `RWNNGP_setup()` and
`get_single_reverse_neighbors()`.

## Worked example

A [reproducible vignette](https://danielturek.github.io/nngp/nngp_RW_sampler/demo.html) fits a
spatial Poisson generalised linear model with a latent NNGP field, walks through neighbor
construction, model specification, sampler assignment and prediction, and benchmarks the
available samplers against NIMBLE's defaults.

## Getting help and contributing

Questions, bug reports and pull requests are welcome — see
[CONTRIBUTING.md](CONTRIBUTING.md) for how to report a problem, where to ask usage questions,
and how to propose a code change. Please also review the
[Code of Conduct](CODE_OF_CONDUCT.md).
* **Maintainer:** Daniel Turek (danielturek@gmail.com).

## Citation

If you use the package, please cite:

Risser, M. D. and Turek, D. (2020). Bayesian inference for high-dimensional nonstationary
Gaussian processes. *Journal of Statistical Computation and Simulation*.
[doi:10.1080/00949655.2020.1792472](https://doi.org/10.1080/00949655.2020.1792472)

If you use the NNGP components specifically (neighbor-structure construction, the latent-field
density and simulator, the `RW_NN_GP` sampler, or `NNGP.pred()`), please also cite:

Ketwaroo, F. R. and Turek, D. (2026). Nearest-neighbor Gaussian processes for Bayesian
hierarchical models in NIMBLE: Scalable densities, samplers and prediction in BayesNSGP.
*Journal of Open Source Software* (submitted).

## License

GPL-3. See [LICENSE](LICENSE).
