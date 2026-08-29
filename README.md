# BayesNSGP
R package for Bayesian analysis of non-stationary Gaussian processes using NIMBLE

The package also provides nearest-neighbor Gaussian process (NNGP) components for use directly
in user-written NIMBLE model code, where the spatial process is retained as a latent field:
neighbor-structure construction, a latent-field density and matching simulation function, a
purpose-built Metropolis-Hastings sampler that updates the field one node at a time, and
posterior prediction at unobserved locations. A
[worked example](https://danielturek.github.io/nngp/nngp_RW_sampler/demo.html) fits a spatial
Poisson model with a latent NNGP field and benchmarks the available samplers.
