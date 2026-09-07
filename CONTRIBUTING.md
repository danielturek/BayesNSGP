# Contributing to BayesNSGP

Contributions are welcome, whether that means reporting a problem, suggesting a feature,
improving documentation or submitting code. This document explains how to do each.

Everyone taking part is expected to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Seeking support

If you have a question about how to use the package — model specification, sampler
assignment, interpreting output — please start with:

1. The help pages: `?nsgpModel`, `?dmnorm_NN_GP`, `?sampler_RW_NN_GP`, or
   `help(package = "BayesNSGP")` for the full index.
2. The [worked example](https://danielturek.github.io/nngp/nngp_RW_sampler/demo.html), which
   covers neighbor construction, model specification, sampler assignment and prediction.
3. For questions about NIMBLE itself rather than this package, the
   [NIMBLE users mailing list](https://groups.google.com/g/nimble-users) is the better venue.

If those don't answer it, open an issue on the
[issue tracker](https://github.com/danielturek/BayesNSGP/issues) with the `question` label. We
aim to respond within two weeks, though this is a research project maintained alongside other
work, so please allow for delays.

## Reporting issues

Please report bugs on the [issue tracker](https://github.com/danielturek/BayesNSGP/issues).
Before opening a new issue, search the existing ones in case it has already been reported.

A useful bug report includes:

* **A minimal reproducible example** — the smallest piece of code that shows the problem,
  using simulated data where possible so it runs without external files.
* **What you expected to happen**, and what happened instead. Paste the exact error message or
  incorrect output rather than describing it.
* **The output of `sessionInfo()`**, which records your R version, platform and the versions of
  `BayesNSGP` and `nimble`.
* **Whether the problem occurs in uncompiled or compiled code**, if you can tell. NIMBLE
  behaviour sometimes differs between the two, and knowing which narrows the cause considerably.

For problems that look like incorrect statistical results rather than errors, please say what
you compared against — an analytic result, a different sampler, or a simulation study.

## Suggesting features

Open an issue describing what you would like to do and why the current functions do not allow
it. For additions to the NNGP components, it helps to say what model structure you have in
mind, since the design is intended to compose with arbitrary hierarchical models.

## Contributing code

For small fixes — typos, documentation corrections, clear bugs — a pull request is welcome
without prior discussion. For anything larger, please open an issue first so the approach can
be agreed before you invest time in it.

To submit a change:

1. Fork the repository and create a branch from `master`.
2. Make your changes. Note that the package source is in the `BayesNSGP/` subdirectory, so
   install your working copy with
   `remotes::install_local("BayesNSGP")` or
   `remotes::install_github("<your-fork>", subdir = "BayesNSGP")`.
3. Documentation is generated with `roxygen2` from comments in the `R/` files. Edit those
   comments rather than the `.Rd` files in `man/`, then regenerate with
   `roxygen2::roxygenise("BayesNSGP")`.
4. Check the package still builds and passes checks: `R CMD build BayesNSGP` followed by
   `R CMD check BayesNSGP_*.tar.gz`.
5. Open a pull request against `master`, describing what the change does and why. Link any
   related issue.

Some conventions worth knowing:

* Functions intended for use inside model code must be written in the subset of R that NIMBLE
  compiles. If you are unsure, check that your function works both uncompiled and after
  `compileNimble()`.
* New user-facing functions need a `roxygen2` block with `@param`, `@return` and, where
  useful, a runnable `@examples` section.
* Please keep changes focused. Several small pull requests are easier to review than one large
  one.

## Maintainer

The package is maintained by Daniel Turek (danielturek@gmail.com). Issues and pull requests are
the preferred route for anything relating to the package, since they leave a record others can
find.
