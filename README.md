
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fable.intermittent <a href="https://github.com/StefanoDamato/fable.intermittent/"><img src="man/figures/logo.png" align="right" height="150" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/StefanoDamato/fable.intermittent/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/StefanoDamato/fable.intermittent/actions/workflows/R-CMD-check.yaml)

[![Coverage
Status](https://coveralls.io/repos/github/StefanoDamato/fable.intermittent/badge.svg?branch=main)](https://coveralls.io/github/StefanoDamato/fable.intermittent?branch=main)

[![CRAN
status](https://www.r-pkg.org/badges/version/fable.intermittent)](https://CRAN.R-project.org/package=fable.intermittent)

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: LGPL (\>=
3)](https://img.shields.io/badge/license-LGPL%20(%3E%3D%203)-yellow.svg)](https://www.gnu.org/licences/lgpl-3.0)
<!-- badges: end -->

The package `fable.intermittent` implements probabilistic methods for
intermittent time series in the [`tidyverts`](https://tidyverts.org/)
framework. The methods follow the
[`fable`](https://fable.tidyverts.org/)/[`fabletools`](https://fabletools.tidyverts.org/)
syntax. To fit the methods, use them as arguments of the `model()`
function. To generate forecasts, pass the fitted models to the
`forecast()` function, following the
[`tidyr`](https://tidyr.tidyverse.org/) pipeline.

The forecasting methods are the following:

| Method | Description |
|----|----|
| `BETANBB()` | Bayesian dynamic negative binomial model with a beta prior on the probability parameter. |
| `EMPSD()` | Empirical resampling baseline that forecasts from the observed distribution. |
| `GAMPOISB()` | Bayesian dynamic Poisson model with a gamma prior on the rate parameter. |
| `HSPES()` | Exponential smoothing model with a hurdle-shifted Poisson forecast distribution. |
| `MARWAL()` | ARMA model with a Markov walk on the occurrence and Gaussian forecast distribution. |
| `NEGBINES()` | Exponential smoothing model with a negative binomial forecast distribution. |
| `NNARMA()` | Non-negative ARMA model with a non-negative Gaussian forecast distribution. |
| `PARAMSD()` | Static count-distribution model that selects among candidate distributions by AIC or BIC. |
| `TWEES()` | Exponential smoothing model with a Tweedie forecast distribution. |
| `VZ()` | Bootstrap method based on Croston decomposition sampling demand sizes and intervals. |
| `WSS()` | Bootstrap method with sampled demand sizes and a Markov-chain for the occurrence. |

The probabilistic forecasts produced by the implemented methods are
[`distributional`](https://github.com/mitchelloharawild/distributional)
objects. Among the predictive distributions used by the methods, there
is the Tweedie distribution, provided by the
[`tweedieDistr`](https://cran.r-project.org/package=tweedieDistr)
package as the `stats`-style `dtweedie()`, `ptweedie()`, `qtweedie()`,
`rtweedie()` functions).

Finally, the package releases four data sets in the
[`tsibble`](https://tsibble.tidyverts.org/) format: `auto`, `pasta`,
`raf`, and `tinyM5`.

## Installation

You can install the **stable** version from CRAN:

``` r
install.packages("fable.intermittent")
```

You can install the **development** version from
[GitHub](https://github.com/StefanoDamato/fable.intermittent):

``` r
# install.packages("devtools")
devtools::install_github("StefanoDamato/fable.intermittent", build_vignettes = TRUE, dependencies = TRUE)
```

## Usage

The package follows the standard `fable` workflow:

1.  Prepare data as a `tsibble`.
2.  Fit forecasting methods with `model()`.
3.  Produce probabilistic forecasts with `forecast()`.

We provide in [this vignette](vignettes/fable.intermittent.Rmd) a simple
usage example; refer to the package documentation for more details on
the methods.

## Acknowledgements

Work on this package was supported by the Swiss National Science
Foundation (SNSF) under grant number 200021_212164 ([“Probabilistic
Forecasting: Global Models, Gaussian Processes and
Hierarchies”](https://data.snf.ch/grants/grant/212164)).

## Contributors

<!-- prettier-ignore-start -->

<!-- markdownlint-disable -->

<table>

<tbody>

<tr>

<td align="center" valign="top" width="20%">

<a href="#">
<img src="https://github.com/StefanoDamato.png" width="100px;" alt="Stefano Damato" style="border-radius:50%;border:1px solid #646464;"/><br />
<sub><b>Stefano Damato</b></sub></a><br /> <sub>(Maintainer)</sub><br />
<a href="mailto:stefano.damato@idsia.ch?subject=[fable.intermittent package]">Email</a>
</td>

<td align="center" valign="top" width="20%">

<a href="#">
<img src="https://github.com/LorenzoZambon.png" width="100px;" alt="Lorenzo Zambon" style="border-radius:50%;border:1px solid #646464;"/><br />
<sub><b>Lorenzo Zambon</b></sub></a><br /> <sub> </sub><br />
<a href="mailto:lorenzo.zambon@idsia.ch?subject=[fable.intermittent package]">Email</a>
</td>

<td align="center" valign="top" width="20%">

<a href="https://dazzimonti.github.io/">
<img src="https://github.com/dazzimonti.png" width="100px;" alt="Dario Azzimonti" style="border-radius:50%;border:1px solid #646464;"/><br />
<sub><b>Dario Azzimonti</b></sub></a><br /> <sub> </sub><br />
<a href="mailto:dario.azzimonti@gmail.com?subject=[fable.intermittent package]">Email</a>
</td>

</tr>

</tbody>

</table>

<!-- markdownlint-restore -->

<!-- prettier-ignore-end -->

## Getting help

If you encounter a bug, please file a minimal reproducible example on
[GitHub](https://github.com/StefanoDamato/fable.intermittent/issues).
