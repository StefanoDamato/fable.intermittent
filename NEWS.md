# fable 0.3.0

## New features

* Released `tinyM5` data set, a subset of the M5 dataset.

* Released `NNARMA` (non-negative autoregressive moving average) model by Sbrana (2026).

* Added occurrence smoothing for `TWEES` to determine the dispersion parameter.

## Minor improvements and bug fixes

* Changed optimisation strategy for damped exponential smoothing models.

* Fixed recursion in `BETANBB` and `GAMPOISB`.

* Implemented non-negative Gaussian forecast distribution for ARMA-based models.

* Increased test coverage.

# fable.intermittent 0.2.0

## New features

* Released `pasta` hierarchical data set.

## Deprecations

* Dropped `dist_tweedie()` and `stats`-like Tweedie functions; moved them to standalone `tweedieDistr` package.

# fable.intermittent 0.1.1

## New features

* Added `report()` and `tidy()` methods for all model classes.

* Details included in `model_sum()` for ES-based models and `PARAMSD()`.

## Minor improvements and bug fixes

* Fixed an error in the Tweedie quantile method (`quantile.dist_tweedie()`).

* `MARWAL()` forecast distributions are now truncated at zero.

# fable.intermittent 0.1.0

* Initial CRAN submission.
