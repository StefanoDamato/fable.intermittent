# fable 0.3.0

## New features

* Released `NNARMA` (non-negative autoregressive moving average) model by Sbrana (2026).

* Added occurrence smoothing for `TWEES` to determine the dispersion parameter.

* Added a static Tweedie distribution to `PARAMSD`. `distr = "auto"` now
  ranks it as a fifth candidate, in a discretised form obtained by rounding to
  the non-negative integers so that its log-likelihood is a probability mass
  comparable with the count distributions; when it wins the model is reported as
  `PARAMSD(tweedie_discrete)`. `distr = "tweedie"` fits the continuous
  Tweedie, for intermittent series that are not counts.

## Deprecations

* Removed `distr = "mixture"` from `PARAMSD`. Use `distr = "auto"` to select
  a single distribution by AIC/BIC.

## Minor improvements and bug fixes

* Fixed the parameter count used by `PARAMSD`'s information criteria.
  It was taken from `distributional::parameters()`, which for the hurdle
  distributions returns the `dist_inflated` wrapper's fields and charged `hsp`
  three parameters instead of two. Counts are now declared per candidate.

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
