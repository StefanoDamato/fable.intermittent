# Contributing to fable.intermittent

Report bugs and request features in the [issue tracker](https://github.com/StefanoDamato/fable.intermittent/issues).
For anything bigger than a bug fix, please open an issue first.

## Setup

```r
install.packages(c("devtools", "testthat", "roxygen2"))
devtools::install_deps(dependencies = TRUE)
devtools::load_all()
```

Use `load_all()` while you work. It compiles `src/` and gives you the internal
functions without installing the package.

## Adding a model

This package has been created following the `fabletools` [vignette on how to build models](https://fabletools.tidyverts.org/articles/extension_models.html)
. Read it and replicate its structure to implement new methods.

Store one method per file, named after the model in lowercase: `TWEES` is in
[R/twees.R](R/twees.R). Copy [R/empsd.R](R/empsd.R) as a template. It is short and
has every required part.

A model called `XXX` needs:

| Function | Notes |
|---|---|
| `XXX()` | Exported constructor. Check arguments here with `rlang::arg_match()`, not in `train_xxx()`, so bad input fails before fitting. Ends with `new_model_definition()`. |
| `train_xxx()` | Not exported, not documented. Rejects multivariate responses, all-`NA` series, and any `NA`. Returns `structure(list(...), class = "XXX")`. |
| `forecast.XXX()` | Must return a `distributional` vector, not a numeric one. |
| `fitted.XXX()`, `residuals.XXX()` | Simple accessors. |
| `generate.XXX()` | Sets `new_data$.sim` and returns `new_data`. |
| `model_sum.XXX()` | Short label for the mable. |
| `tidy.XXX()` | A tibble with `term` and `estimate`. |
| `report.XXX()` | `cat()` a few lines, then `invisible(object)`. |
| `xxx_no_xreg()` | The `xreg` special. No model here takes regressors, so it just calls `abort()`. |

Two rules:

- `model_sum()` labels use `(d)` for damped and `(u)` for undamped, as in
  [R/twees.R:246](R/twees.R#L246). Models that pick a distribution show the choice,
  like `PARAMSD(nbinom)`. Keep labels short: they go in a mable column.
- Do not write hard-coded constants, such as `1e-4`, in the `XXX.R` file. Add a `.XXX_EPSILON` constant in
  [R/utils.R](R/utils.R).

When the model works, add it to `NEWS.md`, to the table in [README.Rmd](README.Rmd)
(then re-knit `README.md`), to the list in
[vignettes/fable.intermittent.Rmd](vignettes/fable.intermittent.Rmd), and to the
`Description` field in `DESCRIPTION` with its DOI.

## Internal functions

None of these are exported. Inside the package, call them directly. From the
console, run `devtools::load_all()` first.

### R helpers, in [R/utils.R](R/utils.R)

| Helper | Returns |
|---|---|
| `crostons_decomp(y)` | `occurrence` (0/1), `demand` (non-zero values), `intervals` (gaps between them). Used by `MARWAL`, `HSPES`, `PARAMSD`, `VZ`, `WSS`. |
| `deseasonalize(y, period, max_prop_zeros, normalize = FALSE)` | `y_deseasonalized` and `seasons`. `normalize` makes the factors sum to `period`. |
| `get_freq(.data, period)` | The seasonal period, as one integer. |
| `make_hurdle_shifted_distr(distr, pzero)` | Shifts `distr` up by one and zero-inflates it. For hurdle models. |
| `fit_nbinom(y)` | `c(size, prob)` by maximum likelihood. Falls back to moment matching if `nloptr` fails. |
| `fit_tweedie(y, discrete = FALSE)` | `c(mean, dispersion, power)`. With `discrete = TRUE` it needs non-negative integers and also fits the mean. |

Two things to watch:

- `deseasonalize()` returns `seasons = NULL` in three cases: `period <= 1`, too many
  zeros, or seasonal factors that are not positive and finite. Check for `NULL`
  instead of assuming you got factors back.
- In `dampedSES()`, `phi` is not a damping parameter. It weights the series mean:
  `mu[t] = alpha * y[t-1] + phi * mean(y) + (1 - alpha - phi) * mu[t-1]`. Also, if
  `alpha + phi > 1`, the function rescales both instead of failing. So the values
  `nloptr` reports may not be the ones the recursion used. Constrain them in the
  optimiser if that matters.

`R/utils.R` also defines two distributions that `distributional` does not have:
`dist_normal_nonneg()` and `dist_tweedie_discrete()`. If you add another, write
`format`, `density`, `cdf`, `quantile`, `generate`, `mean` and `covariance` for it.
`fabletools` and `fable.bayesRecon` need all of them.

### C++ helpers

Four functions in `src/`, available as R functions through
[R/RcppExports.R](R/RcppExports.R):

| Function | Source | Returns |
|---|---|---|
| `dampedSES(y, mu0, alpha, phi)` | [src/ses.cpp](src/ses.cpp) | Smoothed vector |
| `armaDynamic(y, phi, theta, co)` | [src/arma.cpp](src/arma.cpp) | `list(m, v)`: conditional means and innovations |
| `gammaDynamic(y, a0, b0, w)` | [src/bayesian.cpp](src/bayesian.cpp) | `list(a, b)`: Gamma shape and rate |
| `betaDynamic(y, v, a0, b0, w)` | [src/bayesian.cpp](src/bayesian.cpp) | `list(a, b)`: Beta parameters |

After you change a `// [[Rcpp::export]]` signature:

```r
Rcpp::compileAttributes()
devtools::document()
```

Never edit the `RcppExports` files by hand.

## Documentation

Roxygen with markdown, roxygen2 8.1.0. Run `devtools::document()` and commit the new
`NAMESPACE` and `man/`.

- **Constructor:** full block. Title, description, one `@param` per argument,
  `@references` with a `\doi{}`, `@return A model specification.`, and an
  `@examples` that fits and forecasts. Always cite the paper: every model here
  implements a published method.
- **`forecast`, `fitted`, `residuals`, `generate`:** each gets its own title and its
  own `man/*.Rd`.
- **`report`:** use `@rdname XXX` so it appears on the model page.
- **`model_sum`, `tidy`:** only `@export` and `@importFrom`. No title, so no `Rd`
  file. This is on purpose, please do not add titles.
- **Internal functions:** no roxygen. The only exception is the S3 methods for the
  custom distributions, which need `@noRd` plus `@export` or `@exportS3Method`.

Declare imports with `@importFrom` on the function that uses them, not in one central
place. Guard examples that need a suggested package with
`if (requireNamespace("ggtime", quietly = TRUE))`.

## Tests

One file per model: `tests/testthat/test-<model>.R`.
[tests/testthat/helper.R](tests/testthat/helper.R) gives you `test_data`, a list of
six tsibbles (monthly, quarterly, weekly and daily, with different amounts of
intermittency), plus `base_ts()`, `some_na_ts()`, `all_zero_ts()` and `xreg_ts()` for
edge cases.

Loop over `test_data` and check that the model fits, that `fitted()` and
`residuals()` return one row per observation, that `forecast()` returns an `fbl_ts`
with a distribution column, and that `generate()` works. To test error messages, call
`train_xxx()` directly with `specials = list()`.

```r
devtools::test()
devtools::check()
```

## Getting help

I am happy to help. If something here is unclear, or you get stuck while adding a
model, ask in the [issue tracker](https://github.com/StefanoDamato/fable.intermittent/issues)
or write to me at <stefano.damato@supsi.ch>.
