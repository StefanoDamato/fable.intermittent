# Static Tweedie in `STATICDISTR`

Notes on the Tweedie candidate added to `STATICDISTR`, on branch `staticTweedie`.

## What was built

**`dist_tweedie_discrete()` in `R/utils.R`** — a rounded-Tweedie distribution class
following the existing `dist_normal_nonneg` pattern, with `P(Y=0) = F(0.5)` and
`P(Y=k) = F(k+½) − F(k−½)`, plus `format` / `density` / `cdf` / `quantile` /
`generate` / `mean` / `covariance` methods.

Verified as a genuine pmf, not just plumbing: it sums to 1, `cdf` matches
`cumsum(pmf)`, `quantile` inverts the step cdf in both directions, and over
50 000 draws `generate()` reproduces the theoretical mean (1.4837 vs 1.4851) and
`P(0)` (0.2682 vs 0.2673).

**`fit_tweedie(y, discrete = )`** now maximises the *discretised* likelihood when
discretising, rather than plugging continuous estimates into a discrete model —
otherwise the AIC would not be a maximised likelihood. Since the sample mean is
no longer the exact maximiser, all three parameters are optimised. The likelihood
is accumulated over distinct counts weighted by their frequency, which is about
130× cheaper than over the raw series (0.04 s vs 5 s per fit at n = 2000).

**`R/staticdistr.R`** — `"tweedie"` joins the candidate list for both `auto` and
`mixture`, gated by a new `tweedie_discrete = TRUE` argument. `auto` forces
discretisation and rejects `FALSE` with an explanatory error; `mixture` and
explicit `distr = "tweedie"` honour it.

### It works, and the Tweedie genuinely competes

AIC on the test fixtures is now on one scale, and on `ts6` the Tweedie wins
outright:

```
ts6  pois=8647.9  hsp=8161.8  nbinom=7689.9  hsnb=7668.8  tweedie=7665.9  -> selected
```

The default mixture has 5 components and samples only integers;
`tweedie_discrete = FALSE` gives non-integer draws, as intended.

`R CMD check`: 0 errors, 0 warnings. Full test suite green.

## Three things to know

**Explicit `distr = "tweedie"` now defaults to the discretised form**, where in the
previous commit it was continuous. One knob governs all three entry points, which
is probably the right trade, but it is a behaviour change relative to what is
already committed — `tweedie_discrete = FALSE` restores the old result.

**Non-integer data now error under `auto`.** The discretised Tweedie aborts on
non-integer input with a message naming the escape hatch. Since `auto` forces
discretisation, `STATICDISTR(value)` on a non-count series now errors where it
previously returned a (dubious) count fit. That is a hard failure replacing
silent nonsense, but it is a new failure mode.

**`auto` costs about 0.22 s more per series** (0.16 s fit + 0.06 s information
criterion) — roughly 1 s/series on 1825-point pasta SKUs. Caching the pmf across
repeated counts was tried and reverted: `distributional` calls `density()` once
per observation with a length-1 vector, so it was provably dead code. A
`log_likelihood.dist_tweedie_discrete` method could reclaim the 0.06 s, but the
fit dominates.

## Why `P(Y = 0) = F(0.5)`

Because the discretisation is just rounding: `Y = round(X)`.

`Y = 0` exactly when `X` falls in `[0, 0.5)`. A Tweedie with power in (1, 2) has
support `[0, ∞)`, so there is no mass below 0 and the interval's lower edge
contributes nothing:

```
P(Y = 0) = P(0 ≤ X < 0.5) = F(0.5) − F(0⁻) = F(0.5)
```

It is the same `F(k+½) − F(k−½)` rule as every other `k`, with the lower edge
clipped from `−0.5` to `0` — and since `F(−0.5) = F(0⁻) = 0`, that clipping is
free. That is exactly the `lower[kk == 0] <- 0` line in `tweedie_discrete_pmf()`.

Worth noting what this bundles together. `F(0.5)` contains two different pieces
of mass:

- the genuine atom at zero, `P(X = 0) = exp(−μ^(2−p) / (φ(2−p)))` — the
  compound-Poisson "no arrivals" event;
- the continuous mass on `(0, 0.5)` — small positive demands that round down to
  zero.

For μ = 1.2, φ = 1, p = 1.5 the atom alone is 0.1118 while `F(0.5) = 0.3352`, so
most of the zero probability comes from rounding rather than from the atom. That
is a modelling choice, not a derivation: it treats "demand below half a unit" as
an observed zero, which is what you want when the data are counts.
