# IRTistic

**IRT + artistic: one R package for IRT and Rasch measurement, from model to publication-ready plot.**

Most IRT work in R today means stitching several packages together: one to
estimate, another for Rasch fit statistics, another for Wright maps, and
nothing at all in R that reproduces the FACETS workflow for rater data.
`IRTistic` does all of it with one consistent interface and a single
`ggplot2` graphics layer.

## What it does

| Area | Functions |
|---|---|
| **Models** (MML-EM, written from scratch) | `irt_fit()` — Rasch, 2PL, 3PL, graded response, partial credit, generalized partial credit, rating scale; mixed item types; confirmatory multidimensional via a Q-matrix |
| **Many-facet Rasch** (FACETS-style) | `mfrm_fit()` — any number of facets, rating-scale or partial-credit thresholds; `facet_report()` with measures, SEs, fair averages, infit/outfit, separation, reliability, fixed chi-square, category statistics |
| **Person scores** | `irt_scores()` — EAP, MAP, ML, WLE |
| **Fit** | `item_fit()`, `person_fit()` — infit/outfit MNSQ and ZSTD, point-measure correlations, observed vs expected averages |
| **Reliability** | `reliability()` — EAP reliability, person/item separation and strata, Cronbach's alpha |
| **Assumptions** | `local_dependence()` (Yen's Q3), `residual_pca()` (first contrast) |
| **DIF** | `dif_test()` — Mantel-Haenszel with ETS classes; logistic / ordinal logistic regression with Jodoin-Gierl classes |
| **Rater effects** | `mfrm_bias()` (bias/interaction terms), `unexpected()` (large residuals) |
| **Model comparison** | `anova()`, `AIC()`, `BIC()`, `logLik()` |
| **Simulation** | `sim_irt()`, `sim_mfrm()` |
| **Graphics** (all return ggplot objects) | `plot_wright()`, `plot_facets()` (vertical ruler), `plot_icc()`, `plot_info()`, `plot_tcc()`, `plot_fit()` (bubble chart), `plot_ld()`, `plot_pca()`, `plot_dif()`, `plot_categories()`, `plot_bias()`, `theme_irt()` |

## Install

```r
# from the source tarball
install.packages("IRTistic_0.1.0.tar.gz", repos = NULL, type = "source")
# or from a local folder / GitHub checkout
devtools::install("path/to/IRTistic")
```

## Quick tour

```r
library(IRTistic)

# --- IRT / Rasch ---------------------------------------------------------
dat <- sim_irt(500, nitems = 12, model = "rasch")
fit <- irt_fit(dat, model = "rasch")      # or "2pl", "3pl", "grm", "pcm", "gpcm", "rsm"
summary(fit)                              # parameters, latent variance, reliability
item_fit(fit)                             # infit / outfit / ZSTD
irt_scores(fit, "WLE")

plot(fit, "wright")                       # person-item map
plot(fit, "icc")                          # faceted item curves
plot(fit, "info")                         # test information + SE
plot(fit, "fit")                          # fit bubble chart

# --- Many-facet Rasch (FACETS) -------------------------------------------
ratings <- sim_mfrm(n_person = 100, n_rater = 6, n_criterion = 4)
mf <- mfrm_fit(ratings, person = "person",
               facets = c("rater", "criterion"), score = "score")
facet_report(mf)                          # FACETS-style measurement tables
plot(mf, "facets")                        # vertical ruler
plot(mf, "bias", facets = c("rater", "criterion"))
```

Because every plot is a ggplot, customise freely:

```r
plot(fit, "icc") + ggplot2::facet_wrap(~item, ncol = 6) + ggplot2::labs(title = "My test")
```

## Validation

The test suite checks parameter recovery for every model, that Rasch
difficulties and infit/outfit statistics match `eRm`, that the one-facet
many-facet model reproduces the IRT rating scale model's likelihood exactly,
and that DIF screening finds planted DIF. During development the Rasch MML
estimates were also checked against `lme4::glmer()` (equal log-likelihood
and difficulties within 0.02).

## Conventions

* Responses are recoded to 0, 1, ... per item; original codes are kept in
  `fit$labels` and used for expected-score plots.
* MML identification: latent mean 0. Variances of dimensions measured only by
  Rasch-family items are estimated; otherwise fixed at 1.
* In many-facet models non-person facets are centred and the person mean is
  estimated (the FACETS convention). Higher facet measures mean harder to
  score highly on (more difficult criteria, more severe raters).
* Standard errors come from the observed information matrix (numerical
  Hessian of the marginal log-likelihood), treating the latent distribution
  as known.

## Roadmap

Exploratory multidimensional models with rotation, bifactor/testlet models,
S-X² item fit and M2 global fit, IRT linking and equating, anchoring of
facet elements, and DIF via multiple-group calibration.
