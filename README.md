# IRTistic <img src="man/figures/logo.png" align="right" height="120" alt="" />

**IRT + artistic: item response theory and Rasch measurement in R, from model to publication-ready plot and table.**

R has dozens of excellent IRT packages, but a real project usually means
stitching six of them together: one to estimate, one for Rasch fit
statistics, one for Wright maps, one for DIF, one for linking, one for
adaptive testing, then copying numbers into a manuscript by hand. IRTistic
brings that whole workflow under **one consistent interface**, with one
graphics system and one table system.

## Features

| | |
|---|---|
| **Models** | Rasch, 2PL, 3PL, 4PL, graded response, partial credit, generalized partial credit, rating scale; mixed item types; confirmatory multidimensional (`irt_fit()`), all estimated by marginal maximum likelihood with code written from scratch |
| **Many-facet Rasch** | FACETS-style rater models with measurement reports, fair averages, separation, bias/interaction analysis and unexpected responses (`mfrm_fit()`, `facet_report()`, `mfrm_bias()`) |
| **Scoring** | EAP, MAP, ML, WLE (`irt_scores()`); summed-score conversion tables (`score_table()`) |
| **Fit and assumptions** | infit/outfit + ZSTD for every model, S-X² (`item_fit()`, `person_fit()`, `sx2()`); Q3 local dependence, residual PCA (`local_dependence()`, `residual_pca()`) |
| **Reliability and decisions** | EAP/WLE reliability, separation, strata, alpha (`reliability()`); classification accuracy and consistency, Lee and Rudner methods (`classification_accuracy()`) |
| **DIF** | Mantel-Haenszel, logistic/ordinal regression (`dif_test()`); Lord's Wald test with anchor purification, Raju areas, SIDS/UIDS, DTF, any number of groups (`dif_irt()`) |
| **Linking and equating** | mean-mean, mean-sigma, Haebara, Stocking-Lord (`irt_link()`); IRT true-score equating (`equate_true_score()`) |
| **Adaptive testing** | MFI, MPWI, KL, MEPV selection; EAP/MAP/ML/WLE; randomesque, Sympson-Hetter, a-stratification; content balancing; SE and classification stopping; live sessions and Monte Carlo or post-hoc simulation (`cat_design()`, `cat_simulate()`, `cat_start()`) |
| **Graphics** | Wright maps, FACETS rulers, ICC/category curves, information, TCC, fit bubbles, Q3 heatmaps, DIF curves, bias heatmaps, CAT traces, exposure and conditional precision plots; every one a ggplot |
| **Output** | APA-style tables as data frame, Markdown, HTML, LaTeX or CSV (`irt_table()`); a complete HTML analysis report in one call (`irt_report()`) |

## Installation

```r
# install.packages("remotes")
remotes::install_github("Valerie1996Git/IRTistic")
```

## A two-minute tour

```r
library(IRTistic)

dat <- sim_irt(800, nitems = 15, model = "2pl")
fit <- irt_fit(dat, model = "2pl")

irt_table(fit)               # item parameters with SEs, ready for a manuscript
item_fit(fit); sx2(fit)      # infit/outfit and S-X2
plot(fit, "wright")          # person-item map
plot(fit, "icc")             # item curves

group <- rep(c("A", "B"), 400)
dif_irt(fit, group)          # Lord's test with purified anchors

sim <- cat_simulate(cat_design(item_bank(fit), se_target = 0.3), theta = rnorm(500))
summary(sim); plot(sim)      # how would an adaptive version perform?

irt_report(fit, "report.html", group = group)   # everything above, in one HTML file
```

Rater data:

```r
ratings <- sim_mfrm(n_person = 100, n_rater = 6, n_criterion = 4)
mf <- mfrm_fit(ratings, "person", c("rater", "criterion"), "score")
facet_report(mf)
plot(mf, "facets")
```

## Validation

The test suite (over 120 checks) covers parameter recovery for every model,
agreement with eRm for Rasch difficulties and infit/outfit, agreement with
lme4 for Rasch MML, recovery of known linking constants, nominal Type I
error for S-X² and Lord's DIF test, Lee classification accuracy against
simulation, and exact equality of a one-facet many-facet model with the
rating scale model.

## Credits

IRTistic stands on the shoulders of the authors of mirt, TAM, eRm, ltm,
catR, mirtCAT, difR, lordif, DFIT, equateIRT, SNSequate, plink, WrightMap,
ggmirt, irtQ, cacIRT, immer, mfrmr and Linacre's FACETS and Winsteps
programs, and on the psychometric literature they implement. The
[References and credits](https://valerie1996git.github.io/IRTistic/articles/credits.html)
article cites every method and package. The
[Why IRTistic?](https://valerie1996git.github.io/IRTistic/articles/why-irtistic.html)
article maps the package against Choi and Asilkalkan's (2019) review of 45
IRT packages.

## Citation

```r
citation("IRTistic")
```

## Author

Valerie Ofori Aboah
