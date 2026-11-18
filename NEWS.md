# IRTistic 0.1.0

## Updates (November 2026)

New features:

* `dif_test()` can now purify the matching score with `purify = TRUE`, for Mantel-Haenszel and for logistic regression matched on the total score. Flagged items are dropped from the total score and the analysis is repeated until the flagged set is stable. The result shows which items were anchors and the items flagged at each round.
* `irt_report()` can now write Word (`.docx`) reports as well as HTML. The format follows the file extension, or set it with `format = "word"`.
* `irt_table()` can write a table straight to a Word file: `irt_table(x, file = "table.docx")`.

Improvements:

* Clearer error messages for common data problems: an item with only one observed value, an ID or total-score column left in the data, text responses, and persons with no responses.
* `dif_test()` checks for missing group values and says how many persons Mantel-Haenszel drops for missing responses.
* Mantel-Haenszel returns `NA` instead of an error or infinite value when an item has no usable strata.
* Mantel-Haenszel results are checked against `stats::mantelhaen.test()` in the test suite.
* Clearer installation instructions.

## First release

* Models: Rasch, 2PL, 3PL, 4PL, GRM, PCM, GPCM, RSM; mixed item types; confirmatory multidimensional models; many-facet Rasch models.
* Scoring: EAP, MAP, ML, WLE; Lord-Wingersky summed-score tables.
* Diagnostics: infit/outfit, S-X², reliability and separation, Q3, residual PCA, classification accuracy and consistency.
* DIF: Mantel-Haenszel, logistic regression, Lord's Wald test with anchor purification, Raju's areas, SIDS/UIDS.
* Linking and equating: mean-mean, mean-sigma, Haebara, Stocking-Lord; true-score equating.
* Adaptive testing: design, live sessions, Monte Carlo and post-hoc simulation, exposure control, content balancing, classification stopping.
* Graphics, publication tables (`irt_table()`) and one-call HTML reports (`irt_report()`).
