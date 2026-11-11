#' IRTistic: IRT and many-facet Rasch measurement with artful ggplot2 graphics
#'
#' Fit models with [irt_fit()] (Rasch, 2PL, 3PL, GRM, PCM, GPCM, RSM,
#' unidimensional or confirmatory multidimensional) or [mfrm_fit()]
#' (many-facet Rasch). Then score persons ([irt_scores()]), check fit
#' ([item_fit()], [person_fit()], [facet_report()]), reliability
#' ([reliability()]), local dependence ([local_dependence()]),
#' dimensionality ([residual_pca()]), DIF ([dif_test()]) and rater bias
#' ([mfrm_bias()]), and draw everything with `plot()` or the `plot_*()`
#' family, which all return ggplot objects.
#'
#' @keywords internal
#' @importFrom ggplot2 .data
#' @importFrom stats coef setNames
"_PACKAGE"
