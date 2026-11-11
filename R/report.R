#' One-call analysis report
#'
#' Writes a complete, self-contained HTML report of a fitted model: overview
#' and model fit, item parameters, reliability, item fit (infit/outfit and
#' S-X2), local dependence, dimensionality, Wright map, item and test curves,
#' summed-score conversion table and, when a grouping variable is given, an
#' IRT-based DIF analysis. For many-facet models the report follows the
#' FACETS layout: facet summary, a measurement table per facet, rating scale
#' statistics, the variable map, fit, bias/interaction analysis and
#' unexpected responses. Every table and figure in the report can also be
#' produced on its own with the functions it names.
#'
#' Requires the 'rmarkdown' and 'knitr' packages (and pandoc, which ships
#' with RStudio).
#'
#' @param x An `irtfit` or `mfrmfit`.
#' @param file Output HTML file.
#' @param title,author Report title and author.
#' @param group Optional grouping vector for a DIF section (`irtfit` only).
#' @param open Open the report in a browser when done.
#' @return The path to the report, invisibly.
#' @examples
#' \dontrun{
#' fit <- irt_fit(sim_irt(500, 12, "2pl"), "2pl")
#' irt_report(fit, "my-analysis.html", group = rep(c("A", "B"), 250))
#' }
#' @export
irt_report <- function(x, file = "IRTistic-report.html", title = NULL, author = NULL, group = NULL,
                       open = interactive()) {
  for (p in c("rmarkdown", "knitr"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Install the '", p, "' package to build reports.", call. = FALSE)
  if (!inherits(x, c("irtfit", "mfrmfit"))) stop("'x' must be an irtfit or mfrmfit.", call. = FALSE)
  if (is.null(title)) title <- if (inherits(x, "irtfit")) "IRT analysis report" else "Many-facet Rasch analysis report"
  body <- if (inherits(x, "irtfit")) .report_irt(x, !is.null(group)) else .report_mfrm(x)
  yaml <- c("---", sprintf("title: \"%s\"", gsub("\"", "'", title)),
            if (!is.null(author)) sprintf("author: \"%s\"", gsub("\"", "'", author)),
            sprintf("date: \"%s\"", format(Sys.Date(), "%d %B %Y")),
            "output:", "  html_document:", "    toc: true", "    toc_depth: 2", "    df_print: default",
            "---", "",
            "```{r setup, include = FALSE}",
            "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, fig.width = 8, fig.height = 4.8, dpi = 110)",
            "library(IRTistic)",
            "tab <- function(z, ...) irt_table(z, format = 'html', ...)",
            "```", "",
            "<style>body{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;max-width:1000px}",
            "table{font-variant-numeric:tabular-nums;margin-bottom:1em} caption{text-align:left;font-weight:600;color:#0b0b0b}",
            "th,td{padding:4px 12px} th{border-top:2px solid #0b0b0b;border-bottom:1px solid #0b0b0b} h2{margin-top:1.6em}</style>", "")
  rmd <- tempfile(fileext = ".Rmd")
  writeLines(c(yaml, body, "", "---", "",
               sprintf("*Generated with IRTistic %s on %s.*", utils::packageVersion("IRTistic"), format(Sys.time(), "%Y-%m-%d %H:%M"))),
             rmd, useBytes = TRUE)
  env <- new.env(parent = globalenv()); env$x <- x; env$group <- group
  out <- rmarkdown::render(rmd, output_file = basename(file), output_dir = dirname(normalizePath(file, mustWork = FALSE)),
                           envir = env, quiet = TRUE)
  if (open) utils::browseURL(out)
  invisible(out)
}

.chunk <- function(code, opts = "") c(sprintf("```{r%s}", if (nzchar(opts)) paste0(", ", opts) else ""), code, "```", "")

.report_irt <- function(x, has_group) {
  uni <- x$spec$D == 1
  mods <- table(x$model)
  c("## Overview", "",
    .chunk(c(
      "mods <- table(x$model)",
      "cat(sprintf('**%s persons** responded to **%d items** (%s). The model has %d dimension(s) and %d free parameters; estimation %s after %d EM cycles.',",
      "  format(x$nobs, big.mark = ','), ncol(x$data), paste(sprintf('%s: %d', names(mods), mods), collapse = ', '),",
      "  x$spec$D, x$npar, if (x$converged) 'converged' else 'did NOT converge', x$iterations))"), "results = 'asis'"),
    .chunk(c("tab(data.frame(statistic_name = c('Log-likelihood', 'AIC', 'BIC'),",
             "  value = round(c(x$loglik, AIC(x), BIC(x)), 2)), caption = 'Model fit')")),
    "## Item parameters", "", .chunk("tab(x)"),
    "## Reliability and separation", "", .chunk("tab(reliability(x))"),
    "## Item fit", "",
    "Infit and outfit mean-squares near 1 indicate items that fit; values above 1.5 suggest noise (underfit), below 0.5 redundancy (overfit).", "",
    .chunk("tab(item_fit(x))"), .chunk("plot_fit(item_fit(x))"),
    if (uni) c("### S-X\u00b2 item fit", "", .chunk("tab(sx2(x))")),
    "## Local independence", "",
    "Yen's Q3 residual correlations; pairs well above the average Q3 suggest local dependence.", "",
    .chunk(c("ld <- local_dependence(x)", "tab(head(ld$pairs, 10), caption = 'Ten largest residual correlations')")),
    .chunk("plot_ld(ld)", "fig.height = 6.5"),
    if (uni) c("## Dimensionality", "",
               "Principal components of standardised residuals: a first-contrast eigenvalue below about 2 is consistent with unidimensionality.", "",
               .chunk(c("pc <- residual_pca(x)", "cat(sprintf('First contrast eigenvalue: **%.2f** (%.1f%% of residual variance).', pc$eigenvalues[1], 100 * pc$proportion[1]))"), "results = 'asis'"),
               .chunk("plot_pca(pc)")),
    "## Person-item map", "", .chunk(if (uni) "plot_wright(x)" else "plot_wright(x, dim = 1)", "fig.height = 6"),
    "## Item and test functions", "",
    .chunk(sprintf("plot_icc(x, items = seq_len(min(12, ncol(x$data)))%s)", if (uni) "" else ", dim = 1"), "fig.height = 6"),
    .chunk(sprintf("plot_info(x%s)", if (uni) "" else ", dim = 1"), "fig.height = 5.5"),
    if (uni) c(.chunk("plot_tcc(x)"), "## Summed-score conversion table", "", .chunk("tab(score_table(x), digits = 3)")),
    if (has_group && uni) c("## Differential item functioning", "",
                            "Items were calibrated separately in each group, linked by Stocking-Lord on purified anchors, and compared with Lord's Wald test.", "",
                            .chunk(c("dd <- dif_irt(x, group)", "tab(dd)")), .chunk("plot_dif(dd)"),
                            .chunk("if (any(dd$flag)) plot_dif_curves(dd)")))
}

.report_mfrm <- function(x) {
  c("## Overview", "",
    .chunk(c(
      "cat(sprintf('**%s ratings** of **%d %s** by facets %s, %s model with %d categories (%s). Estimation %s.',",
      "  format(x$nobs, big.mark = ','), length(x$levels[[x$person]]), x$person,",
      "  paste(sprintf('%s (%d)', x$facets, x$ctx$nf), collapse = ', '), toupper(x$model), x$ctx$K,",
      "  paste(x$categories, collapse = ', '), if (x$converged) 'converged' else 'did NOT converge'))"), "results = 'asis'"),
    .chunk(c("rep <- facet_report(x)", "tab(rep)")),
    "## Variable map", "", .chunk("plot_facets(x)", "fig.height = 6.5"),
    "## Facet measurement reports", "",
    .chunk(c("for (f in x$facets) print(tab(rep, facet = f))"), "results = 'asis'"),
    "## Rating scale", "", .chunk("tab(rep, facet = 'categories')"), .chunk("plot_categories(x)"),
    "## Fit", "", .chunk("plot_fit(rep)", "fig.height = 4.5"),
    "## Bias / interaction analysis", "",
    .chunk(c("if (length(x$facets) >= 2) {", "  b <- mfrm_bias(x, x$facets[1:2])",
             "  print(tab(head(b, 15), caption = 'Largest bias terms'))", "}"), "results = 'asis'"),
    .chunk("if (length(x$facets) >= 2) plot_bias(b)", "fig.height = 5.5"),
    "## Unexpected responses", "",
    .chunk(c("u <- unexpected(x, 2.5)", "if (nrow(u)) tab(head(u, 20), caption = 'Ratings with |standardised residual| >= 2.5') else cat('None.')"),
           "results = 'asis'"))
}
