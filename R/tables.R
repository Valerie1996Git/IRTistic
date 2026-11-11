# Publication-ready tables ---------------------------------------------------------

#' Publication-ready tables
#'
#' Turns IRTistic results into clean, consistently formatted tables: numbers
#' rounded to a fixed number of decimals, p-values in APA style ("< .001"),
#' estimates optionally combined with their standard errors as "1.23 (0.05)",
#' and a caption and note attached. Tables can be returned as a data frame or
#' rendered as Markdown, HTML or LaTeX (via 'knitr'), and written to `.csv`,
#' `.html`, `.md` or `.tex` files.
#'
#' Supported inputs: `irtfit` (item parameters), `mfrmfit` / `mfrm_report`
#' (facet measurement reports), and the outputs of [item_fit()],
#' [person_fit()], [sx2()], [reliability()], [dif_test()], [dif_irt()],
#' [mfrm_bias()], [classification_accuracy()], [cat_simulate()],
#' [score_table()], [equate_true_score()], as well as any data frame.
#'
#' @param x An object to tabulate.
#' @param digits Number of decimals.
#' @param format `"data.frame"` (default), `"markdown"`, `"html"` or `"latex"`.
#' @param file Optional output file; the format is taken from the extension
#'   (`.csv`, `.html`, `.md`, `.tex`).
#' @param caption Table caption (a sensible default is used).
#' @param ... Further arguments for specific methods: `se` (combine estimates
#'   and SEs, default `TRUE`) for `irtfit`; `facet` for facet reports.
#' @return An `irt_table` (a data frame with `caption` and `note`
#'   attributes), or a 'knitr' table for the other formats. Written files are
#'   returned invisibly.
#' @examples
#' fit <- irt_fit(sim_irt(500, 8, "2pl"), "2pl")
#' irt_table(fit)
#' irt_table(item_fit(fit), format = "markdown")
#' @export
irt_table <- function(x, digits = 2, format = c("data.frame", "markdown", "html", "latex"),
                      file = NULL, caption = NULL, ...) {
  format <- match.arg(format)
  tb <- .as_irt_table(x, digits = digits, ...)
  if (!is.null(caption)) attr(tb, "caption") <- caption
  .emit_table(tb, format, file)
}

.as_irt_table <- function(x, digits, ...) UseMethod(".as_irt_table")

.fmt_num <- function(v, digits) ifelse(is.na(v), "", formatC(v, format = "f", digits = digits))

.fmt_p <- function(p) {
  out <- ifelse(is.na(p), "", ifelse(p < 0.001, "< .001", sub("^0", "", formatC(p, format = "f", digits = 3))))
  out
}

.tbl <- function(df, caption, note = NULL, digits = 2) {
  for (nm in names(df)) {
    v <- df[[nm]]
    if (is.logical(v)) df[[nm]] <- ifelse(is.na(v), "", ifelse(v, "Yes", ""))
    else if (is.numeric(v)) {
      if (grepl("^p$|^p_|_p$|^p\\.", nm)) df[[nm]] <- .fmt_p(v)
      else if (all(is.na(v) | abs(v - round(v)) < 1e-9) && !grepl("measure|theta|estimate|bias|se$", nm))
        df[[nm]] <- ifelse(is.na(v), "", format(round(v), big.mark = ","))
      else df[[nm]] <- .fmt_num(v, digits)
    }
  }
  names(df) <- .pretty_names(names(df))
  rownames(df) <- NULL
  structure(df, caption = caption, note = note, class = c("irt_table", "data.frame"))
}

.pretty_names <- function(n) {
  map <- c(item = "Item", model = "Model", measure = "Measure", se = "SE", n = "N", infit = "Infit MNSQ",
           infit_z = "Infit ZSTD", outfit = "Outfit MNSQ", outfit_z = "Outfit ZSTD", obs_avg = "Obs. avg",
           exp_avg = "Exp. avg", pt_measure = "Pt-measure r", flag = "Flag", person = "Person",
           theta = "\u03b8", raw = "Raw score", chisq = "\u03c7\u00b2", df = "df", p = "p", p_adj = "p (adj.)",
           S_X2 = "S-X\u00b2", RMSEA = "RMSEA", element = "Element", count = "Count", total = "Total score",
           fair_avg = "Fair avg", reliability = "Reliability", separation = "Separation", strata = "Strata",
           index = "Index", focal = "Focal group", delta_b = "\u0394b", signed_area = "Signed area",
           unsigned_area = "Unsigned area", sids = "SIDS", uids = "UIDS", anchor = "Anchor", class = "Class",
           delta_MH = "MH \u0394", alpha_MH = "MH \u03b1", statistic = "Statistic", delta_R2 = "\u0394R\u00b2",
           bias = "Bias", t = "t", obs_score = "Observed", exp_score = "Expected", obs_minus_exp = "Obs - Exp",
           score = "Score", expected_prop = "Expected proportion", score_x = "Form X score",
           score_y = "Form Y equivalent", facet = "Facet", elements = "Elements", mean = "Mean", sd = "SD",
           mean_se = "Mean SE", estimate = "Estimate", term = "Parameter", category = "Category",
           percent = "%", avg_measure = "Avg measure", exp_measure = "Exp measure", threshold = "Threshold",
           threshold_se = "Threshold SE", disordered = "Disordered", group = "Group",
           uniform_chisq = "Uniform \u03c7\u00b2", uniform_p = "Uniform p", nonuniform_chisq = "Non-uniform \u03c7\u00b2",
           nonuniform_p = "Non-uniform p", cut = "Cut", accuracy = "Accuracy", false_positive = "False positive",
           false_negative = "False negative", statistic_name = "Statistic", value = "Value")
  out <- ifelse(n %in% names(map), map[n], n)
  unname(out)
}

#' @export
.as_irt_table.default <- function(x, digits, ...) {
  if (!is.data.frame(x)) stop("Don't know how to tabulate an object of class '", class(x)[1], "'.", call. = FALSE)
  .tbl(as.data.frame(x), caption = NULL, digits = digits)
}

#' @export
.as_irt_table.irtfit <- function(x, digits, se = TRUE, ...) {
  cf <- coef(x)
  terms <- unique(cf$term)
  items <- unique(cf$item)
  out <- data.frame(item = items, model = cf$model[match(items, cf$item)], stringsAsFactors = FALSE)
  for (tm in terms) {
    s <- cf[cf$term == tm, ]
    est <- s$estimate[match(items, s$item)]; sev <- s$se[match(items, s$item)]
    out[[tm]] <- if (se && any(!is.na(sev)))
      ifelse(is.na(est), "", paste0(.fmt_num(est, digits), ifelse(is.na(sev), "", paste0(" (", .fmt_num(sev, digits), ")"))))
    else .fmt_num(est, digits)
  }
  note <- paste0("Marginal maximum likelihood estimates", if (se) "; standard errors in parentheses" else "",
                 ". a = discrimination, b = difficulty / threshold, g = lower asymptote, u = upper asymptote. ",
                 sprintf("N = %s, log-likelihood = %.2f, AIC = %.1f, BIC = %.1f.",
                         format(x$nobs, big.mark = ","), x$loglik, stats::AIC(x), stats::BIC(x)))
  tb <- .tbl(out, caption = "Item parameter estimates", note = note, digits = digits)
  tb
}

#' @export
.as_irt_table.irt_itemfit <- function(x, digits, ...)
  .tbl(as.data.frame(x), "Item fit statistics",
       "MNSQ = mean-square; ZSTD = Wilson-Hilferty standardised mean-square. Productive range 0.5-1.5.", digits)

#' @export
.as_irt_table.irt_personfit <- function(x, digits, ...)
  .tbl(as.data.frame(x), "Person fit statistics", NULL, digits)

#' @export
.as_irt_table.irt_sx2 <- function(x, digits, ...)
  .tbl(as.data.frame(x), "S-X\u00b2 item fit (Orlando & Thissen, 2000)",
       "p-values adjusted with the Benjamini-Hochberg procedure.", digits)

#' @export
.as_irt_table.irt_reliability <- function(x, digits, ...)
  .tbl(as.data.frame(x), "Reliability and separation", NULL, digits)

#' @export
.as_irt_table.irt_dif <- function(x, digits, ...) {
  y <- x; class(y) <- "data.frame"
  .tbl(y, paste0("Differential item functioning (", if (attr(x, "method") == "mh") "Mantel-Haenszel" else "logistic regression", ")"),
       "Class A = negligible, B = moderate, C = large DIF.", digits)
}

#' @export
.as_irt_table.irt_dif_irt <- function(x, digits, ...) {
  y <- x; class(y) <- "data.frame"
  y <- y[, c("item", "focal", "chisq", "df", "p_adj", "delta_b", "signed_area", "sids", "uids", "anchor", "flag")]
  .tbl(y, "IRT-based DIF: Lord's Wald test and Raju's area measures",
       paste0("Reference group: ", attr(x, "reference"),
              ". \u0394b and SIDS are focal minus reference; negative SIDS = harder for the focal group."), digits)
}

#' @export
.as_irt_table.mfrm_bias <- function(x, digits, ...) {
  y <- x; class(y) <- "data.frame"
  .tbl(y, paste("Bias/interaction analysis:", paste(attr(x, "facets"), collapse = " \u00d7 ")),
       "Positive bias = higher scores than expected.", digits)
}

#' @export
.as_irt_table.mfrmfit <- function(x, digits, ...) .as_irt_table(facet_report(x), digits, ...)

#' @export
.as_irt_table.mfrm_report <- function(x, digits, facet = NULL, ...) {
  if (is.null(facet)) return(.tbl(x$summary, "Facet summary", "Separation, strata and reliability as in FACETS; chi-square tests that all elements are equal.", digits))
  if (facet == "categories") return(.tbl(x$categories, "Rating scale category statistics", NULL, digits))
  tab <- x$measures[[facet]]
  if (is.null(tab)) stop("Unknown facet '", facet, "'.", call. = FALSE)
  .tbl(tab[order(-tab$measure), ], paste(facet, "measurement report"),
       "Higher measures = harder to score highly on (except persons). Fair average = expected score for an average element of the other facets.", digits)
}

#' @export
.as_irt_table.irt_classification <- function(x, digits, ...) {
  overall <- data.frame(cut = "Overall", accuracy = x$accuracy, false_positive = NA, false_negative = NA)
  bc <- x$by_cut; bc$cut <- as.character(bc$cut)
  .tbl(rbind(bc, overall), "Classification accuracy and consistency",
       sprintf("%s method. Consistency = %.3f, kappa = %.3f.", if (x$method == "lee") "Lee (2010)" else "Rudner (2005)",
               x$consistency, x$kappa), digits)
}

#' @export
.as_irt_table.cat_sim <- function(x, digits, ...) {
  s <- summary(x)
  rows <- c("Simulees" = s$n, "Mean test length" = s$length[["mean"]], "SD test length" = s$length[["sd"]],
            "Mean final SE" = s$mean_se)
  if (!is.null(s$recovery)) rows <- c(rows, "Bias" = s$recovery[["bias"]], "RMSE" = s$recovery[["rmse"]],
                                      "Correlation with true \u03b8" = s$recovery[["correlation"]])
  if (!is.null(s$classification_accuracy)) rows <- c(rows, "Classification accuracy" = s$classification_accuracy)
  rows <- c(rows, "Maximum exposure rate" = s$exposure[["max"]], "% of bank unused" = s$exposure[["unused_pct"]],
            "Test overlap rate" = s$exposure[["overlap"]])
  df <- data.frame(statistic_name = names(rows), value = .fmt_num(rows, digits), stringsAsFactors = FALSE)
  df$value[1] <- format(s$n, big.mark = ",")
  d <- x$design
  .tbl(df, "CAT simulation results",
       sprintf("Selection: %s; estimator: %s; exposure control: %s.", d$select, d$estimator, d$exposure), digits)
}

.emit_table <- function(tb, format, file) {
  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "csv") { utils::write.csv(.plain_df(tb), file, row.names = FALSE); return(invisible(file)) }
    format <- switch(ext, html = "html", htm = "html", md = "markdown", tex = "latex",
                     stop("Unsupported file extension '.", ext, "'. Use .csv, .html, .md or .tex.", call. = FALSE))
  }
  if (format == "data.frame" && is.null(file)) return(tb)
  if (!requireNamespace("knitr", quietly = TRUE))
    stop("Install the 'knitr' package to render tables as ", format, ".", call. = FALSE)
  kfmt <- switch(format, markdown = "pipe", html = "html", latex = "latex")
  k <- knitr::kable(.plain_df(tb), format = kfmt, caption = attr(tb, "caption"),
                    align = c("l", rep("r", ncol(tb) - 1)), escape = format != "html")
  if (!is.null(attr(tb, "note")))
    k <- structure(c(k, "", paste0(if (format == "latex") "\\textit{Note.} " else "*Note.* ", attr(tb, "note"))),
                   format = attr(k, "format"), class = class(k))
  if (!is.null(file)) {
    txt <- as.character(k)
    if (format == "html") txt <- c(.html_table_head(), txt, "</body></html>")
    writeLines(txt, file, useBytes = TRUE)
    return(invisible(file))
  }
  k
}

.html_table_head <- function() c(
  "<!doctype html><html><head><meta charset='utf-8'>",
  "<style>body{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;color:#0b0b0b;margin:24px}",
  "table{border-collapse:collapse;font-size:14px;font-variant-numeric:tabular-nums}",
  "caption{text-align:left;font-weight:600;padding:6px 0}",
  "th{border-top:2px solid #0b0b0b;border-bottom:1px solid #0b0b0b;padding:6px 10px;text-align:right}",
  "th:first-child,td:first-child{text-align:left}td{padding:4px 10px;text-align:right}",
  "tbody tr:last-child td{border-bottom:2px solid #0b0b0b}</style></head><body>")

#' @export
print.irt_table <- function(x, ...) {
  if (!is.null(attr(x, "caption"))) cat(attr(x, "caption"), "\n\n", sep = "")
  y <- x; class(y) <- "data.frame"
  print(y, row.names = FALSE, right = TRUE)
  if (!is.null(attr(x, "note"))) cat("", strwrap(paste("Note.", attr(x, "note")), width = 90, exdent = 6), sep = "\n")
  invisible(x)
}

.plain_df <- function(tb) {
  y <- tb
  attr(y, "caption") <- NULL; attr(y, "note") <- NULL
  class(y) <- "data.frame"
  y
}
