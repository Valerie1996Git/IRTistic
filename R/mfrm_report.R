.mfrm_rows <- function(fit, theta) {
  ctx <- fit$ctx
  pp <- .mfrm_unpack(ctx, fit$par)
  eta <- .mfrm_eta(ctx, pp$beta)
  P <- do.call(cbind, .mfrm_probs(ctx, matrix(theta[ctx$p] + eta), pp$tau)$P)
  k <- 0:(ctx$K - 1)
  E <- as.vector(P %*% k)
  dev <- outer(-E, k, "+")
  list(E = E, W = rowSums(P * dev^2), C = rowSums(P * dev^4), eta = eta, P = P)
}

.fair_avg <- function(fit, measure, facet) {
  ctx <- fit$ctx; pp <- .mfrm_unpack(ctx, fit$par)
  k <- 0:(ctx$K - 1)
  groups <- seq_len(ctx$S)
  sapply(seq_along(measure), function(i) {
    lin <- if (facet == fit$person) measure[i] else fit$mu - measure[i]
    gs <- if (!is.null(fit$step_facet) && facet == fit$step_facet) i else groups
    mean(vapply(gs, function(s) {
      ct <- c(0, cumsum(pp$tau[s, ]))
      num <- k * lin - ct; p <- exp(num - max(num)); p <- p / sum(p)
      sum(p * k)
    }, numeric(1)))
  }) + min(fit$categories)
}

#' Facet measurement reports (FACETS-style)
#'
#' Produces one measurement table per facet, including persons, in the
#' layout familiar from FACETS: element measure and standard error,
#' observed count and score, observed and fair average, infit and outfit
#' mean-squares with ZSTD, plus separation, reliability and the fixed
#' (all-same) chi-square for each facet.
#'
#' @param fit An `mfrmfit` object.
#' @param person_method Person measure used for fit statistics and the
#'   person report: `"WLE"` (default), `"ML"` or `"EAP"`.
#' @return A list of class `mfrm_report` with `measures` (named list of data
#'   frames, one per facet), `summary` (one row per facet) and `persons`.
#' @export
facet_report <- function(fit, person_method = c("WLE", "ML", "EAP")) {
  person_method <- match.arg(person_method)
  ctx <- fit$ctx
  pm <- .mfrm_person_measures(fit, person_method)
  rw <- .mfrm_rows(fit, pm$theta)
  pp <- .mfrm_unpack(ctx, fit$par)
  se_all <- .mfrm_se(fit)
  x <- ctx$x
  facets_all <- c(fit$person, fit$facets)
  index_list <- c(list(ctx$p), ctx$e); names(index_list) <- facets_all
  tabs <- list(); summ <- list(); pos <- 0
  for (f in facets_all) {
    by <- index_list[[f]]
    fs <- .fit_stats(x, rw$E, rw$W, rw$C, by)
    if (f == fit$person) {
      meas <- pm$theta; se <- pm$se
      ext <- fs$obs_avg == 0 | fs$obs_avg == ctx$K - 1
      fs[ext, c("infit", "infit_z", "outfit", "outfit_z")] <- NA
    } else {
      meas <- pp$beta[[f]]
      n <- length(meas)
      se <- if (is.null(se_all)) rep(NA, n) else se_all[pos + seq_len(n)]
      pos <- pos + n
    }
    tab <- data.frame(element = fit$levels[[f]], measure = meas, se = se,
                      count = fs$n, total = fs$obs_avg * fs$n + min(fit$categories) * fs$n,
                      obs_avg = fs$obs_avg + min(fit$categories),
                      fair_avg = .fair_avg(fit, meas, f),
                      infit = fs$infit, infit_z = fs$infit_z,
                      outfit = fs$outfit, outfit_z = fs$outfit_z,
                      stringsAsFactors = FALSE)
    tab$flag <- .flag_fit(tab$infit, tab$outfit)
    tabs[[f]] <- tab
    ok <- is.finite(tab$se) & tab$se < 5
    sep <- if (sum(ok) > 1) .separation(tab$measure[ok], tab$se[ok]) else c(NA, NA, NA)
    w <- 1 / tab$se[ok]^2
    chi <- sum(w * (tab$measure[ok] - sum(w * tab$measure[ok]) / sum(w))^2)
    summ[[f]] <- data.frame(facet = f, elements = nrow(tab), mean = mean(meas), sd = stats::sd(meas),
                            mean_se = mean(tab$se[ok]), separation = sep[2],
                            strata = sep[3], reliability = sep[1], chisq = chi,
                            df = sum(ok) - 1,
                            p = stats::pchisq(chi, sum(ok) - 1, lower.tail = FALSE))
  }
  # rating scale / category statistics
  if (!is.null(se_all)) tau_se <- matrix(se_all[pos + seq_len(ctx$S * (ctx$K - 1))], ctx$S, byrow = TRUE)
  else tau_se <- matrix(NA, ctx$S, ctx$K - 1)
  cats <- list()
  for (s in seq_len(ctx$S)) {
    rows <- ctx$sg == s
    mlin <- pm$theta[ctx$p] + rw$eta
    cnt <- tabulate(x[rows] + 1L, ctx$K)
    avg_m <- vapply(0:(ctx$K - 1), function(k) mean(mlin[rows & x == k]), numeric(1))
    exp_m <- vapply(0:(ctx$K - 1), function(k) sum(rw$P[rows, k + 1] * mlin[rows]) /
                      sum(rw$P[rows, k + 1]), numeric(1))
    cats[[s]] <- data.frame(group = if (is.null(fit$step_facet)) "all" else fit$levels[[fit$step_facet]][s],
                            category = fit$categories, count = cnt,
                            percent = 100 * cnt / sum(cnt), avg_measure = avg_m,
                            exp_measure = exp_m,
                            threshold = c(NA, pp$tau[s, ]), threshold_se = c(NA, tau_se[s, ]))
    cats[[s]]$disordered <- c(FALSE, FALSE, diff(pp$tau[s, ]) < 0)[seq_len(ctx$K)]
  }
  out <- list(measures = tabs, summary = do.call(rbind, summ),
              categories = do.call(rbind, cats), fit = fit)
  rownames(out$summary) <- NULL
  class(out) <- "mfrm_report"
  out
}

#' @export
print.mfrm_report <- function(x, digits = 2, ...) {
  rnd <- function(d) { num <- vapply(d, is.numeric, logical(1)); d[num] <- lapply(d[num], round, digits); d }
  for (f in names(x$measures)) {
    if (f == x$fit$person) next
    cat("\n== ", f, " measurement report ==\n", sep = "")
    tab <- x$measures[[f]]
    print(rnd(tab[order(-tab$measure), ]), row.names = FALSE)
  }
  cat("\n== Facet summary ==\n")
  print(rnd(x$summary), row.names = FALSE)
  cat("\n== Rating scale (category) statistics ==\n")
  print(rnd(x$categories), row.names = FALSE)
  cat("\nPerson measures: use $measures$", x$fit$person, "\n", sep = "")
  invisible(x)
}

#' @export
print.mfrmfit <- function(x, ...) {
  cat("<mfrmfit> many-facet Rasch (", toupper(x$model),
      if (!is.null(x$step_facet)) paste0(" by ", x$step_facet), ")\n", sep = "")
  cat(sprintf("%d ratings | %s: %d | %s\n", x$nobs, x$person, length(x$levels[[x$person]]),
              paste(sprintf("%s: %d", x$facets, x$ctx$nf), collapse = " | ")))
  cat(sprintf("Categories: %s | person mean = %.3f, SD = %.3f\n",
              paste(x$categories, collapse = ","), x$mu, x$sigma))
  cat(sprintf("log-likelihood = %.3f  AIC = %.2f  BIC = %.2f  (%s in %d cycles)\n",
              x$loglik, stats::AIC(x), stats::BIC(x),
              if (x$converged) "converged" else "NOT converged", x$iterations))
  invisible(x)
}

#' @export
summary.mfrmfit <- function(object, ...) facet_report(object)

#' @export
logLik.mfrmfit <- function(object, ...)
  structure(object$loglik, df = object$npar, nobs = object$nobs, class = "logLik")

#' @export
coef.mfrmfit <- function(object, ...) {
  pp <- .mfrm_unpack(object$ctx, object$par)
  se <- .mfrm_se(object)
  out <- do.call(rbind, lapply(names(pp$beta), function(f)
    data.frame(facet = f, element = object$levels[[f]], estimate = pp$beta[[f]])))
  th <- data.frame(facet = "threshold",
                   element = paste0(rep(if (is.null(object$step_facet)) "tau" else
                     object$levels[[object$step_facet]], each = object$ctx$K - 1), "_",
                     rep(seq_len(object$ctx$K - 1), object$ctx$S)),
                   estimate = as.vector(t(pp$tau)))
  out <- rbind(out, th)
  out$se <- if (is.null(se)) NA else se
  out
}

#' Bias / interaction analysis between facets
#'
#' Estimates a bias (interaction) term for every combination of elements of
#' two (or more) facets, holding all other measures fixed, as in the FACETS
#' bias analysis. Positive bias means observed scores are higher than the
#' model expects (e.g. a rater is more lenient for that criterion).
#'
#' @param fit An `mfrmfit` object.
#' @param facets Two or more facet names (may include the person facet).
#' @param min_count Minimum number of ratings for a combination to be
#'   reported.
#' @return A data frame of class `mfrm_bias`.
#' @export
mfrm_bias <- function(fit, facets, min_count = 2) {
  allf <- c(fit$person, fit$facets)
  if (!all(facets %in% allf)) stop("Unknown facet(s).", call. = FALSE)
  pm <- .mfrm_person_measures(fit, "WLE")
  rw <- .mfrm_rows(fit, pm$theta)
  x <- fit$ctx$x
  keys <- lapply(facets, function(f) as.character(fit$data[[f]]))
  key <- do.call(paste, c(keys, sep = "\r"))
  agg <- function(v) tapply(v, key, sum)
  n <- agg(rep(1, length(x))); res <- agg(x - rw$E); info <- agg(rw$W)
  bias <- res / info; se <- 1 / sqrt(info)
  sp <- do.call(rbind, strsplit(names(n), "\r", fixed = TRUE))
  out <- data.frame(sp, count = as.vector(n), obs_score = as.vector(agg(x)) + min(fit$categories) * as.vector(n),
                    exp_score = as.vector(agg(rw$E)) + min(fit$categories) * as.vector(n),
                    obs_minus_exp = as.vector(res / n),
                    bias = as.vector(bias), se = as.vector(se), stringsAsFactors = FALSE)
  names(out)[seq_along(facets)] <- facets
  out$t <- out$bias / out$se
  out$df <- out$count - 1
  out$p <- 2 * stats::pt(-abs(out$t), pmax(out$df, 1))
  out <- out[out$count >= min_count, ]
  out <- out[order(-abs(out$t)), ]; rownames(out) <- NULL
  attr(out, "facets") <- facets
  class(out) <- c("mfrm_bias", "data.frame")
  out
}

#' @export
print.mfrm_bias <- function(x, n = 15, digits = 2, ...) {
  y <- x; class(y) <- "data.frame"
  cat("Bias/interaction: ", paste(attr(x, "facets"), collapse = " x "), " (",
      sum(abs(y$t) >= 2, na.rm = TRUE), " of ", nrow(y), " terms with |t| >= 2)\n", sep = "")
  num <- vapply(y, is.numeric, logical(1)); y[num] <- lapply(y[num], round, digits)
  print(utils::head(y, n), row.names = FALSE)
  invisible(x)
}

#' Unexpected responses
#'
#' Lists individual ratings whose standardised residual exceeds a cutoff
#' (FACETS Table 4).
#' @param fit An `mfrmfit` object.
#' @param cutoff Absolute standardised residual cutoff.
#' @export
unexpected <- function(fit, cutoff = 3) {
  pm <- .mfrm_person_measures(fit, "WLE")
  rw <- .mfrm_rows(fit, pm$theta)
  z <- (fit$ctx$x - rw$E) / sqrt(rw$W)
  out <- fit$data
  out$expected <- rw$E + min(fit$categories); out$residual <- fit$ctx$x - rw$E; out$std_residual <- z
  out <- out[abs(z) >= cutoff, ]
  out[order(-abs(out$std_residual)), ]
}

#' Simulate many-facet rating data
#'
#' @param n_person,n_rater,n_criterion Numbers of elements.
#' @param ncat Number of rating categories (scored 1..ncat).
#' @param rater_sd,criterion_sd SDs of rater severity and criterion difficulty.
#' @param thresholds Andrich thresholds (length ncat - 1, centred).
#' @param design `"complete"` (every rater rates everyone) or a proportion
#'   of ratings to keep at random.
#' @return A long data frame with columns person, rater, criterion, score
#'   and attribute `"truth"`.
#' @export
sim_mfrm <- function(n_person = 100, n_rater = 5, n_criterion = 4, ncat = 5,
                     rater_sd = 0.5, criterion_sd = 0.7, thresholds = NULL,
                     design = "complete") {
  theta <- stats::rnorm(n_person)
  rater <- stats::rnorm(n_rater, 0, rater_sd); rater <- rater - mean(rater)
  crit <- stats::rnorm(n_criterion, 0, criterion_sd); crit <- crit - mean(crit)
  if (is.null(thresholds)) thresholds <- seq(-1.5, 1.5, length.out = ncat - 1)
  d <- expand.grid(person = seq_len(n_person), rater = seq_len(n_rater),
                   criterion = seq_len(n_criterion))
  if (is.numeric(design)) d <- d[stats::runif(nrow(d)) < design, ]
  lin <- theta[d$person] - rater[d$rater] - crit[d$criterion]
  ct <- c(0, cumsum(thresholds)); k <- 0:(ncat - 1)
  num <- outer(lin, k) - matrix(ct, nrow(d), ncat, byrow = TRUE)
  P <- exp(num - apply(num, 1, max)); P <- P / rowSums(P)
  u <- stats::runif(nrow(d))
  d$score <- rowSums(u > t(apply(P, 1, cumsum))[, -ncat, drop = FALSE]) + 1L
  d$person <- sprintf("P%03d", d$person); d$rater <- sprintf("R%02d", d$rater)
  d$criterion <- sprintf("C%d", d$criterion)
  rownames(d) <- NULL
  attr(d, "truth") <- list(theta = theta, rater = rater, criterion = crit, thresholds = thresholds)
  d
}
