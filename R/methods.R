# Reported parameters ------------------------------------------------------------

.report_item <- function(it, par, D, dimnames) {
  pp <- .decode_item(it, par, D)
  a <- pp$a
  if (D == 1) {
    out <- c()
    if (length(it$idx_a)) out["a"] <- a
    if (it$fam == "dich") {
      out["b"] <- -pp$d / a
      if (length(it$idx_g)) out["g"] <- pp$g
      if (length(it$idx_u)) out["u"] <- pp$u
    } else if (it$model == "rsm") {
      out["b"] <- par[it$idx_d]
    } else {
      b <- -pp$d / a
      names(b) <- paste0("b", seq_along(b))
      if (it$fam == "gpcm") out["b"] <- mean(b)
      out <- c(out, b)
    }
  } else {
    av <- a[it$load]; names(av) <- paste0("a_", dimnames[it$load])
    out <- av
    if (it$fam == "dich") out["d"] <- pp$d else {
      dd <- pp$d; names(dd) <- paste0("d", seq_along(dd)); out <- c(out, dd)
    }
    if (length(it$idx_g)) out["g"] <- pp$g
    if (length(it$idx_u)) out["u"] <- pp$u
  }
  out
}

.jacobian <- function(f, x, idx, h = 1e-6) {
  f0 <- f(x)
  Jm <- matrix(0, length(f0), length(x))
  for (p in idx) {
    up <- x; up[p] <- up[p] + h; dn <- x; dn[p] <- dn[p] - h
    Jm[, p] <- (f(up) - f(dn)) / (2 * h)
  }
  Jm
}

#' Item parameter estimates
#'
#' @param object An `irtfit` object.
#' @param format `"long"` (one row per parameter, with standard errors) or
#'   `"wide"` (one row per item).
#' @param ... Unused.
#' @return A data frame. For unidimensional models parameters are in the
#'   IRT metric (`a` discrimination, `b` difficulty/thresholds, `g`
#'   guessing); for multidimensional models they are slopes and intercepts.
#' @export
coef.irtfit <- function(object, format = c("long", "wide"), ...) {
  format <- match.arg(format)
  spec <- object$spec; par <- object$par; D <- spec$D
  rows <- list()
  for (j in seq_len(spec$J)) {
    it <- spec$items[[j]]
    f <- function(p) .report_item(it, p, D, colnames(object$Q))
    est <- f(par)
    se <- rep(NA_real_, length(est))
    if (!is.null(object$vcov)) {
      Jm <- .jacobian(f, par, it$idx)
      se <- sqrt(pmax(diag(Jm %*% object$vcov %*% t(Jm)), 0))
    }
    rows[[j]] <- data.frame(item = colnames(object$data)[j], model = it$model,
                            term = names(est), estimate = unname(est), se = se,
                            stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  if (!is.null(spec$rsm_tau_idx)) {
    it <- spec$items[[1]]
    ftau <- function(p) { t0 <- p[spec$rsm_tau_idx]; c(t0, -sum(t0)) }
    est <- ftau(par)
    se <- NA_real_
    if (!is.null(object$vcov)) {
      Jm <- .jacobian(ftau, par, spec$rsm_tau_idx)
      se <- sqrt(pmax(diag(Jm %*% object$vcov %*% t(Jm)), 0))
    }
    out <- rbind(out, data.frame(item = "(rating scale)", model = "rsm",
                                 term = paste0("tau", seq_along(est)), estimate = est,
                                 se = se, stringsAsFactors = FALSE))
  }
  rownames(out) <- NULL
  if (format == "wide") {
    terms <- unique(out$term)
    w <- stats::reshape(out[, c("item", "model", "term", "estimate")], idvar = c("item", "model"),
                        timevar = "term", direction = "wide")
    names(w) <- sub("^estimate\\.", "", names(w))
    rownames(w) <- NULL
    return(w[, c("item", "model", intersect(terms, names(w)))])
  }
  out
}

#' @export
logLik.irtfit <- function(object, ...) {
  structure(object$loglik, df = object$npar, nobs = object$nobs, class = "logLik")
}

#' @export
print.irtfit <- function(x, ...) {
  mods <- table(x$model)
  cat("<irtfit> ", paste(sprintf("%s (%d)", names(mods), mods), collapse = ", "),
      " | ", x$nobs, " persons x ", ncol(x$data), " items | ", x$spec$D,
      " dimension(s)\n", sep = "")
  cat(sprintf("log-likelihood = %.3f  AIC = %.2f  BIC = %.2f  (%s in %d cycles)\n",
              x$loglik, stats::AIC(x), stats::BIC(x),
              if (x$converged) "converged" else "NOT converged", x$iterations))
  invisible(x)
}

#' Summarise a fitted IRT model
#'
#' @param object An `irtfit` object.
#' @param ... Unused.
#' @export
summary.irtfit <- function(object, ...) {
  out <- list(fit = object, coef = coef(object, "wide"), Sigma = object$Sigma,
              reliability = reliability(object))
  class(out) <- "summary.irtfit"
  out
}

#' @export
print.summary.irtfit <- function(x, digits = 3, ...) {
  print(x$fit)
  cat("\nItem parameters:\n")
  cf <- x$coef
  num <- vapply(cf, is.numeric, logical(1))
  cf[num] <- lapply(cf[num], round, digits)
  print(cf, row.names = FALSE)
  cat("\nLatent covariance:\n")
  print(round(x$Sigma, digits))
  cat("\nReliability:\n")
  print(x$reliability, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Likelihood-ratio comparison of nested IRT models
#'
#' @param object,... `irtfit` objects fitted to the same data.
#' @export
anova.irtfit <- function(object, ...) {
  fits <- c(list(object), list(...))
  ll <- vapply(fits, function(f) f$loglik, numeric(1))
  df <- vapply(fits, function(f) f$npar, numeric(1))
  n <- object$nobs
  out <- data.frame(model = vapply(fits, function(f) paste(unique(f$model), collapse = "/"),
                                   character(1)),
                    npar = df, logLik = ll, AIC = -2 * ll + 2 * df,
                    BIC = -2 * ll + log(n) * df)
  o <- order(df); out <- out[o, ]
  out$LR <- c(NA, 2 * diff(out$logLik)); out$df_diff <- c(NA, diff(out$npar))
  out$p <- c(NA, stats::pchisq(out$LR[-1], out$df_diff[-1], lower.tail = FALSE))
  rownames(out) <- NULL
  out
}
