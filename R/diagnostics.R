# Expected scores, variances and kurtosis for every person x item ---------------
.expected <- function(fit, theta) {
  X <- fit$data; N <- nrow(X); J <- ncol(X); D <- fit$spec$D
  theta <- matrix(theta, N, D)
  E <- W <- C <- matrix(NA_real_, N, J, dimnames = dimnames(X))
  for (j in seq_len(J)) {
    it <- fit$spec$items[[j]]
    P <- .item_probs(it$fam, .decode_item(it, fit$par, D), theta)
    k <- 0:(it$K - 1)
    E[, j] <- P %*% k
    dev <- outer(-E[, j], k, "+")
    W[, j] <- rowSums(P * dev^2)
    C[, j] <- rowSums(P * dev^4)
  }
  list(E = E, W = W, C = C)
}

.extreme_persons <- function(fit) {
  X <- fit$data
  maxraw <- as.vector((!is.na(X)) %*% vapply(fit$spec$items, function(it) it$K - 1, numeric(1)))
  raw <- rowSums(X, na.rm = TRUE)
  raw == 0 | raw == maxraw
}

.default_theta <- function(fit) {
  rasch <- !any(vapply(fit$spec$items, function(it) length(it$idx_a) > 0, logical(1)))
  if (fit$spec$D == 1 && rasch) "WLE" else "EAP"
}

.theta_matrix <- function(fit, theta) {
  if (is.null(theta)) theta <- .default_theta(fit)
  if (is.character(theta)) {
    sc <- irt_scores(fit, theta)
    cols <- if (fit$spec$D == 1) "theta" else colnames(fit$Q)
    return(as.matrix(sc[, cols, drop = FALSE]))
  }
  as.matrix(theta)
}

.zstd <- function(ms, q) {
  q <- sqrt(pmax(q, 1e-12))
  (ms^(1 / 3) - 1) * (3 / q) + q / 3
}

.fit_stats <- function(x, E, W, C, by) {
  ok <- !is.na(x)
  y <- (x - E)
  z2 <- y^2 / W
  agg <- function(v) tapply(v[ok], by[ok], sum)
  n <- agg(rep(1, length(x)))
  outfit <- agg(z2) / n
  infit <- agg(y^2) / agg(W)
  q_out <- agg(C / W^2) / n^2 - 1 / n
  q_in <- agg(C - W^2) / agg(W)^2
  data.frame(n = as.vector(n), infit = as.vector(infit),
             infit_z = as.vector(.zstd(infit, q_in)), outfit = as.vector(outfit),
             outfit_z = as.vector(.zstd(outfit, q_out)),
             obs_avg = as.vector(agg(x) / n), exp_avg = as.vector(agg(E) / n))
}

.flag_fit <- function(infit, outfit, lower = 0.5, upper = 1.5) {
  ifelse(pmax(infit, outfit) > upper, "underfit",
         ifelse(pmin(infit, outfit) < lower, "overfit", "productive"))
}

.item_measures <- function(fit) {
  cf <- coef(fit)
  items <- colnames(fit$data)
  m <- se <- rep(NA_real_, length(items))
  if (fit$spec$D == 1) for (j in seq_along(items)) {
    s <- cf[cf$item == items[j], ]
    if ("b" %in% s$term) { m[j] <- s$estimate[s$term == "b"]; se[j] <- s$se[s$term == "b"] }
    else { bb <- grepl("^b[0-9]+$", s$term); m[j] <- mean(s$estimate[bb]) }
  }
  data.frame(item = items, measure = m, se = se)
}

#' Item fit statistics (infit / outfit)
#'
#' Computes Rasch-style information-weighted (infit) and unweighted (outfit)
#' mean-square statistics, their Wilson-Hilferty standardised values
#' (ZSTD), point-measure correlations and observed versus expected average
#' scores. They are available for every model in the package.
#'
#' @param fit An `irtfit` object.
#' @param theta Person measures to condition on: a method name passed to
#'   [irt_scores()] or a numeric vector/matrix. Defaults to WLE for
#'   unidimensional Rasch-family models and EAP otherwise.
#' @param range Mean-square range regarded as productive for measurement.
#' @return A data frame with one row per item.
#' @export
item_fit <- function(fit, theta = NULL, range = c(0.5, 1.5)) {
  th <- .theta_matrix(fit, theta)
  ex <- .expected(fit, th)
  X <- fit$data
  X[.extreme_persons(fit), ] <- NA  # extreme scores carry no fit information
  by <- matrix(seq_len(ncol(X)), nrow(X), ncol(X), byrow = TRUE)
  fs <- .fit_stats(X, ex$E, ex$W, ex$C, by)
  pm <- vapply(seq_len(ncol(X)), function(j)
    suppressWarnings(stats::cor(X[, j], th[, which(fit$Q[j, ] != 0)[1]], use = "complete.obs")),
    numeric(1))
  out <- cbind(.item_measures(fit), fs, pt_measure = pm)
  out$flag <- .flag_fit(out$infit, out$outfit, range[1], range[2])
  out
}

#' Person fit statistics
#'
#' @inheritParams item_fit
#' @return A data frame with person measures and infit/outfit statistics.
#' @export
person_fit <- function(fit, theta = NULL, range = c(0.5, 1.5)) {
  th <- .theta_matrix(fit, theta)
  ex <- .expected(fit, th)
  X <- fit$data
  by <- matrix(seq_len(nrow(X)), nrow(X), ncol(X))
  fs <- .fit_stats(X, ex$E, ex$W, ex$C, by)
  out <- data.frame(person = if (!is.null(rownames(X))) rownames(X) else seq_len(nrow(X)))
  if (ncol(th) == 1) out$theta <- th[, 1] else out <- cbind(out, th)
  out <- cbind(out, fs[, c("n", "infit", "infit_z", "outfit", "outfit_z")])
  ext <- .extreme_persons(fit)
  out[ext, c("infit", "infit_z", "outfit", "outfit_z")] <- NA
  out$flag <- .flag_fit(out$infit, out$outfit, range[1], range[2])
  out
}

.separation <- function(est, se) {
  ok <- is.finite(est) & is.finite(se)
  obs <- stats::var(est[ok]); mse <- mean(se[ok]^2)
  true <- max(obs - mse, 0)
  G <- sqrt(true / mse)
  c(reliability = true / obs, separation = G, strata = (4 * G + 1) / 3)
}

#' Reliability and separation indices
#'
#' Reports marginal (EAP) reliability, Rasch person and item separation
#' and reliability (based on WLE person measures and item SEs, for
#' unidimensional models), and Cronbach's alpha.
#'
#' @param fit An `irtfit` object.
#' @return A data frame.
#' @export
reliability <- function(fit) {
  rows <- list()
  D <- fit$spec$D; dn <- colnames(fit$Q)
  eap <- fit$posterior %*% fit$theta
  psd2 <- fit$posterior %*% fit$theta^2 - eap^2
  for (d in seq_len(D)) {
    v <- stats::var(eap[, d]); r <- v / (v + mean(psd2[, d]))
    G <- sqrt(r / (1 - r))
    rows[[length(rows) + 1]] <- data.frame(index = paste0("Person EAP", if (D > 1) paste0(" (", dn[d], ")")),
                                           reliability = r, separation = G, strata = (4 * G + 1) / 3)
  }
  if (D == 1) {
    w <- irt_scores(fit, "WLE")
    keep <- !(w$raw == 0 | abs(w$se) > 5)
    s <- .separation(w$theta[keep], w$se[keep])
    rows[[length(rows) + 1]] <- data.frame(index = "Person WLE", t(s))
    im <- .item_measures(fit)
    if (all(is.finite(im$se))) rows[[length(rows) + 1]] <- data.frame(index = "Item", t(.separation(im$measure, im$se)))
  }
  for (d in seq_len(D)) {
    Xd <- fit$data[, fit$Q[, d] != 0, drop = FALSE]
    Xd <- Xd[stats::complete.cases(Xd), , drop = FALSE]
    k <- ncol(Xd)
    alpha <- if (k > 1 && nrow(Xd) > 2)
      k / (k - 1) * (1 - sum(apply(Xd, 2, stats::var)) / stats::var(rowSums(Xd))) else NA
    rows[[length(rows) + 1]] <- data.frame(index = paste0("Cronbach alpha", if (D > 1) paste0(" (", dn[d], ")")),
                                           reliability = alpha, separation = NA, strata = NA)
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

.std_residuals <- function(fit, theta = NULL) {
  th <- .theta_matrix(fit, theta)
  ex <- .expected(fit, th)
  (fit$data - ex$E) / sqrt(ex$W)
}

#' Local dependence (Yen's Q3)
#'
#' @inheritParams item_fit
#' @param cutoff Flag pairs whose Q3 exceeds the average Q3 by this amount
#'   (Christensen et al.'s relative criterion).
#' @return A list with the Q3 matrix, the average off-diagonal Q3 and a
#'   data frame of item pairs sorted by Q3.
#' @export
local_dependence <- function(fit, theta = NULL, cutoff = 0.2) {
  R <- .std_residuals(fit, theta)
  q3 <- suppressWarnings(stats::cor(R, use = "pairwise.complete.obs"))
  ut <- upper.tri(q3)
  avg <- mean(q3[ut], na.rm = TRUE)
  pairs <- data.frame(item1 = rownames(q3)[row(q3)[ut]], item2 = colnames(q3)[col(q3)[ut]],
                      Q3 = q3[ut])
  pairs$Q3_star <- pairs$Q3 - avg
  pairs$flag <- pairs$Q3_star > cutoff
  pairs <- pairs[order(-pairs$Q3), ]; rownames(pairs) <- NULL
  out <- list(Q3 = q3, mean_Q3 = avg, pairs = pairs)
  class(out) <- "irt_ld"
  out
}

#' @export
print.irt_ld <- function(x, n = 10, ...) {
  cat(sprintf("Yen's Q3: mean off-diagonal = %.3f; %d pair(s) flagged\n",
              x$mean_Q3, sum(x$pairs$flag, na.rm = TRUE)))
  print(utils::head(x$pairs, n), digits = 3, row.names = FALSE)
  invisible(x)
}

#' Principal components analysis of standardised residuals
#'
#' A Winsteps-style check of unidimensionality: after removing the Rasch/IRT
#' dimension, the first residual contrast should look like noise (an
#' eigenvalue below about 2 items' worth of variance is a common rule of
#' thumb).
#'
#' @inheritParams item_fit
#' @return A list with eigenvalues, variance proportions and item loadings
#'   on the first contrasts.
#' @export
residual_pca <- function(fit, theta = NULL) {
  R <- .std_residuals(fit, theta)
  Rm <- suppressWarnings(stats::cor(R, use = "pairwise.complete.obs"))
  Rm[is.na(Rm)] <- 0
  e <- eigen(Rm, symmetric = TRUE)
  ncon <- min(3, ncol(R))
  load <- e$vectors[, seq_len(ncon), drop = FALSE] %*% diag(sqrt(pmax(e$values[seq_len(ncon)], 0)), ncon)
  load <- sweep(load, 2, sign(colSums(load)) + (colSums(load) == 0), "*")
  out <- list(eigenvalues = e$values, proportion = e$values / sum(e$values),
              loadings = data.frame(item = colnames(R), measure = .item_measures(fit)$measure,
                                    setNames(as.data.frame(load), paste0("contrast", seq_len(ncon)))))
  class(out) <- "irt_pca"
  out
}

#' @export
print.irt_pca <- function(x, ...) {
  cat("Residual PCA eigenvalues (first 5):\n")
  print(round(utils::head(x$eigenvalues, 5), 3))
  cat(sprintf("First contrast: eigenvalue %.2f (%.1f%% of residual variance)\n",
              x$eigenvalues[1], 100 * x$proportion[1]))
  invisible(x)
}
