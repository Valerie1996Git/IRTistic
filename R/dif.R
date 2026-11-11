#' Differential item functioning
#'
#' Screens items for DIF between two or more groups.
#'
#' * `"mh"`: Mantel-Haenszel for dichotomous items (two groups), matched on
#'   total score, with the ETS delta effect size and A/B/C classification.
#' * `"lr"`: logistic regression (binary items) or proportional-odds
#'   ordinal logistic regression (polytomous items), matched on the latent
#'   trait (EAP from `fit`) or on the total score. Reports uniform,
#'   non-uniform and total DIF likelihood-ratio tests and the change in
#'   Nagelkerke R-squared with the Jodoin-Gierl classification.
#'
#' @param fit An `irtfit` object (or a response matrix / data frame when
#'   `method = "mh"` or `match = "total"`).
#' @param group A vector of group memberships, one per person.
#' @param method `"mh"` or `"lr"`.
#' @param match `"theta"` (default when `fit` is a model) or `"total"`.
#' @param reference Reference group level for MH (defaults to the first).
#' @param p.adjust Multiple-comparison adjustment passed to [stats::p.adjust()].
#' @return A data frame of class `irt_dif`, one row per item.
#' @export
dif_test <- function(fit, group, method = c("lr", "mh"), match = NULL,
                     reference = NULL, p.adjust = "BH") {
  method <- match.arg(method)
  if (inherits(fit, "irtfit")) { X <- fit$data } else { X <- .prep_data(fit)$X; fit <- NULL }
  group <- as.factor(group)
  if (length(group) != nrow(X)) stop("'group' must have one value per person.", call. = FALSE)
  if (is.null(match)) match <- if (!is.null(fit) && method == "lr") "theta" else "total"
  out <- if (method == "mh") .dif_mh(X, group, reference) else .dif_lr(X, group, fit, match)
  out$p_adj <- stats::p.adjust(out$p, p.adjust)
  attr(out, "method") <- method
  class(out) <- c("irt_dif", "data.frame")
  out
}

.dif_mh <- function(X, group, reference) {
  if (nlevels(group) != 2) stop("Mantel-Haenszel needs exactly two groups.", call. = FALSE)
  if (any(apply(X, 2, max, na.rm = TRUE) > 1))
    stop("Mantel-Haenszel is implemented for dichotomous items; use method = 'lr'.", call. = FALSE)
  if (is.null(reference)) reference <- levels(group)[1]
  ref <- group == reference
  cc <- stats::complete.cases(X)
  X <- X[cc, , drop = FALSE]; ref <- ref[cc]
  total <- rowSums(X)
  res <- lapply(seq_len(ncol(X)), function(j) {
    num <- den <- A <- EA <- VA <- 0
    for (s in unique(total)) {
      k <- total == s
      a <- sum(X[k & ref, j] == 1); b <- sum(X[k & ref, j] == 0)
      c <- sum(X[k & !ref, j] == 1); d <- sum(X[k & !ref, j] == 0)
      n <- a + b + c + d
      if (n < 2 || (a + b) == 0 || (c + d) == 0) next
      num <- num + a * d / n; den <- den + b * c / n
      A <- A + a; EA <- EA + (a + b) * (a + c) / n
      VA <- VA + (a + b) * (c + d) * (a + c) * (b + d) / (n^2 * (n - 1))
    }
    alpha <- num / den
    chi <- (abs(A - EA) - 0.5)^2 / VA
    delta <- -2.35 * log(alpha)
    p <- stats::pchisq(chi, 1, lower.tail = FALSE)
    cls <- if (abs(delta) >= 1.5 && p < 0.05) "C" else if (abs(delta) >= 1 && p < 0.05) "B" else "A"
    data.frame(item = colnames(X)[j], alpha_MH = alpha, delta_MH = delta, statistic = chi,
               df = 1, p = p, class = cls)
  })
  do.call(rbind, res)
}

.nagelkerke <- function(ll, ll0, n) {
  cs <- 1 - exp(2 * (ll0 - ll) / n)
  cs / (1 - exp(2 * ll0 / n))
}

.dif_lr <- function(X, group, fit, match) {
  if (match == "theta") {
    if (is.null(fit)) stop("match = 'theta' requires a fitted model.", call. = FALSE)
    if (fit$spec$D > 1) stop("Use match = 'total' for multidimensional models.", call. = FALSE)
    m <- as.vector(fit$posterior %*% fit$theta)
  } else m <- as.vector(scale(rowSums(X, na.rm = TRUE)))
  res <- lapply(seq_len(ncol(X)), function(j) {
    ok <- !is.na(X[, j])
    d <- data.frame(y = X[ok, j], m = m[ok], g = droplevels(group[ok]))
    binary <- length(unique(d$y)) == 2
    d$y <- if (binary) as.integer(d$y == max(d$y)) else factor(d$y, ordered = TRUE)
    fitm <- function(f) {
      if (binary) stats::glm(f, family = stats::binomial, data = d)
      else MASS::polr(f, data = d, Hess = FALSE)
    }
    m0 <- fitm(y ~ 1)
    m1 <- fitm(y ~ m); m2 <- fitm(y ~ m + g); m3 <- fitm(y ~ m * g)
    ll <- vapply(list(m0, m1, m2, m3), function(z) as.numeric(stats::logLik(z)), numeric(1))
    dfs <- vapply(list(m1, m2, m3), function(z) attr(stats::logLik(z), "df"), numeric(1))
    n <- nrow(d)
    r2 <- .nagelkerke(ll[2:4], ll[1], n)
    dR2 <- r2[3] - r2[1]
    stat <- 2 * (ll[4] - ll[2]); df <- dfs[3] - dfs[1]
    cls <- if (dR2 >= 0.07) "C" else if (dR2 >= 0.035) "B" else "A"
    p <- stats::pchisq(stat, df, lower.tail = FALSE)
    if (p >= 0.05) cls <- "A"
    data.frame(item = colnames(X)[j],
               uniform_chisq = 2 * (ll[3] - ll[2]),
               uniform_p = stats::pchisq(2 * (ll[3] - ll[2]), dfs[2] - dfs[1], lower.tail = FALSE),
               nonuniform_chisq = 2 * (ll[4] - ll[3]),
               nonuniform_p = stats::pchisq(2 * (ll[4] - ll[3]), dfs[3] - dfs[2], lower.tail = FALSE),
               statistic = stat, df = df, p = p, delta_R2 = dR2, class = cls)
  })
  do.call(rbind, res)
}

#' @export
print.irt_dif <- function(x, digits = 3, ...) {
  cat("DIF analysis (", attr(x, "method"), "): ", sum(x$class != "A"),
      " item(s) with B/C-level DIF\n", sep = "")
  y <- x; class(y) <- "data.frame"
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], round, digits)
  print(y, row.names = FALSE)
  invisible(x)
}
