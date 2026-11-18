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
#' @param purify Iteratively purify the matching score (MH, and LR with
#'   `match = "total"`). Items flagged for DIF are removed from the total
#'   score and the analysis is repeated until the set of flagged items no
#'   longer changes (Holland & Thayer, 1988; Clauser & Mazor, 1998). The
#'   studied item is always kept in its own matching score. Not available
#'   with `match = "theta"`; use [dif_irt()] for purified IRT-based DIF.
#' @param alpha Significance level used to flag items during purification.
#' @param max_iter Maximum number of purification rounds.
#' @return A data frame of class `irt_dif`, one row per item. With
#'   `purify = TRUE`, the column `anchor` shows which items were in the final
#'   matching score, and `attr(, "purification")` holds the number of rounds,
#'   whether the flagged set converged, and the flagged items at each round.
#' @references
#' Holland, P. W., & Thayer, D. T. (1988). Differential item performance and
#' the Mantel-Haenszel procedure. In H. Wainer & H. I. Braun (Eds.),
#' *Test validity* (pp. 129-145). Erlbaum.
#'
#' Clauser, B. E., & Mazor, K. M. (1998). Using statistical procedures to
#' identify differentially functioning test items. *Educational Measurement:
#' Issues and Practice, 17*(1), 31-44.
#'
#' Zumbo, B. D. (1999). *A handbook on the theory and methods of differential
#' item functioning (DIF)*. Directorate of Human Resources Research and
#' Evaluation, Department of National Defense.
#' @examples
#' set.seed(4)
#' d <- sim_irt(1000, 12, "rasch")
#' g <- rep(c("ref", "foc"), 500)
#' th <- attr(d, "truth")$theta[, 1]
#' foc <- g == "foc"
#' d$I03[foc] <- rbinom(500, 1, plogis(th[foc] - 1.2))   # harder for focal
#' d$I08[foc] <- rbinom(500, 1, plogis(th[foc] + 1.0))   # easier for focal
#' dif_test(d, g, method = "mh", purify = TRUE)
#' @export
dif_test <- function(fit, group, method = c("lr", "mh"), match = NULL,
                     reference = NULL, p.adjust = "BH", purify = FALSE, alpha = 0.05,
                     max_iter = 10) {
  method <- match.arg(method)
  if (inherits(fit, "irtfit")) { X <- fit$data } else { X <- .prep_data(fit)$X; fit <- NULL }
  if (is.null(colnames(X))) colnames(X) <- paste0("Item", seq_len(ncol(X)))
  if (any(is.na(group))) stop("'group' has missing values; drop those persons first.", call. = FALSE)
  group <- droplevels(as.factor(group))
  if (length(group) != nrow(X)) stop("'group' must have one value per person.", call. = FALSE)
  if (nlevels(group) < 2) stop("'group' must have at least two levels.", call. = FALSE)
  if (is.null(match)) match <- if (!is.null(fit) && method == "lr") "theta" else "total"
  match <- match.arg(match, c("theta", "total"))
  if (method == "mh" && match == "theta")
    stop("Mantel-Haenszel matches on the total score; drop match = 'theta'.", call. = FALSE)
  if (purify && match == "theta")
    stop("Purification needs match = 'total'. For purified IRT-based DIF use dif_irt().", call. = FALSE)
  run <- function(anchor) {
    o <- if (method == "mh") .dif_mh(X, group, reference, anchor) else .dif_lr(X, group, fit, match, anchor)
    o$p_adj <- stats::p.adjust(o$p, p.adjust)
    o
  }
  items <- colnames(X)
  anchor <- rep(TRUE, ncol(X))
  out <- run(anchor)
  pur <- NULL
  if (purify) {
    history <- list(); converged <- FALSE
    flagged <- !is.na(out$p_adj) & out$p_adj < alpha
    history[[1]] <- items[flagged]
    for (iter in seq_len(max_iter)) {
      new_anchor <- !flagged
      if (sum(new_anchor) < 2) {
        warning("Purification left fewer than two anchor items; results use the previous matching score.",
                call. = FALSE)
        break
      }
      if (identical(new_anchor, anchor)) { converged <- TRUE; break }
      anchor <- new_anchor
      out <- run(anchor)
      flagged <- !is.na(out$p_adj) & out$p_adj < alpha
      history[[iter + 1]] <- items[flagged]
    }
    if (!converged)
      warning("Purification did not converge in ", max_iter, " rounds.", call. = FALSE)
    out$anchor <- anchor
    pur <- list(rounds = length(history), converged = converged, alpha = alpha,
                flagged = stats::setNames(history, paste("round", seq_along(history))))
  }
  attr(out, "method") <- method
  attr(out, "match") <- match
  attr(out, "purification") <- pur
  class(out) <- c("irt_dif", "data.frame")
  out
}

.dif_mh <- function(X, group, reference, anchor = rep(TRUE, ncol(X))) {
  if (nlevels(group) != 2) stop("Mantel-Haenszel needs exactly two groups.", call. = FALSE)
  if (any(apply(X, 2, max, na.rm = TRUE) > 1))
    stop("Mantel-Haenszel is implemented for dichotomous items; use method = 'lr'.", call. = FALSE)
  if (is.null(reference)) reference <- levels(group)[1]
  if (!reference %in% levels(group)) stop("'reference' is not a level of 'group'.", call. = FALSE)
  ref <- group == reference
  cc <- stats::complete.cases(X)
  if (!all(cc))
    message("Mantel-Haenszel uses complete cases: ", sum(!cc), " person(s) with missing responses dropped.")
  X <- X[cc, , drop = FALSE]; ref <- ref[cc]
  base <- rowSums(X[, anchor, drop = FALSE])
  res <- lapply(seq_len(ncol(X)), function(j) {
    total <- if (anchor[j]) base else base + X[, j]   # studied item always in its own matching score
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
    if (VA <= 0 || den == 0 || num == 0) {
      alpha <- NA_real_; chi <- NA_real_; delta <- NA_real_; p <- NA_real_
    } else {
      alpha <- num / den
      chi <- max(abs(A - EA) - 0.5, 0)^2 / VA
      delta <- -2.35 * log(alpha)
      p <- stats::pchisq(chi, 1, lower.tail = FALSE)
    }
    cls <- if (is.na(p)) NA_character_ else if (abs(delta) >= 1.5 && p < 0.05) "C" else
      if (abs(delta) >= 1 && p < 0.05) "B" else "A"
    data.frame(item = colnames(X)[j], alpha_MH = alpha, delta_MH = delta, statistic = chi,
               df = 1, p = p, class = cls)
  })
  do.call(rbind, res)
}

.nagelkerke <- function(ll, ll0, n) {
  cs <- 1 - exp(2 * (ll0 - ll) / n)
  cs / (1 - exp(2 * ll0 / n))
}

.dif_lr <- function(X, group, fit, match, anchor = rep(TRUE, ncol(X))) {
  if (match == "theta") {
    if (is.null(fit)) stop("match = 'theta' requires a fitted model.", call. = FALSE)
    if (fit$spec$D > 1) stop("Use match = 'total' for multidimensional models.", call. = FALSE)
    m_theta <- as.vector(fit$posterior %*% fit$theta)
  }
  base <- rowSums(X[, anchor, drop = FALSE], na.rm = TRUE)
  res <- lapply(seq_len(ncol(X)), function(j) {
    m <- if (match == "theta") m_theta else {
      tot <- if (anchor[j]) base else base + ifelse(is.na(X[, j]), 0, X[, j])
      as.vector(scale(tot))
    }
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
  lab <- c(mh = "Mantel-Haenszel", lr = "logistic regression")[[attr(x, "method")]]
  cat("DIF analysis (", lab, ", matched on ", attr(x, "match"), "): ",
      sum(x$class %in% c("B", "C")), " item(s) with B/C-level DIF\n", sep = "")
  pur <- attr(x, "purification")
  if (!is.null(pur))
    cat(sprintf("Purified matching score: %d anchor item(s) after %d round(s)%s\n",
                sum(x$anchor), pur$rounds, if (pur$converged) "" else " (did not converge)"))
  y <- x; class(y) <- "data.frame"
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], round, digits)
  print(y, row.names = FALSE)
  invisible(x)
}
