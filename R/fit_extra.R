#' Orlando-Thissen S-X2 item fit
#'
#' Compares, within each summed-score group, the observed proportion of
#' responses in each category with the model-implied proportion computed by
#' the Lord-Wingersky recursion (Orlando & Thissen, 2000; Kang & Chen, 2008
#' for polytomous items). Unlike fit statistics that condition on estimated
#' traits, S-X2 conditions on an observed quantity and so has a well-behaved
#' chi-square reference distribution.
#'
#' @param fit A unidimensional `irtfit`. Persons with missing responses are
#'   dropped (S-X2 requires complete data).
#' @param min_expected Adjacent score groups are merged until every expected
#'   cell count is at least this value.
#' @param p.adjust Multiple-comparison adjustment.
#' @return A data frame with `S_X2`, `df`, `p`, adjusted p and RMSEA per item.
#' @references
#' Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for
#' dichotomous item response theory models. *Applied Psychological
#' Measurement, 24*(1), 50-64.
#'
#' Kang, T., & Chen, T. T. (2008). Performance of the generalized S-X2 item
#' fit index for polytomous IRT models. *Journal of Educational Measurement,
#' 45*(4), 391-406.
#' @examples
#' fit <- irt_fit(sim_irt(800, 10, "2pl"), "2pl", se = "none")
#' sx2(fit)
#' @export
sx2 <- function(fit, min_expected = 1, p.adjust = "BH") {
  if (fit$spec$D != 1) stop("S-X2 is implemented for unidimensional models.", call. = FALSE)
  X <- fit$data
  cc <- stats::complete.cases(X)
  if (!all(cc)) message("Dropping ", sum(!cc), " person(s) with missing responses.")
  X <- X[cc, , drop = FALSE]
  bank <- item_bank(fit)
  q <- .quad(0, sqrt(fit$Sigma[1, 1]), 81)
  items <- bank$items
  total <- rowSums(X)
  fall <- .lord_wingersky(items, q$theta)
  Ns <- tabulate(total + 1, ncol(fall))
  res <- lapply(seq_along(items), function(j) {
    it <- items[[j]]; K <- it$K
    fminus <- .lord_wingersky(items[-j], q$theta)
    P <- .bp(it, q$theta)
    S <- ncol(fall) - 1
    E <- matrix(0, S + 1, K); O <- matrix(0, S + 1, K)
    denom <- as.vector(q$w %*% fall)
    for (k in seq_len(K)) {
      sidx <- (k - 1) + seq_len(ncol(fminus))           # score s = (k-1) + r
      num <- as.vector(q$w %*% (fminus * P[, k]))
      E[sidx, k] <- num / denom[sidx]
    }
    for (s in 0:S) if (Ns[s + 1]) O[s + 1, ] <- tabulate(X[total == s, j] + 1, K)
    E <- E * Ns
    # drop score groups with a single feasible category (no information)
    keep <- rowSums(E > 1e-10) > 1 & Ns > 0
    E <- E[keep, , drop = FALSE]; O <- O[keep, , drop = FALSE]
    # merge adjacent groups until expected counts are large enough
    grp <- integer(nrow(E)); g <- 1; acc <- rep(0, K)
    for (r in seq_len(nrow(E))) {
      grp[r] <- g; acc <- acc + E[r, ]
      if (all(acc[acc > 1e-10] >= min_expected)) { g <- g + 1; acc <- rep(0, K) }
    }
    if (any(acc > 0) && g > 1) grp[grp == g] <- g - 1
    Eg <- rowsum(E, grp); Og <- rowsum(O, grp)
    ok <- Eg > 1e-10
    chi <- sum((Og[ok] - Eg[ok])^2 / Eg[ok])
    npar <- length(fit$spec$items[[j]]$idx) - length(fit$spec$items[[j]]$idx_tau)
    df <- sum(rowSums(ok) - 1) - npar
    data.frame(item = it$name, S_X2 = chi, df = df,
               p = if (df > 0) stats::pchisq(chi, df, lower.tail = FALSE) else NA_real_,
               RMSEA = if (df > 0) sqrt(max(chi - df, 0) / (df * (nrow(X) - 1))) else NA_real_)
  })
  out <- do.call(rbind, res)
  out$p_adj <- stats::p.adjust(out$p, p.adjust)
  class(out) <- c("irt_sx2", "data.frame")
  out
}

#' Classification accuracy and consistency
#'
#' Estimates how accurately and how consistently examinees are classified
#' into performance levels.
#'
#' * `method = "lee"` (default) uses the IRT model and the Lord-Wingersky
#'   recursion to obtain, for every trait value, the distribution of summed
#'   scores and hence of observed classifications (Lee, 2010). It needs
#'   cut scores on the summed-score scale (`score_cuts`) or on the trait
#'   scale (`theta_cuts`, converted through the test characteristic curve).
#' * `method = "rudner"` uses each person's trait estimate and standard
#'   error with cuts on the trait scale (Rudner, 2001, 2005).
#'
#' @param fit A unidimensional `irtfit` (or `irt_bank` for `"lee"`).
#' @param score_cuts Summed-score cuts: a person is placed at or above a
#'   level when the score is at least the cut.
#' @param theta_cuts Cuts on the trait scale.
#' @param method `"lee"` or `"rudner"`.
#' @param person_method Trait estimator for `"rudner"`.
#' @return A list of class `irt_classification` with overall accuracy,
#'   consistency, Cohen's kappa, the false positive/negative rates at each
#'   cut, and the joint true x observed classification table.
#' @references
#' Lee, W.-C. (2010). Classification consistency and accuracy for complex
#' assessments using item response theory. *Journal of Educational
#' Measurement, 47*(1), 1-17.
#'
#' Rudner, L. M. (2005). Expected classification accuracy. *Practical
#' Assessment, Research & Evaluation, 10*(13).
#' @examples
#' fit <- irt_fit(sim_irt(800, 20, "2pl"), "2pl", se = "none")
#' classification_accuracy(fit, theta_cuts = c(-0.5, 1))
#' @export
classification_accuracy <- function(fit, score_cuts = NULL, theta_cuts = NULL,
                                    method = c("lee", "rudner"), person_method = "EAP") {
  method <- match.arg(method)
  if (method == "rudner") {
    if (!inherits(fit, "irtfit")) stop("Rudner's method needs a fitted model.", call. = FALSE)
    if (is.null(theta_cuts)) stop("Rudner's method needs 'theta_cuts'.", call. = FALSE)
    sc <- irt_scores(fit, person_method)
    cuts <- c(-Inf, sort(theta_cuts), Inf); L <- length(cuts) - 1
    pr <- vapply(seq_len(L), function(c) stats::pnorm((cuts[c + 1] - sc$theta) / sc$se) -
                   stats::pnorm((cuts[c] - sc$theta) / sc$se), numeric(nrow(sc)))
    obs <- findInterval(sc$theta, sort(theta_cuts)) + 1
    # rows = true level, columns = observed level
    joint <- vapply(seq_len(L), function(o) colSums(pr[obs == o, , drop = FALSE]), numeric(L)) / nrow(sc)
    cons <- mean(rowSums(pr^2))
    chance <- sum(colMeans(pr)^2)
    return(.class_out(joint, cons, chance, method, theta_cuts))
  }
  bank <- if (inherits(fit, "irtfit")) item_bank(fit) else fit
  q <- .quad(bank$mu, bank$sigma, 121)
  f <- .lord_wingersky(bank$items, q$theta)            # P(score | theta)
  offset <- sum(vapply(bank$items, function(it) min(it$scores), numeric(1)))
  scores <- offset + seq_len(ncol(f)) - 1
  tcc <- Reduce(`+`, lapply(bank$items, function(it) .be(it, q$theta)))
  if (is.null(score_cuts)) {
    if (is.null(theta_cuts)) stop("Supply 'score_cuts' or 'theta_cuts'.", call. = FALSE)
    score_cuts <- vapply(sort(theta_cuts), function(t)
      ceiling(sum(vapply(bank$items, function(it) .be(it, t), numeric(1))) - 1e-9), numeric(1))
    true_level <- findInterval(q$theta, sort(theta_cuts)) + 1
  } else true_level <- findInterval(tcc, sort(score_cuts)) + 1
  score_cuts <- sort(score_cuts)
  L <- length(score_cuts) + 1
  obs_level <- findInterval(scores, score_cuts) + 1
  pc <- vapply(seq_len(L), function(c) rowSums(f[, obs_level == c, drop = FALSE]), numeric(nrow(f)))
  joint <- t(vapply(seq_len(L), function(tl) colSums(pc[true_level == tl, , drop = FALSE] * q$w[true_level == tl]),
                    numeric(L)))
  cons <- sum(q$w * rowSums(pc^2))
  chance <- sum(colSums(pc * q$w)^2)
  out <- .class_out(joint, cons, chance, method, if (is.null(theta_cuts)) NULL else sort(theta_cuts))
  out$score_cuts <- score_cuts
  out
}

.class_out <- function(joint, cons, chance, method, theta_cuts) {
  L <- nrow(joint)
  dimnames(joint) <- list(true = paste0("L", seq_len(L)), observed = paste0("L", seq_len(L)))
  acc <- sum(diag(joint))
  cuts <- data.frame(cut = seq_len(L - 1),
                     accuracy = vapply(seq_len(L - 1), function(k) {
                       below <- seq_len(k); sum(joint[below, below]) + sum(joint[-below, -below]) }, numeric(1)),
                     false_positive = vapply(seq_len(L - 1), function(k) sum(joint[seq_len(k), -seq_len(k)]), numeric(1)),
                     false_negative = vapply(seq_len(L - 1), function(k) sum(joint[-seq_len(k), seq_len(k)]), numeric(1)))
  structure(list(accuracy = acc, consistency = cons, kappa = (cons - chance) / (1 - chance),
                 by_cut = cuts, joint = joint, method = method, theta_cuts = theta_cuts),
            class = "irt_classification")
}

#' @export
print.irt_classification <- function(x, digits = 3, ...) {
  cat(sprintf("Classification accuracy and consistency (%s method)\n", if (x$method == "lee") "Lee" else "Rudner"))
  if (!is.null(x$score_cuts)) cat("  summed-score cuts:", paste(x$score_cuts, collapse = ", "), "\n")
  if (!is.null(x$theta_cuts)) cat("  theta cuts:", paste(round(x$theta_cuts, 3), collapse = ", "), "\n")
  cat(sprintf("  accuracy = %.3f | consistency = %.3f | kappa = %.3f\n\n", x$accuracy, x$consistency, x$kappa))
  print(x$by_cut, digits = digits, row.names = FALSE)
  cat("\nJoint distribution (rows = true level, columns = observed level):\n")
  print(round(x$joint, digits))
  invisible(x)
}
