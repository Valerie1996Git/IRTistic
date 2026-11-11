#' Estimate person measures
#'
#' @param fit An `irtfit` object.
#' @param method `"EAP"` (expected a posteriori, default), `"MAP"` (maximum a
#'   posteriori), `"ML"` (maximum likelihood) or `"WLE"` (Warm's weighted
#'   likelihood).
#' @param bounds Range searched for ML/WLE/MAP estimates.
#' @return A data frame with one row per person: the measure(s), their
#'   standard error(s), raw score and number of responses. For ML, persons
#'   with extreme (zero or perfect) scores are flagged in `extreme` and
#'   placed at the bounds.
#' @examples
#' dat <- sim_irt(200, 10, "rasch")
#' fit <- irt_fit(dat, "rasch")
#' head(irt_scores(fit, "WLE"))
#' @export
irt_scores <- function(fit, method = c("EAP", "MAP", "ML", "WLE"), bounds = c(-6, 6)) {
  method <- toupper(method)
  method <- match.arg(method)
  X <- fit$data; D <- fit$spec$D; N <- nrow(X)
  dn <- colnames(fit$Q)
  raw <- rowSums(X, na.rm = TRUE)
  nresp <- rowSums(!is.na(X))
  maxraw <- as.vector((!is.na(X)) %*% vapply(fit$spec$items, function(it) it$K - 1, numeric(1)))
  extreme <- raw == 0 | raw == maxraw

  if (method == "EAP") {
    est <- fit$posterior %*% fit$theta
    sd <- sqrt(pmax(fit$posterior %*% fit$theta^2 - est^2, 0))
  } else if (D == 1) {
    res <- .score_uni(fit, method, bounds)
    est <- matrix(res$est); sd <- matrix(res$se)
  } else {
    res <- .score_multi(fit, method, bounds)
    est <- res$est; sd <- res$se
  }
  out <- data.frame(person = if (!is.null(rownames(X))) rownames(X) else seq_len(N))
  for (d in seq_len(D)) {
    out[[if (D == 1) "theta" else dn[d]]] <- est[, d]
    out[[if (D == 1) "se" else paste0("se_", dn[d])]] <- sd[, d]
  }
  out$raw <- raw; out$n <- nresp
  if (method %in% c("ML", "WLE")) out$extreme <- extreme
  attr(out, "method") <- method
  out
}

# vectorised per-person log-likelihood and information (unidimensional)
.loglik_at <- function(fit, th) {
  X <- fit$data; ll <- numeric(length(th)); tm <- matrix(th)
  for (j in seq_len(ncol(X))) {
    it <- fit$spec$items[[j]]
    P <- .item_probs(it$fam, .decode_item(it, fit$par, 1), tm)
    ok <- !is.na(X[, j])
    ll[ok] <- ll[ok] + log(P[cbind(which(ok), X[ok, j] + 1L)])
  }
  ll
}

.info_at <- function(fit, th, h = 1e-4) {
  X <- fit$data; info <- numeric(length(th))
  for (j in seq_len(ncol(X))) {
    it <- fit$spec$items[[j]]
    pp <- .decode_item(it, fit$par, 1)
    P <- .item_probs(it$fam, pp, matrix(th))
    dP <- (.item_probs(it$fam, pp, matrix(th + h)) - .item_probs(it$fam, pp, matrix(th - h))) / (2 * h)
    ok <- !is.na(X[, j])
    info[ok] <- info[ok] + rowSums(dP^2 / P)[ok]
  }
  info
}

.score_uni <- function(fit, method, bounds) {
  s2 <- fit$Sigma[1, 1]
  obj <- function(th) {
    v <- .loglik_at(fit, th)
    if (method == "MAP") v <- v - th^2 / (2 * s2)
    if (method == "WLE") v <- v + 0.5 * log(.info_at(fit, th))
    v
  }
  grid <- seq(bounds[1], bounds[2], length.out = 121)
  vals <- vapply(grid, function(g) obj(rep(g, nrow(fit$data))), numeric(nrow(fit$data)))
  th <- grid[max.col(vals, ties.method = "first")]
  h <- 1e-3
  for (k in 1:30) {
    f0 <- obj(th); fp <- obj(th + h); fm <- obj(th - h)
    g1 <- (fp - fm) / (2 * h); g2 <- (fp - 2 * f0 + fm) / h^2
    step <- ifelse(g2 < 0, -g1 / g2, sign(g1) * 0.1)
    step <- pmax(pmin(step, 0.5), -0.5)
    th <- pmin(pmax(th + step, bounds[1]), bounds[2])
    if (max(abs(step)) < 1e-6) break
  }
  info <- .info_at(fit, th)
  if (method == "MAP") info <- info + 1 / s2
  list(est = th, se = 1 / sqrt(info))
}

.score_multi <- function(fit, method, bounds) {
  X <- fit$data; D <- fit$spec$D; N <- nrow(X)
  Si <- solve(fit$Sigma)
  items <- lapply(seq_len(ncol(X)), function(j) {
    it <- fit$spec$items[[j]]; list(fam = it$fam, p = .decode_item(it, fit$par, D)) })
  infomat <- function(th, x) {
    I <- matrix(0, D, D); tm <- matrix(th, 1)
    for (j in which(!is.na(x))) {
      P <- .item_probs(items[[j]]$fam, items[[j]]$p, tm)
      a <- items[[j]]$p$a
      gr <- vapply(seq_len(D), function(d) {
        e <- numeric(D); e[d] <- 1e-4
        (.item_probs(items[[j]]$fam, items[[j]]$p, tm + e) -
           .item_probs(items[[j]]$fam, items[[j]]$p, tm - e)) / 2e-4 }, numeric(length(P)))
      gr <- matrix(gr, ncol = D)
      I <- I + crossprod(gr / sqrt(as.vector(P)))
    }
    I
  }
  est <- matrix(0, N, D); se <- matrix(NA, N, D)
  start <- fit$posterior %*% fit$theta
  for (i in seq_len(N)) {
    x <- X[i, ]
    obj <- function(th) {
      tm <- matrix(th, 1); v <- 0
      for (j in which(!is.na(x)))
        v <- v + log(.item_probs(items[[j]]$fam, items[[j]]$p, tm)[x[j] + 1])
      if (method == "MAP") v <- v - 0.5 * sum((th %*% Si) * th)
      if (method == "WLE") v <- v + 0.5 * log(max(det(infomat(th, x)), 1e-300))
      -v
    }
    o <- stats::optim(start[i, ], obj, method = "L-BFGS-B", lower = bounds[1], upper = bounds[2])
    est[i, ] <- o$par
    I <- infomat(o$par, x); if (method == "MAP") I <- I + Si
    se[i, ] <- tryCatch(sqrt(diag(solve(I))), error = function(e) rep(NA, D))
  }
  list(est = est, se = se)
}
