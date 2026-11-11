# Item banks ------------------------------------------------------------------------
# A light, model-agnostic representation of calibrated unidimensional items, used by
# linking, IRT-based DIF, summed-score tables and computerized adaptive testing.

#' Build an item bank from a fitted model
#'
#' Extracts the calibrated items of a unidimensional [irt_fit()] model into an
#' `irt_bank`, the object used by [cat_design()], [irt_link()] and friends.
#'
#' @param fit A unidimensional `irtfit` object.
#' @param content Optional vector of content areas (one per item) for content
#'   balancing in adaptive tests.
#' @return An object of class `irt_bank`.
#' @seealso [make_bank()] to build a bank from known parameters.
#' @examples
#' fit <- irt_fit(sim_irt(500, 20, "2pl"), "2pl", se = "none")
#' bank <- item_bank(fit)
#' bank
#' @export
item_bank <- function(fit, content = NULL) {
  if (!inherits(fit, "irtfit")) stop("'fit' must be an irtfit object.", call. = FALSE)
  if (fit$spec$D != 1) stop("Item banks require a unidimensional model.", call. = FALSE)
  items <- lapply(seq_len(ncol(fit$data)), function(j) {
    it <- fit$spec$items[[j]]
    pp <- .decode_item(it, fit$par, 1)
    list(name = colnames(fit$data)[j], model = it$model, fam = it$fam, a = pp$a[1],
         d = pp$d, g = pp$g, u = if (is.null(pp$u)) 1 else pp$u, K = it$K, scores = fit$labels[[j]])
  })
  .new_bank(items, content, mu = 0, sigma = sqrt(fit$Sigma[1, 1]))
}

#' Build an item bank from item parameters
#'
#' Creates an `irt_bank` from known item parameters in the usual IRT metric,
#' for example parameters published with an operational test, so that
#' adaptive tests and scoring tables can be run without re-estimating.
#'
#' @param b Item difficulties: a vector (dichotomous items) or a matrix with
#'   one row per item and one column per threshold (polytomous items; pad with
#'   `NA` when items have different numbers of categories). For `"grm"` these
#'   are the boundary locations, for `"pcm"`/`"gpcm"` the step difficulties.
#' @param a Discriminations (default 1).
#' @param g Lower asymptotes for 3PL/4PL items (default 0).
#' @param u Upper asymptotes for 4PL items (default 1).
#' @param model One model per item or one for all: `"rasch"`, `"2pl"`,
#'   `"3pl"`, `"4pl"`, `"grm"`, `"pcm"` or `"gpcm"`.
#' @param names Item names.
#' @param content Optional content areas.
#' @param D Scaling constant; use 1.702 for parameters reported in the
#'   normal-ogive metric.
#' @return An object of class `irt_bank`.
#' @examples
#' bank <- make_bank(b = rnorm(100), a = rlnorm(100, 0, 0.3), model = "2pl")
#' bank
#' @export
make_bank <- function(b, a = 1, g = 0, u = 1, model = "2pl", names = NULL, content = NULL, D = 1) {
  b <- as.matrix(b)
  J <- nrow(b)
  model <- tolower(rep(model, length.out = J))
  bad <- setdiff(model, c("rasch", "2pl", "3pl", "4pl", "grm", "pcm", "gpcm"))
  if (length(bad)) stop("Unsupported model(s) in make_bank(): ", paste(bad, collapse = ", "), call. = FALSE)
  a <- rep(a, length.out = J) * D
  g <- rep(g, length.out = J); u <- rep(u, length.out = J)
  if (is.null(names)) names <- sprintf("I%03d", seq_len(J))
  items <- lapply(seq_len(J), function(j) {
    bj <- b[j, !is.na(b[j, ])]
    m <- model[j]
    if (m == "rasch" && length(bj) > 1) m <- "pcm"
    aj <- if (m %in% c("rasch", "pcm")) 1 * D else a[j]
    if (m %in% c("rasch", "2pl", "3pl", "4pl") && length(bj) != 1)
      stop("Dichotomous item ", names[j], " needs exactly one difficulty.", call. = FALSE)
    list(name = names[j], model = m, fam = .family_of(m), a = aj, d = -aj * bj,
         g = if (m %in% c("3pl", "4pl")) g[j] else 0, u = if (m == "4pl") u[j] else 1,
         K = length(bj) + 1L, scores = 0:length(bj))
  })
  .new_bank(items, content, mu = 0, sigma = 1)
}

.new_bank <- function(items, content, mu, sigma) {
  if (!is.null(content)) {
    if (length(content) != length(items)) stop("'content' needs one value per item.", call. = FALSE)
    for (j in seq_along(items)) items[[j]]$content <- as.character(content[j])
  }
  names(items) <- vapply(items, `[[`, character(1), "name")
  structure(list(items = items, mu = mu, sigma = sigma), class = "irt_bank")
}

#' @export
print.irt_bank <- function(x, ...) {
  m <- table(vapply(x$items, `[[`, character(1), "model"))
  cat("<irt_bank> ", length(x$items), " items: ",
      paste(sprintf("%s (%d)", names(m), m), collapse = ", "), "\n", sep = "")
  ct <- vapply(x$items, function(it) if (is.null(it$content)) NA_character_ else it$content, character(1))
  if (!all(is.na(ct))) cat("Content areas:", paste(sprintf("%s (%d)", names(table(ct)), table(ct)), collapse = ", "), "\n")
  invisible(x)
}

#' @export
length.irt_bank <- function(x) length(x$items)

#' @export
`[.irt_bank` <- function(x, i) {
  x$items <- x$items[i]
  x
}

#' @export
as.data.frame.irt_bank <- function(x, ...) {
  maxk <- max(vapply(x$items, `[[`, numeric(1), "K"))
  rows <- lapply(x$items, function(it) {
    b <- rep(NA_real_, maxk - 1)
    bb <- -it$d / it$a
    b[seq_along(bb)] <- bb
    out <- data.frame(item = it$name, model = it$model, a = it$a,
                      b = if (it$K == 2) bb else mean(bb), g = it$g, u = if (is.null(it$u)) 1 else it$u,
                      stringsAsFactors = FALSE)
    if (maxk > 2) out <- cbind(out, stats::setNames(as.list(b), paste0("b", seq_len(maxk - 1))))
    if (!is.null(it$content)) out$content <- it$content
    out
  })
  do.call(rbind, rows)
}

# evaluation helpers --------------------------------------------------------------
.bp <- function(it, theta) .item_probs(it$fam, list(a = it$a, d = it$d, g = it$g, u = it$u, K = it$K), matrix(theta))

.bi <- function(it, theta, h = 1e-4) {
  P <- .bp(it, theta)
  dP <- (.bp(it, theta + h) - .bp(it, theta - h)) / (2 * h)
  rowSums(dP^2 / P)
}

.be <- function(it, theta) as.vector(.bp(it, theta) %*% it$scores)

# theta* = A * theta + B
.bank_transform <- function(bank, A, B) {
  bank$items <- lapply(bank$items, function(it) {
    it$d <- it$d - (it$a / A) * B
    it$a <- it$a / A
    it
  })
  bank$mu <- A * bank$mu + B
  bank$sigma <- A * bank$sigma
  bank
}

# Lord-Wingersky recursion: Q x (max summed score + 1) matrix of P(score = s | theta)
.lord_wingersky <- function(items, theta) {
  f <- matrix(1, length(theta), 1)
  for (it in items) {
    P <- .bp(it, theta)
    K <- it$K
    nf <- matrix(0, length(theta), ncol(f) + K - 1)
    for (k in seq_len(K)) nf[, k - 1 + seq_len(ncol(f))] <- nf[, k - 1 + seq_len(ncol(f))] + f * P[, k]
    f <- nf
  }
  f
}

.quad <- function(mu = 0, sigma = 1, n = 81) {
  th <- seq(mu - 6 * sigma, mu + 6 * sigma, length.out = n)
  w <- stats::dnorm(th, mu, sigma); list(theta = th, w = w / sum(w))
}

#' Summed-score conversion table
#'
#' Uses the Lord-Wingersky recursion to give, for every possible summed score,
#' the EAP trait estimate, its posterior SD and the model-implied proportion of
#' examinees obtaining that score. This is the usual way to report IRT scores
#' for tests scored by number-correct.
#'
#' @param x An `irtfit` (unidimensional) or `irt_bank`.
#' @return A data frame with columns `score`, `theta`, `se` and `expected_prop`.
#' @examples
#' fit <- irt_fit(sim_irt(500, 10, "2pl"), "2pl", se = "none")
#' score_table(fit)
#' @export
score_table <- function(x) {
  bank <- if (inherits(x, "irtfit")) item_bank(x) else x
  q <- .quad(bank$mu, bank$sigma, 121)
  f <- .lord_wingersky(bank$items, q$theta)
  marg <- as.vector(q$w %*% f)
  post <- sweep(f * q$w, 2, marg, "/")
  eap <- as.vector(q$theta %*% post)
  sd <- sqrt(pmax(as.vector(q$theta^2 %*% post) - eap^2, 0))
  offset <- sum(vapply(bank$items, function(it) min(it$scores), numeric(1)))
  data.frame(score = offset + seq_len(ncol(f)) - 1, theta = eap, se = sd, expected_prop = marg)
}
