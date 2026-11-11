#' Simulate item response data
#'
#' @param n Number of persons.
#' @param nitems Number of items (ignored if `b` is given).
#' @param model One model for all items or one per item (see [irt_fit()]).
#' @param ncat Number of categories for polytomous items.
#' @param a Slopes: a vector (unidimensional) or a matrix items x dimensions.
#'   Drawn from a log-normal(0, 0.25) if `NULL`.
#' @param b Difficulties: a vector for dichotomous items or a matrix
#'   items x (ncat - 1) of thresholds. Drawn if `NULL`.
#' @param g Guessing parameters for 3PL/4PL items.
#' @param u Upper asymptotes for 4PL items.
#' @param Q Optional loading matrix for multidimensional data.
#' @param Sigma Latent covariance matrix.
#' @param theta Optional matrix of true abilities.
#' @return A data frame of responses with attribute `"truth"` holding the
#'   generating parameters and abilities.
#' @examples
#' dat <- sim_irt(300, nitems = 8, model = "grm", ncat = 4)
#' head(dat)
#' @export
sim_irt <- function(n, nitems = 10, model = "2pl", ncat = 4, a = NULL, b = NULL,
                    g = NULL, u = NULL, Q = NULL, Sigma = NULL, theta = NULL) {
  model <- tolower(model)
  if (!is.null(b)) nitems <- NROW(b)
  if (length(model) == 1) model <- rep(model, nitems)
  if (is.null(Q)) Q <- matrix(1, nitems, 1)
  D <- ncol(Q)
  if (is.null(Sigma)) Sigma <- diag(D)
  if (is.null(theta)) theta <- MASS::mvrnorm(n, rep(0, D), Sigma)
  theta <- matrix(theta, n, D)
  if (is.null(a)) a <- matrix(stats::rlnorm(nitems * D, 0, 0.25), nitems, D)
  a <- matrix(a, nitems, D) * (Q != 0)
  fixed <- model %in% c("rasch", "pcm", "rsm")
  a[fixed, ] <- (Q[fixed, , drop = FALSE] != 0) * 1
  K <- ifelse(model %in% .dich_models, 2, ncat)
  if (is.null(b)) {
    b <- matrix(NA, nitems, max(K) - 1)
    loc <- stats::rnorm(nitems)
    steps <- sort(stats::rnorm(max(K) - 1, 0, 0.8))
    steps <- steps - mean(steps)
    for (j in seq_len(nitems)) {
      if (K[j] == 2) b[j, 1] <- loc[j] else
        b[j, seq_len(K[j] - 1)] <- loc[j] + if (model[j] == "rsm") steps else
          sort(stats::rnorm(K[j] - 1, 0, 0.8))
    }
  }
  b <- matrix(b, nitems)
  if (is.null(g)) g <- ifelse(model %in% c("3pl", "4pl"), stats::runif(nitems, 0.1, 0.25), 0)
  if (is.null(u)) u <- ifelse(model == "4pl", stats::runif(nitems, 0.85, 0.97), 1)
  u <- rep(u, length.out = nitems)
  g <- rep(g, length.out = nitems)
  X <- matrix(NA_integer_, n, nitems)
  for (j in seq_len(nitems)) {
    fam <- .family_of(model[j])
    bj <- b[j, !is.na(b[j, ])]
    aj <- a[j, ]
    slope <- if (D == 1) aj[1] else 1
    d <- -bj * slope
    if (D > 1) d <- -bj
    P <- .item_probs(fam, list(a = aj, d = d, g = g[j], u = u[j], K = K[j]), theta)
    rnd <- stats::runif(n)
    X[, j] <- rowSums(rnd > t(apply(P, 1, cumsum))[, -K[j], drop = FALSE])
  }
  colnames(X) <- sprintf("I%02d", seq_len(nitems))
  out <- as.data.frame(X)
  attr(out, "truth") <- list(theta = theta, a = a, b = b, g = g, u = u, model = model, Q = Q)
  out
}
