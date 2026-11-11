#' Fit an item response theory model
#'
#' Estimates unidimensional or confirmatory multidimensional IRT models by
#' marginal maximum likelihood (MML) using an EM algorithm with fixed
#' quadrature. Dichotomous and polytomous items can be mixed in one test.
#'
#' @param data A matrix or data frame of item responses, persons in rows and
#'   items in columns. Responses may be any ordered numeric codes; they are
#'   recoded to consecutive integers starting at 0. `NA` is treated as missing
#'   (not presented) and is ignored in the likelihood.
#' @param model Character, either one model for all items or one per item.
#'   One of `"rasch"`, `"2pl"`, `"3pl"`, `"4pl"`, `"grm"` (graded response),
#'   `"pcm"` (partial credit), `"gpcm"` (generalized partial credit) or
#'   `"rsm"` (rating scale). `"rasch"` applied to a polytomous item is
#'   treated as `"pcm"`.
#' @param Q Optional item-by-dimension loading matrix (0/1) for confirmatory
#'   multidimensional models. Between- and within-item multidimensionality
#'   are both allowed. Defaults to a single dimension.
#' @param quadpts Number of quadrature points per dimension. Defaults to 61
#'   for one dimension, 21 for two, 11 for three and 7 otherwise.
#' @param priors Optional list of priors used for stabilisation:
#'   `g = c(alpha, beta)` a Beta prior on 3PL guessing (default `c(5, 17)`
#'   when any 3PL item is present) and `a = c(meanlog, sdlog)` a log-normal
#'   prior on slopes (default none).
#' @param se Standard error method: `"hessian"` (numerical observed
#'   information, default), `"crossprod"` (outer product of person gradients)
#'   or `"none"`.
#' @param maxit Maximum number of EM cycles.
#' @param tol Convergence tolerance on the absolute change in log-likelihood
#'   and on the maximum parameter change.
#' @param verbose Print progress every 25 cycles.
#'
#' @return An object of class `irtfit`. Use [summary()], [coef()],
#'   [irt_scores()], [item_fit()], [plot()] and friends on it.
#'
#' @details The latent distribution is multivariate normal with mean zero.
#'   Variances of dimensions measured only by fixed-slope (Rasch family)
#'   items are estimated; otherwise variances are fixed at 1 to identify
#'   the slopes. Correlations between dimensions are always estimated.
#'   Standard errors treat the latent covariance matrix as known.
#'
#' @examples
#' set.seed(1)
#' dat <- sim_irt(500, nitems = 10, model = "2pl")
#' fit <- irt_fit(dat, model = "2pl")
#' summary(fit)
#' @export
irt_fit <- function(data, model = "rasch", Q = NULL, quadpts = NULL, priors = NULL,
                    se = c("hessian", "crossprod", "none"), maxit = 500,
                    tol = 1e-4, verbose = FALSE) {
  cl <- match.call()
  se <- match.arg(se)
  pd <- .prep_data(data)
  X <- pd$X; K <- pd$K; J <- ncol(X); N <- nrow(X)
  model <- tolower(model)
  if (length(model) == 1) model <- rep(model, J)
  if (length(model) != J) stop("'model' must have length 1 or ncol(data).", call. = FALSE)
  bad <- setdiff(model, .all_models)
  if (length(bad)) stop("Unknown model(s): ", paste(bad, collapse = ", "), call. = FALSE)
  promote <- model == "rasch" & K > 2
  if (any(promote)) model[promote] <- "pcm"
  dich_bad <- model %in% c("2pl", "3pl", "4pl") & K > 2
  if (any(dich_bad))
    stop("Items ", paste(colnames(X)[dich_bad], collapse = ", "),
         " are polytomous; use 'grm', 'gpcm' or 'pcm' for them.", call. = FALSE)
  if (any(model == "rsm") && !all(model == "rsm"))
    stop("The rating scale model must be used for all items.", call. = FALSE)

  if (is.null(Q)) Q <- matrix(1, J, 1)
  Q <- as.matrix(Q)
  if (nrow(Q) != J) stop("'Q' must have one row per item.", call. = FALSE)
  if (any(rowSums(Q != 0) == 0)) stop("Every item must load on a dimension in 'Q'.", call. = FALSE)
  D <- ncol(Q)
  if (is.null(colnames(Q))) colnames(Q) <- if (D == 1) "theta" else paste0("F", seq_len(D))
  rownames(Q) <- colnames(X)

  if (is.null(priors)) priors <- list()
  if (any(model %in% c("3pl", "4pl")) && !"g" %in% names(priors)) priors$g <- c(5, 17)
  if (any(model == "4pl") && !"u" %in% names(priors)) priors$u <- c(17, 2)

  spec <- .build_spec(model, K, Q, priors)

  # quadrature grid
  if (is.null(quadpts)) quadpts <- c(61, 21, 11, 7)[min(D, 4)]
  if (D > 4) warning("More than 4 dimensions with fixed quadrature will be slow.", call. = FALSE)
  rng <- if (D == 1) 6 else 5
  nodes <- seq(-rng, rng, length.out = quadpts)
  theta <- as.matrix(expand.grid(rep(list(nodes), D)))
  colnames(theta) <- colnames(Q)

  # which dimension variances are free (only fixed-slope items load on it)
  var_free <- vapply(seq_len(D), function(dd) {
    on <- Q[, dd] != 0
    !any(vapply(spec$items[on], function(it) length(it$idx_a) > 0, logical(1)))
  }, logical(1))
  Sigma <- diag(D)

  # data indexes
  idx <- lapply(seq_len(J), function(j) { v <- X[, j] + 1L; v[is.na(v)] <- K[j] + 1L; v })
  onehot <- lapply(seq_len(J), function(j) {
    m <- matrix(0, N, K[j]); ok <- !is.na(X[, j])
    m[cbind(which(ok), X[ok, j] + 1L)] <- 1; m })
  ctx <- list(spec = spec, theta = theta, idx = idx, onehot = onehot, N = N, J = J, D = D)

  par <- .start_values(spec, X)
  ll_old <- -Inf
  blocks <- if (any(model == "rsm")) list(seq_len(spec$npar)) else
    lapply(spec$items, function(it) it$idx)
  block_items <- if (any(model == "rsm")) list(seq_len(J)) else as.list(seq_len(J))
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    es <- .estep(ctx, par, Sigma)
    par_old <- par
    for (b in seq_along(blocks)) {
      par[blocks[[b]]] <- .mstep_block(ctx, par, blocks[[b]], block_items[[b]], es$r)
    }
    # latent covariance update
    if (D > 1 || any(var_free)) {
      w <- colSums(es$post) / N
      S <- crossprod(theta * sqrt(w))
      fixed <- !var_free
      if (any(fixed)) {
        sc <- rep(1, D); sc[fixed] <- 1 / sqrt(diag(S)[fixed])
        S <- diag(sc, D) %*% S %*% diag(sc, D)
      }
      Sigma <- (S + t(S)) / 2
    }
    dll <- es$ll - ll_old
    if (verbose && iter %% 25 == 0)
      message(sprintf("Iteration %d: log-likelihood = %.4f", iter, es$ll))
    if (is.finite(dll) && abs(dll) < tol && max(abs(par - par_old)) < tol * 10) {
      converged <- TRUE; break
    }
    ll_old <- es$ll
  }
  if (!converged) warning("EM did not converge in ", maxit, " cycles.", call. = FALSE)
  es <- .estep(ctx, par, Sigma)

  vc <- NULL
  if (se != "none") vc <- .vcov_irt(ctx, par, Sigma, es, se)

  dimnames(Sigma) <- list(colnames(Q), colnames(Q))
  fit <- list(call = cl, data = X, labels = pd$labels, model = model, Q = Q, spec = spec,
              par = par, vcov = vc, Sigma = Sigma, var_free = var_free,
              theta = theta, prior = es$prior, posterior = es$post,
              loglik = es$ll + .log_prior(spec, par), npar = spec$npar +
                sum(var_free) + D * (D - 1) / 2,
              nobs = N, iterations = iter, converged = converged, se_method = se)
  class(fit) <- "irtfit"
  fit
}

# E-step -----------------------------------------------------------------------
.prior_weights <- function(theta, Sigma) {
  Si <- solve(Sigma)
  lw <- -0.5 * rowSums((theta %*% Si) * theta)
  w <- exp(lw - max(lw)); w / sum(w)
}

.item_logp <- function(ctx, par, j) {
  it <- ctx$spec$items[[j]]
  log(.item_probs(it$fam, .decode_item(it, par, ctx$D), ctx$theta))
}

.estep <- function(ctx, par, Sigma) {
  prior <- .prior_weights(ctx$theta, Sigma)
  L <- matrix(0, ctx$N, nrow(ctx$theta))
  for (j in seq_len(ctx$J)) {
    LP <- rbind(t(.item_logp(ctx, par, j)), 0)
    L <- L + LP[ctx$idx[[j]], , drop = FALSE]
  }
  L <- sweep(L, 2, log(prior), "+")
  mx <- apply(L, 1, max)
  E <- exp(L - mx)
  rs <- rowSums(E)
  post <- E / rs
  r <- lapply(seq_len(ctx$J), function(j) crossprod(post, ctx$onehot[[j]]))
  list(post = post, r = r, ll = sum(mx + log(rs)), prior = prior)
}

# M-step for one block ------------------------------------------------------------
.mstep_block <- function(ctx, par, block, items, r) {
  spec <- ctx$spec
  sub_spec <- spec; sub_spec$items <- spec$items[items]
  fn <- function(p) {
    full <- par; full[block] <- p
    v <- 0
    for (j in items) v <- v + sum(r[[j]] * .item_logp(ctx, full, j))
    v <- v + .log_prior(sub_spec, full)
    if (!is.finite(v)) return(1e10)
    -v
  }
  opt <- stats::optim(par[block], fn, method = "BFGS",
                      control = list(maxit = 100, reltol = 1e-10))
  opt$par
}

# Gradient of the marginal log-likelihood (Fisher identity) -----------------------
.dlogp <- function(ctx, par, j, h = 1e-5) {
  it <- ctx$spec$items[[j]]
  ids <- it$idx
  out <- array(0, c(nrow(ctx$theta), it$K, length(ids)))
  for (p in seq_along(ids)) {
    up <- par; up[ids[p]] <- up[ids[p]] + h
    dn <- par; dn[ids[p]] <- dn[ids[p]] - h
    out[, , p] <- (.item_logp(ctx, up, j) - .item_logp(ctx, dn, j)) / (2 * h)
  }
  out
}

.marg_grad <- function(ctx, par, Sigma) {
  es <- .estep(ctx, par, Sigma)
  g <- numeric(length(par))
  for (j in seq_len(ctx$J)) {
    ids <- ctx$spec$items[[j]]$idx
    dl <- .dlogp(ctx, par, j)
    for (p in seq_along(ids)) g[ids[p]] <- g[ids[p]] + sum(es$r[[j]] * dl[, , p])
  }
  h <- 1e-5
  for (p in seq_along(par)) {
    up <- par; up[p] <- up[p] + h; dn <- par; dn[p] <- dn[p] - h
    g[p] <- g[p] + (.log_prior(ctx$spec, up) - .log_prior(ctx$spec, dn)) / (2 * h)
  }
  g
}

.vcov_irt <- function(ctx, par, Sigma, es, method) {
  np <- length(par)
  if (method == "hessian") {
    H <- matrix(0, np, np)
    h <- 1e-4
    for (p in seq_len(np)) {
      up <- par; up[p] <- up[p] + h; dn <- par; dn[p] <- dn[p] - h
      H[, p] <- (.marg_grad(ctx, up, Sigma) - .marg_grad(ctx, dn, Sigma)) / (2 * h)
    }
    H <- (H + t(H)) / 2
    vc <- tryCatch(solve(-H), error = function(e) NULL)
    if (!is.null(vc) && all(diag(vc) > 0)) return(vc)
    warning("Observed information not positive definite; using cross-product SEs.",
            call. = FALSE)
  }
  # outer product of per-person gradients
  G <- matrix(0, ctx$N, np)
  for (j in seq_len(ctx$J)) {
    ids <- ctx$spec$items[[j]]$idx
    dl <- .dlogp(ctx, par, j)
    for (p in seq_along(ids)) {
      M <- rbind(t(dl[, , p]), 0)
      G[, ids[p]] <- G[, ids[p]] + rowSums(es$post * M[ctx$idx[[j]], , drop = FALSE])
    }
  }
  tryCatch(solve(crossprod(G)), error = function(e) {
    warning("Could not invert information matrix; SEs unavailable.", call. = FALSE)
    NULL })
}
