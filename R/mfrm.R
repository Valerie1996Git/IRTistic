#' Fit a many-facet Rasch model
#'
#' Estimates a many-facet Rasch measurement (MFRM) model of the kind fitted
#' by Linacre's FACETS program, for example persons rated by raters on
#' several criteria. The log-odds of scoring in category k rather than k-1 is
#' \deqn{\theta_n - \sum_f \beta_{f,e} - \tau_{s,k}}
#' where \eqn{\theta_n} is the person measure, \eqn{\beta_{f,e}} the measure
#' of element e of facet f (item difficulty, rater severity, ...) and
#' \eqn{\tau_{s,k}} the Andrich thresholds, shared by all observations
#' (rating scale model) or specific to each element of a `step_facet`
#' (partial credit model).
#'
#' Estimation is by marginal maximum likelihood with persons drawn from
#' \eqn{N(\mu, \sigma^2)}. As in FACETS, the non-person facets are centred
#' at zero and the person mean is estimated. Higher facet measures mean
#' "harder to score highly on" (more difficult items, more severe raters).
#'
#' @param data A long-format data frame with one row per rating.
#' @param person Name of the person column.
#' @param facets Character vector naming the facet columns (e.g.
#'   `c("rater", "criterion")`).
#' @param score Name of the score column (ordered numeric categories).
#' @param model `"rsm"` (rating scale, default) or `"pcm"` (partial credit).
#' @param step_facet For `"pcm"`, the facet whose elements get their own
#'   thresholds (defaults to the first facet).
#' @param quadpts Number of quadrature points.
#' @param maxit,tol Control for the EM warm-up (at most 10 cycles) that is
#'   followed by direct quasi-Newton maximisation of the marginal likelihood.
#' @param se Compute standard errors (`TRUE`).
#' @return An object of class `mfrmfit`.
#' @examples
#' set.seed(1)
#' ratings <- sim_mfrm(n_person = 40, n_rater = 4, n_criterion = 3, ncat = 4)
#' fit <- mfrm_fit(ratings, "person", c("rater", "criterion"), "score")
#' fit
#' \donttest{
#' facet_report(fit)
#' plot(fit, "facets")
#' }
#' @export
mfrm_fit <- function(data, person, facets, score, model = c("rsm", "pcm"),
                     step_facet = NULL, quadpts = 61, maxit = 500, tol = 1e-5,
                     se = TRUE) {
  cl <- match.call()
  model <- match.arg(tolower(model), c("rsm", "pcm"))
  data <- as.data.frame(data)
  cols <- c(person, facets, score)
  miss <- setdiff(cols, names(data))
  if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  keep <- stats::complete.cases(data[, cols])
  if (any(!keep)) message("Dropping ", sum(!keep), " row(s) with missing values.")
  data <- data[keep, cols]
  if (model == "pcm") {
    if (is.null(step_facet)) step_facet <- facets[1]
    if (!step_facet %in% facets) stop("'step_facet' must be one of 'facets'.", call. = FALSE)
  } else step_facet <- NULL

  pf <- factor(data[[person]])
  ef <- lapply(facets, function(f) factor(data[[f]])); names(ef) <- facets
  u <- sort(unique(as.numeric(data[[score]])))
  x <- match(as.numeric(data[[score]]), u) - 1L
  K <- length(u)
  if (K < 2) stop("Scores must have at least two categories.", call. = FALSE)
  sg <- if (is.null(step_facet)) rep(1L, length(x)) else as.integer(ef[[step_facet]])
  S <- max(sg)
  nf <- vapply(ef, nlevels, integer(1))
  if (any(nf < 2)) stop("Every facet needs at least two elements.", call. = FALSE)
  if (model == "pcm") {
    tabs <- table(sg, x)
    if (any(tabs == 0)) warning("Some categories are unobserved for some '", step_facet,
                                "' elements; thresholds may be unstable.", call. = FALSE)
  }

  nodes <- matrix(seq(-6, 6, length.out = quadpts))
  ctx <- list(x = x, p = as.integer(pf), N = nlevels(pf), e = lapply(ef, as.integer),
              nf = nf, K = K, S = S, sg = sg, nodes = nodes[, 1], facets = facets)
  ctx$npar_f <- sum(nf - 1); ctx$npar_t <- S * max(K - 2, 0)
  par <- numeric(ctx$npar_f + ctx$npar_t)
  # sensible starting thresholds from category frequencies
  if (K > 2) {
    tab <- tabulate(x + 1L, K) + 0.5
    t0 <- -log(tab[-1] / tab[-K]); t0 <- t0 - mean(t0)
    par[ctx$npar_f + seq_len(ctx$npar_t)] <- rep(t0[-(K - 1)], S)
  }
  mu <- 0; sigma <- 1
  ll_old <- -Inf
  for (iter in seq_len(min(maxit, 3))) {  # EM warm-up
    es <- .mfrm_estep(ctx, par, mu, sigma)
    W <- es$post[ctx$p, , drop = FALSE]
    cache <- new.env()
    obj <- function(p) {
      if (!identical(cache$p, p)) { cache$p <- p; cache$v <- .mfrm_obj(ctx, p, W, grad = TRUE) }
      cache$v
    }
    # generalised EM: a few quasi-Newton steps per cycle are enough
    opt <- stats::optim(par, function(p) -obj(p)$f, function(p) -obj(p)$g,
                        method = "BFGS", control = list(maxit = 4, reltol = 1e-12))
    par_old <- par; par <- opt$par
    mu <- sum(es$post %*% ctx$nodes) / ctx$N
    sigma <- sqrt(sum(es$post %*% (ctx$nodes - mu)^2) / ctx$N)
    if (abs(es$ll - ll_old) < tol * 1e3 && max(abs(par - par_old)) < 1e-2) break
    ll_old <- es$ll
  }
  # finish by direct quasi-Newton maximisation of the marginal likelihood,
  # with the gradient from Fisher's identity (much faster than slow EM tails)
  np <- length(par)
  mcache <- new.env()
  marg <- function(pf) {
    if (!identical(mcache$p, pf)) {
      mcache$p <- pf
      p <- pf[seq_len(np)]; m <- pf[np + 1]; s <- exp(pf[np + 2])
      e2 <- .mfrm_estep(ctx, p, m, s)
      g <- .mfrm_obj(ctx, p, e2$post[ctx$p, , drop = FALSE], grad = TRUE)$g
      z <- (ctx$nodes - m) / s
      cs <- colSums(e2$post)
      gm <- sum(cs * (z / s - sum(e2$prior * z / s)))
      gs <- sum(cs * (z^2 - sum(e2$prior * z^2)))
      mcache$v <- list(f = e2$ll, g = c(g, gm, gs))
    }
    mcache$v
  }
  opt <- stats::optim(c(par, mu, log(sigma)), function(p) -marg(p)$f, function(p) -marg(p)$g,
                      method = "BFGS", control = list(maxit = 500, reltol = 1e-12))
  par <- opt$par[seq_len(np)]; mu <- opt$par[np + 1]; sigma <- exp(opt$par[np + 2])
  converged <- opt$convergence == 0
  if (!converged) warning("Optimisation did not converge.", call. = FALSE)
  es <- .mfrm_estep(ctx, par, mu, sigma)

  vc <- NULL
  if (se && length(par)) {
    np <- length(par); H <- matrix(0, np, np); h <- 1e-4
    mg <- function(p) {
      e2 <- .mfrm_estep(ctx, p, mu, sigma)
      .mfrm_obj(ctx, p, e2$post[ctx$p, , drop = FALSE], grad = TRUE)$g
    }
    for (k in seq_len(np)) {
      up <- par; up[k] <- up[k] + h; dn <- par; dn[k] <- dn[k] - h
      H[, k] <- (mg(up) - mg(dn)) / (2 * h)
    }
    H <- (H + t(H)) / 2
    vc <- tryCatch(solve(-H), error = function(e) {
      warning("Information matrix is singular; SEs unavailable.", call. = FALSE); NULL })
  }

  fit <- list(call = cl, data = data, person = person, facets = facets, score = score,
              model = model, step_facet = step_facet, ctx = ctx, par = par, vcov = vc,
              levels = c(stats::setNames(list(levels(pf)), person), lapply(ef, levels)),
              categories = u, mu = mu, sigma = sigma, posterior = es$post,
              loglik = es$ll, npar = length(par) + 2, nobs = length(x),
              iterations = iter + opt$counts[["function"]], converged = converged)
  class(fit) <- "mfrmfit"
  fit
}

# Parameter maps ---------------------------------------------------------------
.mfrm_unpack <- function(ctx, par) {
  beta <- list(); pos <- 0
  for (f in seq_along(ctx$nf)) {
    b <- par[pos + seq_len(ctx$nf[f] - 1)]; pos <- pos + ctx$nf[f] - 1
    beta[[f]] <- c(b, -sum(b))
  }
  names(beta) <- ctx$facets
  tau <- matrix(0, ctx$S, ctx$K - 1)
  if (ctx$K > 2) for (s in seq_len(ctx$S)) {
    t0 <- par[pos + seq_len(ctx$K - 2)]; pos <- pos + ctx$K - 2
    tau[s, ] <- c(t0, -sum(t0))
  }
  list(beta = beta, tau = tau)
}

.mfrm_eta <- function(ctx, beta) {
  eta <- numeric(length(ctx$x))
  for (f in seq_along(beta)) eta <- eta - beta[[f]][ctx$e[[f]]]
  eta
}

# Category probabilities: a list of K matrices (R x Q) plus the log-probability
# of each observed rating. Written with plain matrices for speed.
.mfrm_probs <- function(ctx, lin, tau, x = ctx$x, sg = ctx$sg) {
  K <- ctx$K
  ct <- if (K == 2) cbind(0, tau) else cbind(0, t(apply(tau, 1, cumsum)))
  ctr <- ct[sg, , drop = FALSE]
  A <- lapply(seq_len(K), function(k) (k - 1) * lin - ctr[, k])
  shift <- if (max(abs(lin)) * (K - 1) > 500) Reduce(pmax, A) else 0
  E <- lapply(A, function(a) exp(a - shift))
  den <- Reduce(`+`, E)
  logden <- log(den) + shift
  obs <- A[[1]]
  for (k in 2:K) { r <- x == k - 1; if (any(r)) obs[r, ] <- A[[k]][r, , drop = FALSE] }
  list(P = lapply(E, function(e) e / den), logobs = obs - logden)
}

.mfrm_estep <- function(ctx, par, mu, sigma) {
  pp <- .mfrm_unpack(ctx, par)
  lin <- outer(.mfrm_eta(ctx, pp$beta), ctx$nodes, "+")
  LP <- .mfrm_probs(ctx, lin, pp$tau)$logobs
  L <- rowsum(LP, ctx$p, reorder = TRUE)
  prior <- stats::dnorm(ctx$nodes, mu, sigma); prior <- prior / sum(prior)
  L <- sweep(L, 2, log(prior), "+")
  mx <- apply(L, 1, max); E <- exp(L - mx); rs <- rowSums(E)
  list(post = E / rs, ll = sum(mx + log(rs)), prior = prior)
}

.mfrm_obj <- function(ctx, par, W, grad = FALSE) {
  pp <- .mfrm_unpack(ctx, par)
  lin <- outer(.mfrm_eta(ctx, pp$beta), ctx$nodes, "+")
  pr <- .mfrm_probs(ctx, lin, pp$tau)
  f <- sum(W * pr$logobs)
  if (!grad) return(list(f = f))
  P <- pr$P; K <- ctx$K
  Ek <- P[[2]]; if (K > 2) for (k in 3:K) Ek <- Ek + (k - 1) * P[[k]]
  geta <- ctx$x - rowSums(W * Ek)
  g <- numeric(length(par)); pos <- 0
  for (f_i in seq_along(ctx$nf)) {
    gb <- -as.vector(rowsum(geta, ctx$e[[f_i]], reorder = TRUE))
    n <- ctx$nf[f_i]
    g[pos + seq_len(n - 1)] <- gb[-n] - gb[n]; pos <- pos + n - 1
  }
  if (K > 2) {
    gt <- matrix(0, ctx$S, K - 1)
    Gv <- P[[K]]
    for (v in (K - 1):1) {
      # P(X >= v): categories v..K-1
      if (v < K - 1) Gv <- Gv + P[[v + 1]]
      term <- -(ctx$x >= v) + rowSums(W * Gv)
      gt[, v] <- as.vector(rowsum(term, ctx$sg, reorder = TRUE))
    }
    for (s in seq_len(ctx$S)) {
      g[pos + seq_len(K - 2)] <- gt[s, -(K - 1)] - gt[s, K - 1]; pos <- pos + K - 2
    }
  }
  list(f = f, g = g)
}

# Standard errors for all elements (including the constrained last one)
.mfrm_se <- function(fit) {
  ctx <- fit$ctx; np <- length(fit$par)
  if (is.null(fit$vcov)) return(NULL)
  # Jacobian from free parameters to full element vector
  M <- matrix(0, sum(ctx$nf) + ctx$S * (ctx$K - 1), np)
  r <- 0; cpos <- 0
  for (f in seq_along(ctx$nf)) {
    n <- ctx$nf[f]
    M[r + seq_len(n - 1), cpos + seq_len(n - 1)] <- diag(n - 1)
    M[r + n, cpos + seq_len(n - 1)] <- -1
    r <- r + n; cpos <- cpos + n - 1
  }
  if (ctx$K > 2) for (s in seq_len(ctx$S)) {
    n <- ctx$K - 1
    M[r + seq_len(n - 1), cpos + seq_len(n - 1)] <- diag(n - 1)
    M[r + n, cpos + seq_len(n - 1)] <- -1
    r <- r + n; cpos <- cpos + n - 1
  }
  sqrt(pmax(diag(M %*% fit$vcov %*% t(M)), 0))
}

# Person measures by row-vectorised Newton steps ------------------------------
.mfrm_person_measures <- function(fit, method = c("WLE", "ML", "EAP")) {
  method <- match.arg(method)
  ctx <- fit$ctx
  if (method == "EAP") {
    est <- as.vector(fit$posterior %*% ctx$nodes)
    sd <- sqrt(pmax(as.vector(fit$posterior %*% ctx$nodes^2) - est^2, 0))
    return(list(theta = est, se = sd))
  }
  pp <- .mfrm_unpack(ctx, fit$par)
  eta <- .mfrm_eta(ctx, pp$beta)
  stats_at <- function(th) {
    P <- do.call(cbind, .mfrm_probs(ctx, matrix(th[ctx$p] + eta), pp$tau)$P)
    k <- 0:(ctx$K - 1)
    E <- as.vector(P %*% k); V <- as.vector(P %*% k^2) - E^2
    V3 <- as.vector(P %*% k^3) - 3 * E * as.vector(P %*% k^2) + 2 * E^3
    list(g = as.vector(rowsum(ctx$x - E, ctx$p, reorder = TRUE)),
         I = as.vector(rowsum(V, ctx$p, reorder = TRUE)),
         J = as.vector(rowsum(V3, ctx$p, reorder = TRUE)))
  }
  th <- rep(fit$mu, ctx$N)
  for (it in 1:100) {
    s <- stats_at(th)
    g <- s$g + if (method == "WLE") s$J / (2 * s$I) else 0
    step <- pmax(pmin(g / s$I, 1), -1)
    th <- pmin(pmax(th + step, -10), 10)
    if (max(abs(step)) < 1e-7) break
  }
  list(theta = th, se = 1 / sqrt(stats_at(th)$I))
}
