# Item response functions and parameter handling -----------------------------
#
# Internally every item is written in slope-intercept form:
#   z = theta %*% a            (a is a length-D slope vector, 0 where not loaded)
#   dich : P(1) = g + (1 - g) * plogis(z + d)
#   grm  : P(X >= k) = plogis(z + d_k),  d_1 > d_2 > ... > d_{K-1}
#   gpcm : P(X = k) propto exp(k * z + sum_{v <= k} d_v)
# Reported IRT parameters are derived from these (e.g. b = -d / a).

.dich_models <- c("rasch", "2pl", "3pl")
.poly_models <- c("grm", "pcm", "gpcm", "rsm")
.all_models  <- c(.dich_models, .poly_models)

.family_of <- function(model) {
  switch(model,
         rasch = , `2pl` = , `3pl` = "dich",
         grm = "grm",
         pcm = , gpcm = , rsm = "gpcm")
}

.slope_free <- function(model) model %in% c("2pl", "3pl", "grm", "gpcm")

# Category probabilities: returns Q x K matrix
.item_probs <- function(fam, pars, theta) {
  z <- as.vector(theta %*% pars$a)
  K <- pars$K
  if (fam == "dich") {
    p1 <- pars$g + (1 - pars$g) * stats::plogis(z + pars$d)
    P <- cbind(1 - p1, p1)
  } else if (fam == "grm") {
    cum <- stats::plogis(outer(z, pars$d, "+"))
    P <- cbind(1, cum) - cbind(cum, 0)
  } else {
    k <- 0:(K - 1)
    num <- outer(z, k) + matrix(c(0, cumsum(pars$d)), length(z), K, byrow = TRUE)
    num <- num - apply(num, 1, max)
    e <- exp(num)
    P <- e / rowSums(e)
  }
  P[P < 1e-300] <- 1e-300
  P
}

# Build the parameter specification ------------------------------------------

.build_spec <- function(model, K, Q, priors) {
  J <- length(model)
  D <- ncol(Q)
  items <- vector("list", J)
  pos <- 0L
  rsm_tau_idx <- NULL
  if (any(model == "rsm")) {
    Kr <- unique(K[model == "rsm"])
    if (length(Kr) > 1)
      stop("All 'rsm' items must have the same number of categories.", call. = FALSE)
    if (Kr > 2) {
      rsm_tau_idx <- pos + seq_len(Kr - 2)
      pos <- pos + (Kr - 2)
    }
  }
  for (j in seq_len(J)) {
    m <- model[j]
    load <- Q[j, ] != 0
    n_a <- if (.slope_free(m)) sum(load) else 0L
    n_d <- if (m == "rsm") 1L else K[j] - 1L
    n_g <- if (m == "3pl") 1L else 0L
    idx_a <- if (n_a) pos + seq_len(n_a) else integer(0); pos <- pos + n_a
    idx_d <- pos + seq_len(n_d); pos <- pos + n_d
    idx_g <- if (n_g) pos + 1L else integer(0); pos <- pos + n_g
    own <- c(idx_a, idx_d, idx_g)
    items[[j]] <- list(model = m, fam = .family_of(m), K = K[j], load = load,
                       idx_a = idx_a, idx_d = idx_d, idx_g = idx_g,
                       idx_tau = if (m == "rsm") rsm_tau_idx else integer(0),
                       idx = sort(c(own, if (m == "rsm") rsm_tau_idx)))
  }
  list(items = items, npar = pos, J = J, D = D, Q = Q, rsm_tau_idx = rsm_tau_idx,
       priors = priors)
}

# Decode the free parameter vector for item j into list(a, d, g, K)
.decode_item <- function(it, par, D) {
  a <- numeric(D)
  if (length(it$idx_a)) a[it$load] <- par[it$idx_a] else a[it$load] <- 1
  g <- 0
  if (it$fam == "dich") {
    d <- par[it$idx_d]
    if (length(it$idx_g)) g <- stats::plogis(par[it$idx_g])
  } else if (it$fam == "grm") {
    raw <- par[it$idx_d]
    d <- if (length(raw) > 1) raw[1] - c(0, cumsum(exp(raw[-1]))) else raw
  } else if (it$model == "rsm") {
    tau_free <- par[it$idx_tau]
    tau <- c(tau_free, -sum(tau_free))
    if (it$K == 2) tau <- 0
    d <- -(par[it$idx_d] + tau)
  } else {
    d <- par[it$idx_d]
  }
  list(a = a, d = d, g = g, K = it$K)
}

.log_prior <- function(spec, par) {
  pr <- spec$priors
  lp <- 0
  for (it in spec$items) {
    if (length(it$idx_g) && !is.null(pr$g)) {
      gg <- stats::plogis(par[it$idx_g])
      lp <- lp + stats::dbeta(gg, pr$g[1], pr$g[2], log = TRUE)
    }
    if (length(it$idx_a) && !is.null(pr$a)) {
      aa <- par[it$idx_a]
      lp <- lp + sum(stats::dlnorm(pmax(aa, 1e-8), pr$a[1], pr$a[2], log = TRUE))
    }
  }
  lp
}

# Starting values ---------------------------------------------------------------
.start_values <- function(spec, X) {
  par <- numeric(spec$npar)
  for (j in seq_len(spec$J)) {
    it <- spec$items[[j]]
    x <- X[, j]; x <- x[!is.na(x)]
    tab <- tabulate(x + 1L, nbins = it$K) + 0.5
    if (length(it$idx_a)) par[it$idx_a] <- 1
    if (it$fam == "dich") {
      p <- tab[2] / sum(tab)
      par[it$idx_d] <- stats::qlogis(p) * 1.2
      if (length(it$idx_g)) { par[it$idx_g] <- stats::qlogis(0.2)
        par[it$idx_d] <- stats::qlogis(min(max((p - 0.2) / 0.8, 0.05), 0.95)) }
    } else if (it$fam == "grm") {
      cp <- rev(cumsum(rev(tab)))[-1] / sum(tab)
      dd <- sort(stats::qlogis(pmin(pmax(cp, 0.01), 0.99)), decreasing = TRUE) * 1.2
      dd <- dd - c(0, seq_len(length(dd) - 1) * 1e-3)
      par[it$idx_d] <- if (length(dd) > 1) c(dd[1], log(pmax(-diff(dd), 0.05))) else dd
    } else if (it$model == "rsm") {
      par[it$idx_d] <- -stats::qlogis(min(max(mean(x) / (it$K - 1), 0.02), 0.98))
    } else {
      par[it$idx_d] <- log(tab[-1] / tab[-it$K])
    }
  }
  par
}

# Prepare response data ---------------------------------------------------------
.prep_data <- function(data) {
  if (is.data.frame(data)) data <- as.matrix(data.frame(lapply(data, function(v) {
    if (is.factor(v)) as.integer(v) - 1L else as.numeric(v)
  }), check.names = FALSE))
  if (!is.matrix(data)) stop("'data' must be a matrix or data frame.", call. = FALSE)
  if (is.null(colnames(data))) colnames(data) <- paste0("Item", seq_len(ncol(data)))
  storage.mode(data) <- "double"
  J <- ncol(data)
  X <- matrix(NA_integer_, nrow(data), J, dimnames = list(rownames(data), colnames(data)))
  labels <- vector("list", J); names(labels) <- colnames(data)
  for (j in seq_len(J)) {
    v <- data[, j]
    u <- sort(unique(v[!is.na(v)]))
    if (length(u) < 2)
      stop("Item '", colnames(data)[j], "' has fewer than two observed categories.",
           call. = FALSE)
    X[, j] <- match(v, u) - 1L
    labels[[j]] <- u
  }
  allna <- rowSums(!is.na(X)) == 0
  if (any(allna)) stop("Some persons have no responses; remove them first.", call. = FALSE)
  list(X = X, labels = labels, K = vapply(labels, length, integer(1)))
}
