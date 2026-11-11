# Computerized adaptive testing -----------------------------------------------------

#' Design a computerized adaptive test
#'
#' Collects everything that defines an adaptive test: the item bank, how the
#' first item is chosen, the item selection rule, the trait estimator,
#' exposure control, content balancing and stopping rules. The design is then
#' used either live, one item at a time ([cat_start()], [cat_next()],
#' [cat_answer()]), or in a simulation study ([cat_simulate()]).
#'
#' @param bank An `irt_bank` (from [item_bank()] or [make_bank()]) or a
#'   unidimensional `irtfit`.
#' @param select Item selection rule: `"MFI"` (maximum Fisher information at
#'   the current estimate), `"MPWI"` (maximum posterior-weighted
#'   information), `"KL"` (Kullback-Leibler information over a shrinking
#'   interval, Chang & Ying 1996), `"MEPV"` (minimum expected posterior
#'   variance) or `"random"`.
#' @param estimator Interim and final trait estimator: `"EAP"`, `"MAP"`,
#'   `"ML"` or `"WLE"`. ML and WLE fall back to EAP until the response pattern
#'   is mixed.
#' @param prior_mean,prior_sd Normal prior for EAP/MAP (defaults to the bank's
#'   latent distribution).
#' @param start_theta Starting trait value for selecting the first item.
#' @param min_items,max_items Minimum and maximum test length.
#' @param se_target Stop when the standard error falls to this value.
#' @param cut Optional cut score(s) on the trait scale for classification
#'   testing: the test stops once the confidence interval of the estimate
#'   excludes every cut (subject to `min_items`).
#' @param cut_level Confidence level used with `cut`.
#' @param exposure Exposure control: `"none"`, `"randomesque"` (choose at
#'   random among the `randomesque_k` best items, Kingsbury & Zara 1989) or
#'   `"sympson-hetter"` (probabilistic filter with exposure parameters
#'   `sh_k`, Sympson & Hetter 1985; see [cat_sympson_hetter()]).
#' @param randomesque_k Number of candidate items for randomesque selection.
#' @param sh_k Sympson-Hetter exposure control parameters, one per item.
#' @param strata Number of a-stratification strata (Chang & Ying 1999):
#'   items are grouped by discrimination and the test moves from low- to
#'   high-discrimination strata. `NULL` (default) for none.
#' @param content_target Named vector of target proportions per content area
#'   (requires content areas in the bank). Items are drawn from the area
#'   furthest below its target (Kingsbury & Zara 1989).
#' @param exclude Items never to administer (names or indices).
#' @param grid Number of quadrature points for interim computations.
#' @return An object of class `cat_design`.
#' @references
#' Chang, H.-H., & Ying, Z. (1996). A global information approach to
#' computerized adaptive testing. *Applied Psychological Measurement, 20*(3),
#' 213-229.
#'
#' Chang, H.-H., & Ying, Z. (1999). a-stratified multistage computerized
#' adaptive testing. *Applied Psychological Measurement, 23*(3), 211-222.
#'
#' Kingsbury, G. G., & Zara, A. R. (1989). Procedures for selecting items for
#' computerized adaptive tests. *Applied Measurement in Education, 2*(4),
#' 359-375.
#'
#' Sympson, J. B., & Hetter, R. D. (1985). Controlling item-exposure rates in
#' computerized adaptive testing. *Proceedings of the 27th Annual Meeting of
#' the Military Testing Association*, 973-977.
#'
#' Magis, D., & Barrada, J. R. (2017). Computerized adaptive testing with R:
#' Recent updates of the package catR. *Journal of Statistical Software,
#' 76*(Code Snippet 1), 1-19.
#' @examples
#' bank <- make_bank(b = rnorm(200), a = rlnorm(200, 0, 0.3), model = "2pl")
#' design <- cat_design(bank, max_items = 20, se_target = 0.35)
#' design
#' @export
cat_design <- function(bank, select = c("MFI", "MPWI", "KL", "MEPV", "random"),
                       estimator = c("EAP", "MAP", "ML", "WLE"),
                       prior_mean = NULL, prior_sd = NULL, start_theta = 0,
                       min_items = 1, max_items = 20, se_target = NULL,
                       cut = NULL, cut_level = 0.95,
                       exposure = c("none", "randomesque", "sympson-hetter"),
                       randomesque_k = 5, sh_k = NULL, strata = NULL,
                       content_target = NULL, exclude = NULL, grid = 121) {
  if (inherits(bank, "irtfit")) bank <- item_bank(bank)
  if (!inherits(bank, "irt_bank")) stop("'bank' must be an irt_bank or irtfit.", call. = FALSE)
  select <- match.arg(toupper(select), c("MFI", "MPWI", "KL", "MEPV", "RANDOM"))
  estimator <- match.arg(toupper(estimator), c("EAP", "MAP", "ML", "WLE"))
  exposure <- match.arg(exposure)
  J <- length(bank$items)
  max_items <- min(max_items, J)
  if (min_items > max_items) stop("'min_items' exceeds 'max_items'.", call. = FALSE)
  if (is.null(se_target) && is.null(cut) && max_items == J)
    message("No stopping rule other than administering the whole bank.")
  if (exposure == "sympson-hetter") {
    if (is.null(sh_k)) stop("Sympson-Hetter needs 'sh_k' (see cat_sympson_hetter()).", call. = FALSE)
    sh_k <- rep(sh_k, length.out = J)
  }
  content <- vapply(bank$items, function(it) if (is.null(it$content)) NA_character_ else it$content, character(1))
  if (!is.null(content_target)) {
    if (all(is.na(content))) stop("'content_target' needs content areas in the bank.", call. = FALSE)
    if (is.null(names(content_target)) || !all(stats::na.omit(content) %in% names(content_target)))
      stop("'content_target' must be named with every content area in the bank.", call. = FALSE)
    content_target <- content_target / sum(content_target)
  }
  excl <- rep(FALSE, J)
  if (!is.null(exclude)) excl[if (is.character(exclude)) match(exclude, names(bank$items)) else exclude] <- TRUE
  pm <- if (is.null(prior_mean)) bank$mu else prior_mean
  ps <- if (is.null(prior_sd)) bank$sigma else prior_sd
  th <- seq(pm - 6 * ps, pm + 6 * ps, length.out = grid)
  # precompute probabilities and information on the grid
  P <- lapply(bank$items, function(it) .bp(it, th))
  info <- vapply(bank$items, function(it) .bi(it, th), numeric(grid))
  a <- vapply(bank$items, `[[`, numeric(1), "a")
  stratum <- if (is.null(strata)) rep(1L, J) else {
    if (strata < 2) stop("'strata' must be at least 2.", call. = FALSE)
    as.integer(cut(rank(a, ties.method = "first"), strata, labels = FALSE))
  }
  structure(list(bank = bank, select = select, estimator = estimator, prior_mean = pm, prior_sd = ps,
                 start_theta = start_theta, min_items = min_items, max_items = max_items,
                 se_target = se_target, cut = cut, cut_level = cut_level, exposure = exposure,
                 randomesque_k = randomesque_k, sh_k = sh_k, strata = strata, stratum = stratum,
                 content = content, content_target = content_target, exclude = excl,
                 theta = th, P = P, info = info,
                 prior = { w <- stats::dnorm(th, pm, ps); w / sum(w) }),
            class = "cat_design")
}

#' @export
print.cat_design <- function(x, ...) {
  cat("<cat_design> ", length(x$bank$items), "-item bank\n", sep = "")
  cat("  selection: ", x$select, if (!is.null(x$strata)) sprintf(" within %d a-strata", x$strata),
      " | estimator: ", x$estimator, sprintf(" (prior N(%.2f, %.2f^2))", x$prior_mean, x$prior_sd), "\n", sep = "")
  cat("  exposure control: ", x$exposure,
      if (x$exposure == "randomesque") sprintf(" (k = %d)", x$randomesque_k), "\n", sep = "")
  if (!is.null(x$content_target))
    cat("  content balancing:", paste(sprintf("%s %.0f%%", names(x$content_target), 100 * x$content_target),
                                     collapse = ", "), "\n")
  stops <- c(sprintf("%d-%d items", x$min_items, x$max_items),
             if (!is.null(x$se_target)) sprintf("SE <= %.3f", x$se_target),
             if (!is.null(x$cut)) sprintf("classification at cut(s) %s (%.0f%% CI)",
                                          paste(round(x$cut, 2), collapse = ", "), 100 * x$cut_level))
  cat("  stopping:", paste(stops, collapse = "; "), "\n")
  invisible(x)
}

# --- state handling ----------------------------------------------------------------

#' Run an adaptive test live, one item at a time
#'
#' `cat_start()` creates a test session from a design, `cat_next()` returns
#' the next item to administer, `cat_answer()` records a response and updates
#' the trait estimate, and `cat_done()` says whether a stopping rule has been
#' met. Sessions are ordinary R lists, so they can be stored between page
#' loads of a web application (for example in 'shiny').
#'
#' @param design A [cat_design()].
#' @param session A session from `cat_start()` (or a previous `cat_answer()`).
#' @param item Item name or index being answered.
#' @param response The response, on the item's original score scale.
#' @return `cat_start()` and `cat_answer()` return the updated session;
#'   `cat_next()` returns the name of the next item (or `NA` when the test
#'   is finished); `cat_done()` returns `TRUE`/`FALSE` with the reason as an
#'   attribute.
#' @examples
#' bank <- make_bank(b = seq(-2, 2, length.out = 40), a = 1.2, model = "2pl")
#' design <- cat_design(bank, max_items = 10)
#' s <- cat_start(design)
#' while (!cat_done(s)) {
#'   it <- cat_next(s)
#'   s <- cat_answer(s, it, response = rbinom(1, 1, 0.6))   # replace with a real answer
#' }
#' s
#' @export
cat_start <- function(design) {
  structure(list(design = design, items = character(0), idx = integer(0), responses = numeric(0),
                 logpost = log(design$prior), loglik = numeric(length(design$theta)),
                 theta = design$start_theta, se = design$prior_sd,
                 history = data.frame(step = integer(0), item = character(0), response = numeric(0),
                                      theta = numeric(0), se = numeric(0))),
            class = "cat_session")
}

#' @rdname cat_start
#' @export
cat_next <- function(session) {
  if (cat_done(session)) return(NA_character_)
  d <- session$design
  j <- .cat_select(d, session)
  names(d$bank$items)[j]
}

#' @rdname cat_start
#' @export
cat_answer <- function(session, item, response) {
  d <- session$design
  j <- if (is.character(item)) match(item, names(d$bank$items)) else as.integer(item)
  if (is.na(j)) stop("Unknown item '", item, "'.", call. = FALSE)
  if (j %in% session$idx) stop("Item '", names(d$bank$items)[j], "' was already administered.", call. = FALSE)
  it <- d$bank$items[[j]]
  k <- match(response, it$scores)
  if (is.na(k)) stop("Response ", response, " is not a valid score for item '", it$name, "'.", call. = FALSE)
  lp <- log(d$P[[j]][, k])
  session$loglik <- session$loglik + lp
  session$logpost <- session$logpost + lp
  session$idx <- c(session$idx, j)
  session$items <- c(session$items, it$name)
  session$responses <- c(session$responses, response)
  est <- .cat_estimate(d, session)
  session$theta <- est[1]; session$se <- est[2]
  session$history <- rbind(session$history, data.frame(step = length(session$idx), item = it$name,
                                                       response = response, theta = est[1], se = est[2]))
  session
}

#' @rdname cat_start
#' @export
cat_done <- function(session) {
  d <- session$design; n <- length(session$idx)
  reason <- NULL
  if (n >= d$max_items) reason <- "max_items"
  else if (!any(.cat_available(d, session))) reason <- "bank_exhausted"
  else if (n >= d$min_items && n > 0) {
    if (!is.null(d$se_target) && session$se <= d$se_target) reason <- "se_target"
    if (is.null(reason) && !is.null(d$cut)) {
      z <- stats::qnorm(1 - (1 - d$cut_level) / 2)
      if (all(abs(session$theta - d$cut) > z * session$se)) reason <- "classified"
    }
  }
  structure(!is.null(reason), reason = reason)
}

#' @export
print.cat_session <- function(x, ...) {
  n <- length(x$idx)
  dn <- cat_done(x)
  cat(sprintf("<cat_session> %d item(s) administered | theta = %.3f, SE = %.3f | %s\n", n, x$theta, x$se,
              if (dn) paste("finished:", attr(dn, "reason")) else "in progress"))
  if (n) print(utils::tail(x$history, 10), row.names = FALSE, digits = 3)
  invisible(x)
}

.cat_available <- function(d, s) {
  av <- !d$exclude
  av[s$idx] <- FALSE
  av
}

.cat_grid_idx <- function(d, theta) {
  i <- findInterval(theta, d$theta, all.inside = TRUE)
  i + (abs(d$theta[i + 1] - theta) < abs(d$theta[i] - theta))
}

.cat_estimate <- function(d, s) {
  th <- d$theta
  post <- exp(s$logpost - max(s$logpost)); post <- post / sum(post)
  eap <- sum(th * post); psd <- sqrt(sum((th - eap)^2 * post))
  est <- d$estimator
  pos <- vapply(seq_along(s$idx), function(i) {
    sc <- d$bank$items[[s$idx[i]]]$scores
    if (s$responses[i] == min(sc)) -1L else if (s$responses[i] == max(sc)) 1L else 0L }, integer(1))
  mixed <- !(all(pos == -1L) || all(pos == 1L))
  if (est == "EAP" || (est %in% c("ML", "WLE") && !mixed)) return(c(eap, psd))
  obj <- function(t) {
    v <- 0
    for (i in seq_along(s$idx)) {
      it <- d$bank$items[[s$idx[i]]]
      v <- v + log(.bp(it, t)[1, match(s$responses[i], it$scores)])
    }
    if (est == "MAP") v <- v + stats::dnorm(t, d$prior_mean, d$prior_sd, log = TRUE)
    if (est == "WLE") v <- v + 0.5 * log(sum(vapply(s$idx, function(j) .bi(d$bank$items[[j]], t), numeric(1))))
    v
  }
  start <- th[which.max(if (est == "MAP") s$logpost else s$loglik)]
  o <- stats::optimize(obj, c(max(min(th), start - 1.5), min(max(th), start + 1.5)), maximum = TRUE)
  t <- o$maximum
  I <- sum(vapply(s$idx, function(j) .bi(d$bank$items[[j]], t), numeric(1)))
  if (est == "MAP") I <- I + 1 / d$prior_sd^2
  c(t, 1 / sqrt(I))
}

.cat_rank <- function(d, s) {
  av <- .cat_available(d, s)
  n <- length(s$idx)
  # a-stratification: current stratum by position in the test
  if (!is.null(d$strata)) {
    cur <- min(d$strata, 1 + floor(n * d$strata / d$max_items))
    in_st <- av & d$stratum == cur
    if (any(in_st)) av <- in_st
  }
  # content balancing: restrict to the most under-represented area
  if (!is.null(d$content_target)) {
    adm <- table(factor(d$content[s$idx], levels = names(d$content_target)))
    prop <- if (n) adm / n else adm * 0
    gap <- d$content_target - prop
    for (area in names(sort(gap, decreasing = TRUE))) {
      cand <- av & d$content %in% area
      if (any(cand)) { av <- cand; break }
    }
  }
  cand <- which(av)
  if (length(cand) == 1) return(cand)
  if (d$select == "RANDOM") return(cand[sample.int(length(cand))])
  theta <- if (n == 0) d$start_theta else s$theta
  post <- exp(s$logpost - max(s$logpost)); post <- post / sum(post)
  crit <- switch(d$select,
    MFI = d$info[.cat_grid_idx(d, theta), cand],
    MPWI = as.vector(post %*% d$info[, cand, drop = FALSE]),
    KL = {
      delta <- 3 / sqrt(max(n, 1))
      gi <- .cat_grid_idx(d, theta)
      win <- abs(d$theta - theta) <= delta
      vapply(cand, function(j) {
        P <- d$P[[j]]; p0 <- matrix(P[gi, ], sum(win), ncol(P), byrow = TRUE)
        sum(p0 * log(p0 / P[win, , drop = FALSE]))
      }, numeric(1))
    },
    MEPV = -vapply(cand, function(j) {
      P <- d$P[[j]]
      pk <- as.vector(post %*% P)
      sum(vapply(seq_len(ncol(P)), function(k) {
        pp <- post * P[, k]; pp <- pp / sum(pp)
        m <- sum(pp * d$theta); pk[k] * sum(pp * (d$theta - m)^2)
      }, numeric(1)))
    }, numeric(1)))
  cand[order(crit, decreasing = TRUE)]
}

.cat_select <- function(d, s) {
  ord <- .cat_rank(d, s)
  if (length(ord) == 1 || d$select == "RANDOM") return(ord[1])
  if (d$exposure == "randomesque") {
    top <- ord[seq_len(min(d$randomesque_k, length(ord)))]
    return(top[sample.int(length(top), 1)])
  }
  if (d$exposure == "sympson-hetter") {
    for (j in ord) if (stats::runif(1) <= d$sh_k[j]) return(j)
    return(ord[length(ord)])
  }
  ord[1]
}

# --- simulation --------------------------------------------------------------------

#' Simulate adaptive tests
#'
#' Runs a Monte Carlo (or post-hoc) simulation of a CAT design. Simulees either
#' have known true trait values, from which responses are generated with the
#' bank's item parameters, or bring a complete matrix of real responses
#' (post-hoc simulation) from which the adaptive test reads the answers.
#'
#' @param design A [cat_design()].
#' @param theta True trait values (Monte Carlo simulation).
#' @param responses Optional matrix of responses with one column per bank
#'   item (post-hoc simulation). Overrides `theta`.
#' @param seed Optional random seed.
#' @param progress Print progress.
#' @return An object of class `cat_sim` with per-person results (`results`),
#'   the full administration history (`history`), item exposure rates
#'   (`exposure`) and the design.
#' @examples
#' bank <- make_bank(b = rnorm(150), a = rlnorm(150, 0, 0.3), model = "2pl")
#' design <- cat_design(bank, max_items = 25, se_target = 0.3, exposure = "randomesque")
#' sim <- cat_simulate(design, theta = rnorm(200), seed = 1)
#' summary(sim)
#' @export
cat_simulate <- function(design, theta = NULL, responses = NULL, seed = NULL, progress = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  items <- design$bank$items
  if (!is.null(responses)) {
    responses <- as.matrix(responses)
    if (ncol(responses) != length(items)) stop("'responses' needs one column per bank item.", call. = FALSE)
    N <- nrow(responses); theta <- rep(NA_real_, N)
  } else {
    if (is.null(theta)) stop("Supply 'theta' or 'responses'.", call. = FALSE)
    N <- length(theta)
  }
  res <- vector("list", N); hist <- vector("list", N)
  for (i in seq_len(N)) {
    di <- design
    if (!is.null(responses)) di$exclude <- di$exclude | is.na(responses[i, ])
    s <- cat_start(di)
    while (!cat_done(s)) {
      j <- .cat_select(di, s)
      if (!is.null(responses)) {
        r <- responses[i, j]
      } else {
        p <- .bp(items[[j]], theta[i])
        r <- items[[j]]$scores[sum(stats::runif(1) > cumsum(p)) + 1]
      }
      s <- cat_answer(s, j, r)
    }
    dn <- cat_done(s)
    res[[i]] <- data.frame(person = i, true_theta = theta[i], theta = s$theta, se = s$se,
                           n_items = length(s$idx), stop = attr(dn, "reason"))
    hist[[i]] <- cbind(person = i, s$history)
    if (progress && i %% 100 == 0) message(i, " / ", N)
  }
  results <- do.call(rbind, res); history <- do.call(rbind, hist)
  if (!is.null(design$cut)) {
    results$class <- findInterval(results$theta, sort(design$cut))
    if (all(is.finite(results$true_theta))) results$true_class <- findInterval(results$true_theta, sort(design$cut))
  }
  cnt <- table(factor(history$item, levels = names(items)))
  exposure <- data.frame(item = names(items), exposure = as.vector(cnt) / N,
                         a = vapply(items, `[[`, numeric(1), "a"),
                         b = vapply(items, function(it) mean(-it$d / it$a), numeric(1)))
  structure(list(results = results, history = history, exposure = exposure, design = design),
            class = "cat_sim")
}

#' Calibrate Sympson-Hetter exposure parameters
#'
#' Iteratively simulates the CAT and adjusts each item's exposure control
#' parameter so that no item is administered to more than `r_max` of
#' examinees (Sympson & Hetter, 1985).
#'
#' @param design A [cat_design()] (its exposure setting is ignored).
#' @param theta Simulated trait values representative of the population.
#' @param r_max Target maximum exposure rate.
#' @param iterations Number of adjustment cycles.
#' @return A numeric vector of exposure parameters, for `cat_design(..., sh_k = )`.
#' @export
cat_sympson_hetter <- function(design, theta, r_max = 0.25, iterations = 8) {
  J <- length(design$bank$items)
  k <- rep(1, J)
  N <- length(theta)
  for (cycle in seq_len(iterations)) {
    sel <- numeric(J)
    for (i in seq_len(N)) {
      s <- cat_start(design)
      while (!cat_done(s)) {
        ord <- .cat_rank(design, s)
        j <- ord[length(ord)]
        for (cand in ord) {             # every item passed to the filter counts as selected
          sel[cand] <- sel[cand] + 1
          if (stats::runif(1) <= k[cand]) { j <- cand; break }
        }
        p <- .bp(design$bank$items[[j]], theta[i])
        r <- design$bank$items[[j]]$scores[sum(stats::runif(1) > cumsum(p)) + 1]
        s <- cat_answer(s, j, r)
      }
    }
    psel <- sel / N
    k <- ifelse(psel > r_max, r_max / psel, 1)
  }
  k
}

#' @export
summary.cat_sim <- function(object, ...) {
  r <- object$results; e <- object$exposure$exposure
  out <- list(
    n = nrow(r),
    length = c(mean = mean(r$n_items), sd = stats::sd(r$n_items), min = min(r$n_items), max = max(r$n_items)),
    stop = table(r$stop),
    exposure = c(max = max(e), unused_pct = 100 * mean(e == 0),
                 # Chen, Ankenmann & Spray (2003) test overlap rate
                 overlap = { n <- nrow(r); L <- mean(r$n_items); n / (n - 1) * sum(e^2) / L - 1 / (n - 1) },
                 chisq = sum((e - mean(r$n_items) / length(e))^2 / (mean(r$n_items) / length(e)))),
    mean_se = mean(r$se))
  if (all(is.finite(r$true_theta))) {
    out$recovery <- c(bias = mean(r$theta - r$true_theta), rmse = sqrt(mean((r$theta - r$true_theta)^2)),
                      correlation = stats::cor(r$theta, r$true_theta))
    if (!is.null(r$class) && !is.null(r$true_class)) out$classification_accuracy <- mean(r$class == r$true_class)
  }
  class(out) <- "summary.cat_sim"
  out
}

#' @export
print.summary.cat_sim <- function(x, ...) {
  cat("CAT simulation:", x$n, "simulees\n")
  cat(sprintf("  test length: mean %.1f (SD %.1f), range %d-%d\n", x$length[1], x$length[2], x$length[3], x$length[4]))
  cat("  stopping reasons:", paste(sprintf("%s %d", names(x$stop), x$stop), collapse = ", "), "\n")
  cat(sprintf("  mean final SE: %.3f\n", x$mean_se))
  if (!is.null(x$recovery))
    cat(sprintf("  recovery: bias %.3f, RMSE %.3f, r = %.3f\n", x$recovery[1], x$recovery[2], x$recovery[3]))
  if (!is.null(x$classification_accuracy))
    cat(sprintf("  classification accuracy: %.1f%%\n", 100 * x$classification_accuracy))
  cat(sprintf("  exposure: max %.3f, %.1f%% of bank unused, test overlap %.3f, chi-square %.1f\n",
              x$exposure[1], x$exposure[2], x$exposure[3], x$exposure[4]))
  invisible(x)
}

#' @export
print.cat_sim <- function(x, ...) print(summary(x))

#' @export
as.data.frame.cat_sim <- function(x, ...) x$results

#' Conditional CAT performance
#'
#' Bias, RMSE, mean SE and mean test length within bins of the true trait.
#' @param x A `cat_sim` with known true trait values.
#' @param bins Number of equal-count bins.
#' @export
cat_conditional <- function(x, bins = 10) {
  r <- x$results
  if (!all(is.finite(r$true_theta))) stop("Conditional results need true theta values.", call. = FALSE)
  br <- unique(stats::quantile(r$true_theta, seq(0, 1, length.out = bins + 1)))
  b <- cut(r$true_theta, br, include.lowest = TRUE)
  agg <- function(f) as.vector(tapply(seq_len(nrow(r)), b, function(i) f(r[i, ])))
  data.frame(bin = levels(b), theta = agg(function(z) mean(z$true_theta)), n = agg(nrow),
             bias = agg(function(z) mean(z$theta - z$true_theta)),
             rmse = agg(function(z) sqrt(mean((z$theta - z$true_theta)^2))),
             mean_se = agg(function(z) mean(z$se)), mean_length = agg(function(z) mean(z$n_items)))
}
