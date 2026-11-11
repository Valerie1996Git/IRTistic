#' IRT-based differential item functioning
#'
#' Calibrates the items separately in each group, places every focal group
#' on the reference group's scale by characteristic-curve linking on anchor
#' items, and then compares the item parameters with Lord's Wald chi-square
#' test and the item response functions with Raju's area measures and
#' signed/unsigned expected-score differences weighted by the focal group's
#' trait distribution (SIDS/UIDS). Works for dichotomous and polytomous
#' items and for any number of groups. Anchors can be purified iteratively.
#'
#' @param data A response matrix/data frame, or an `irtfit` (its data and
#'   item models are reused).
#' @param group A vector of group memberships, one per person.
#' @param model Item model(s) when `data` is not a fitted model.
#' @param reference Reference group (defaults to the first level).
#' @param anchors Names of anchor items assumed DIF-free. When `NULL`, all
#'   items start as anchors and, with `purify = TRUE`, flagged items are
#'   removed and the linking repeated until the anchor set is stable.
#' @param purify Iteratively purify the anchor set.
#' @param alpha Significance level used for flagging and purification.
#' @param link_method Linking method passed to [irt_link()].
#' @param p.adjust Multiple-comparison adjustment.
#' @return A data frame of class `irt_dif_irt`, one row per item and focal
#'   group, with Lord's chi-square (`chisq`, `df`, `p`, `p_adj`), the mean
#'   difficulty difference (`delta_b`, focal minus reference; positive means
#'   harder for the focal group), Raju's signed and unsigned areas, SIDS and
#'   UIDS (expected item score differences, focal minus reference, in score
#'   points), and a `flag`. The differential test functioning (expected total
#'   score difference) is in `attr(, "DTF")`.
#' @examples
#' \donttest{
#' set.seed(3)
#' d <- sim_irt(1000, 10, "2pl")
#' g <- rep(c("ref", "foc"), 500)
#' th <- attr(d, "truth")$theta[, 1]
#' d$I04[g == "foc"] <- rbinom(500, 1, plogis(1.2 * (th[g == "foc"] - 1)))
#' res <- dif_irt(d, g, model = "2pl", reference = "ref")
#' res
#' plot_dif_curves(res)
#' }
#' @export
dif_irt <- function(data, group, model = "2pl", reference = NULL, anchors = NULL,
                    purify = TRUE, alpha = 0.05, link_method = "stocking-lord", p.adjust = "BH") {
  if (inherits(data, "irtfit")) {
    model <- data$model
    X <- data$data
    for (j in seq_len(ncol(X))) X[, j] <- data$labels[[j]][X[, j] + 1]
  } else {
    X <- data.matrix(as.data.frame(data))
    if (is.null(colnames(X))) colnames(X) <- paste0("Item", seq_len(ncol(X)))
  }
  group <- as.factor(group)
  if (length(group) != nrow(X)) stop("'group' must have one value per person.", call. = FALSE)
  if (is.null(reference)) reference <- levels(group)[1]
  if (!reference %in% levels(group)) stop("'reference' is not a level of 'group'.", call. = FALSE)
  focal <- setdiff(levels(group), reference)
  items <- colnames(X)
  # check categories are comparable across groups
  for (j in items) {
    cats <- lapply(levels(group), function(l) sort(unique(stats::na.omit(X[group == l, j]))))
    if (length(unique(vapply(cats, paste, character(1), collapse = ","))) > 1)
      stop("Item '", j, "' does not have the same observed categories in every group; ",
           "collapse sparse categories first.", call. = FALSE)
  }
  fits <- lapply(levels(group), function(l) irt_fit(X[group == l, , drop = FALSE], model = model))
  names(fits) <- levels(group)
  user_anchors <- !is.null(anchors)
  res_all <- list(); links <- list(); dtf <- c()
  for (fg in focal) {
    anc <- if (user_anchors) anchors else items
    for (iter in 1:10) {
      lk <- irt_link(fits[[fg]], fits[[reference]], anchors = anc, method = link_method)
      res <- .dif_irt_stats(fits[[reference]], fits[[fg]], lk)
      res$p_adj <- stats::p.adjust(res$p, p.adjust)
      if (!purify || user_anchors) break
      new_anc <- items[!(res$p_adj < alpha)]
      if (length(new_anc) < 2) {
        warning("Purification left fewer than two anchors for group '", fg,
                "'; using the previous anchor set.", call. = FALSE)
        break
      }
      if (setequal(new_anc, anc)) break
      anc <- new_anc
    }
    res$focal <- fg
    res$anchor <- res$item %in% lk$anchors
    res$flag <- res$p_adj < alpha
    res_all[[fg]] <- res; links[[fg]] <- lk
    dtf[fg] <- sum(res$sids)
  }
  out <- do.call(rbind, res_all); rownames(out) <- NULL
  out <- out[, c("item", "focal", "chisq", "df", "p", "p_adj", "delta_b", "signed_area",
                 "unsigned_area", "sids", "uids", "anchor", "flag")]
  attr(out, "fits") <- fits; attr(out, "links") <- links; attr(out, "reference") <- reference
  attr(out, "DTF") <- dtf
  class(out) <- c("irt_dif_irt", "data.frame")
  out
}

# reported parameters of one item with their covariance matrix (unidimensional)
.item_par_vcov <- function(fit, j) {
  it <- fit$spec$items[[j]]
  f <- function(p) {
    r <- .report_item(it, p, 1, "theta")
    if (it$fam == "gpcm" && it$model != "rsm" && any(grepl("^b[0-9]", names(r)))) r <- r[names(r) != "b"]
    r
  }
  est <- f(fit$par)
  V <- matrix(NA, length(est), length(est))
  if (!is.null(fit$vcov)) {
    Jm <- .jacobian(f, fit$par, it$idx)
    V <- Jm %*% fit$vcov %*% t(Jm)
  }
  dimnames(V) <- list(names(est), names(est))
  list(est = est, V = V)
}

.dif_irt_stats <- function(fit_r, fit_f, lk) {
  A <- lk$A; B <- lk$B
  bank_r <- item_bank(fit_r)
  bank_f <- lk$bank
  th <- seq(-6, 6, length.out = 601); dth <- th[2] - th[1]
  wf <- stats::dnorm(th, bank_f$mu, bank_f$sigma); wf <- wf / sum(wf)
  rows <- lapply(seq_len(ncol(fit_r$data)), function(j) {
    r <- .item_par_vcov(fit_r, j); f <- .item_par_vcov(fit_f, j)
    tr <- ifelse(grepl("^b", names(f$est)), A, ifelse(names(f$est) == "a", 1 / A, 1))
    shift <- ifelse(grepl("^b", names(f$est)), B, 0)
    fe <- f$est * tr + shift
    fV <- diag(tr, length(tr)) %*% f$V %*% diag(tr, length(tr))
    keep <- names(r$est)
    dlt <- fe[keep] - r$est[keep]
    S <- r$V + fV
    chi <- tryCatch(as.numeric(t(dlt) %*% solve(S, dlt)), error = function(e) NA_real_)
    bterms <- grepl("^b", keep)
    er <- .be(bank_r$items[[j]], th); ef <- .be(bank_f$items[[j]], th)
    data.frame(item = colnames(fit_r$data)[j], chisq = chi, df = length(dlt),
               p = stats::pchisq(chi, length(dlt), lower.tail = FALSE),
               delta_b = mean(dlt[bterms]),
               signed_area = sum(er - ef) * dth, unsigned_area = sum(abs(er - ef)) * dth,
               sids = sum(wf * (ef - er)), uids = sum(wf * abs(ef - er)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' @export
print.irt_dif_irt <- function(x, digits = 3, ...) {
  y <- x; class(y) <- "data.frame"
  cat("IRT-based DIF (Lord's Wald test, Raju areas) | reference group: ",
      attr(x, "reference"), "\n", sep = "")
  for (fg in names(attr(x, "links"))) {
    lk <- attr(x, "links")[[fg]]
    cat(sprintf("  %s: %d anchor(s); theta_ref = %.3f * theta_%s + %.3f; DTF = %+.3f score points\n",
                fg, length(lk$anchors), lk$A, fg, lk$B, attr(x, "DTF")[fg]))
  }
  cat("  ", sum(y$flag), " item(s) flagged\n\n", sep = "")
  num <- vapply(y, is.numeric, logical(1)); y[num] <- lapply(y[num], round, digits)
  print(y, row.names = FALSE)
  invisible(x)
}

#' Plot group-specific item response functions from an IRT DIF analysis
#'
#' Draws the reference and focal expected item score curves on the common
#' (linked) scale for each item, shading the area between them, with the
#' Lord's test and SIDS printed in each panel.
#'
#' @param x Output of [dif_irt()].
#' @param items Items to plot. Defaults to the flagged items (or all items if
#'   none are flagged).
#' @param focal Focal group to show (defaults to the first).
#' @param theta_range Trait range.
#' @param ncol Number of facet columns.
#' @export
plot_dif_curves <- function(x, items = NULL, focal = NULL, theta_range = c(-4, 4), ncol = NULL) {
  if (!inherits(x, "irt_dif_irt")) stop("'x' must come from dif_irt().", call. = FALSE)
  ref <- attr(x, "reference")
  if (is.null(focal)) focal <- names(attr(x, "links"))[1]
  lk <- attr(x, "links")[[focal]]
  br <- item_bank(attr(x, "fits")[[ref]]); bf <- lk$bank
  xs <- x[x$focal == focal, ]
  if (is.null(items)) items <- if (any(xs$flag)) xs$item[xs$flag] else xs$item
  th <- .theta_seq(theta_range)
  df <- do.call(rbind, lapply(items, function(j) data.frame(
    item = j, theta = th, ref = .be(br$items[[j]], th), foc = .be(bf$items[[j]], th))))
  lab <- xs[match(items, xs$item), ]
  lab$text <- sprintf("chi2(%d) = %.1f, p = %s\nSIDS = %+.2f", lab$df, lab$chisq,
                      format.pval(lab$p_adj, digits = 2, eps = 0.001), lab$sids)
  df$item <- factor(df$item, levels = items); lab$item <- factor(lab$item, levels = items)
  long <- rbind(data.frame(df[, c("item", "theta")], value = df$ref, group = ref),
                data.frame(df[, c("item", "theta")], value = df$foc, group = focal))
  long$group <- factor(long$group, levels = c(ref, focal))
  cols <- stats::setNames(c(.pal[["blue"]], .pal[["orange"]]), c(ref, focal))
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = df, ggplot2::aes(.data$theta, ymin = pmin(.data$ref, .data$foc),
                                                 ymax = pmax(.data$ref, .data$foc)),
                         fill = .pal[["orange"]], alpha = 0.15) +
    ggplot2::geom_line(data = long, ggplot2::aes(.data$theta, .data$value, colour = .data$group),
                       linewidth = 0.8) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = theta_range[1], y = Inf, label = .data$text),
                       hjust = 0, vjust = 1.2, size = 2.8, colour = .ink[["secondary"]], lineheight = 0.9) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::facet_wrap(~item, ncol = ncol, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.35))) +
    ggplot2::labs(x = expression(theta ~ "(reference scale)"), y = "Expected item score",
                  title = "Group item response functions",
                  subtitle = sprintf("%s vs %s after %s linking on %d anchor item(s); shading = DIF",
                                     focal, ref, lk$method, length(lk$anchors))) +
    theme_irt()
}
