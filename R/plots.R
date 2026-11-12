# Graphics -------------------------------------------------------------------------
# All plotting functions return ggplot objects, so they can be extended with the
# usual `+ theme(...)`, `+ labs(...)`, `+ facet_wrap(...)` grammar.

.pal <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
          magenta = "#e87ba4", green = "#008300", violet = "#6250d6", red = "#e34948")
.ink <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#8a8984", grid = "#e7e6e2",
          band = "#eef3fb", neutral = "#d9d8d3")

#' Categorical palette used by IRTistic
#' @param n Number of colours (at most 8).
#' @export
irt_palette <- function(n = 8) {
  if (n > 8) stop("At most 8 categorical colours; group the rest or use facets.", call. = FALSE)
  unname(.pal[seq_len(n)])
}

.seq_ramp <- function(n) grDevices::colorRampPalette(c("#a9c8f0", "#2a78d6", "#0d2f5c"))(n)
.div_ramp <- c(low = "#2a78d6", mid = "#f2f1ee", high = "#e34948")

#' A clean ggplot2 theme for measurement graphics
#' @param base_size Base font size.
#' @export
theme_irt <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = .ink[["grid"]], linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = .ink[["secondary"]]),
      axis.title = ggplot2::element_text(colour = .ink[["secondary"]]),
      plot.title = ggplot2::element_text(colour = .ink[["primary"]], face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = .ink[["secondary"]]),
      plot.caption = ggplot2::element_text(colour = .ink[["muted"]]),
      strip.text = ggplot2::element_text(colour = .ink[["primary"]], face = "bold", hjust = 0),
      legend.position = "top", legend.justification = "left",
      legend.title = ggplot2::element_text(colour = .ink[["secondary"]]),
      plot.title.position = "plot")
}

.theta_seq <- function(range = c(-4, 4), n = 161) seq(range[1], range[2], length.out = n)

.check_uni <- function(fit, dim) {
  if (fit$spec$D > 1 && is.null(dim))
    stop("Multidimensional model: choose a dimension with 'dim' (others held at 0).", call. = FALSE)
  if (is.null(dim)) 1 else if (is.character(dim)) match(dim, colnames(fit$Q)) else dim
}

.grid_matrix <- function(fit, th, d) {
  m <- matrix(0, length(th), fit$spec$D); m[, d] <- th; m
}

.select_items <- function(fit, items) {
  nm <- colnames(fit$data)
  if (is.null(items)) return(seq_along(nm))
  if (is.character(items)) { j <- match(items, nm); if (anyNA(j)) stop("Unknown item(s).", call. = FALSE); j } else items
}

#' Item characteristic and category probability curves
#'
#' @param fit An `irtfit` object.
#' @param items Items to show (names or indices); all by default.
#' @param type `"category"` for category probability curves or
#'   `"expected"` for expected item scores.
#' @param theta_range Range of the latent trait.
#' @param dim For multidimensional models, the dimension to vary.
#' @param ncol Number of facet columns.
#' @return A ggplot object faceted by item.
#' @export
plot_icc <- function(fit, items = NULL, type = c("category", "expected"),
                     theta_range = c(-4, 4), dim = NULL, ncol = NULL) {
  type <- match.arg(type)
  d <- .check_uni(fit, dim)
  js <- .select_items(fit, items)
  th <- .theta_seq(theta_range)
  tm <- .grid_matrix(fit, th, d)
  rows <- lapply(js, function(j) {
    it <- fit$spec$items[[j]]
    P <- .item_probs(it$fam, .decode_item(it, fit$par, fit$spec$D), tm)
    lab <- fit$labels[[j]]
    if (type == "expected") data.frame(item = colnames(fit$data)[j], theta = th,
                                        value = as.vector(P %*% lab), category = "E")
    else data.frame(item = colnames(fit$data)[j], theta = rep(th, it$K),
                    value = as.vector(P), category = factor(rep(lab, each = length(th))))
  })
  df <- do.call(rbind, rows)
  df$item <- factor(df$item, levels = colnames(fit$data)[js])
  dich <- all(vapply(fit$spec$items[js], function(it) it$K == 2, logical(1)))
  if (type == "category" && dich) { df <- df[df$category == levels(df$category)[2], ]; type <- "icc" }
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$theta, .data$value))
  if (type == "category") {
    K <- nlevels(df$category)
    p <- p + ggplot2::geom_line(ggplot2::aes(colour = .data$category), linewidth = 0.7) +
      ggplot2::scale_colour_manual(values = .seq_ramp(K), name = "Category") +
      ggplot2::labs(y = "Probability")
  } else {
    p <- p + ggplot2::geom_line(colour = .pal[["blue"]], linewidth = 0.7) +
      ggplot2::labs(y = if (type == "icc") "P(correct)" else "Expected score")
  }
  p + ggplot2::facet_wrap(~item, ncol = ncol) +
    ggplot2::labs(x = expression(theta), title = switch(type, category = "Category probability curves",
      icc = "Item characteristic curves", expected = "Expected item score curves")) +
    theme_irt()
}

.item_info <- function(fit, j, tm, d, h = 1e-4) {
  it <- fit$spec$items[[j]]; pp <- .decode_item(it, fit$par, fit$spec$D)
  P <- .item_probs(it$fam, pp, tm)
  e <- matrix(0, nrow(tm), ncol(tm)); e[, d] <- h
  dP <- (.item_probs(it$fam, pp, tm + e) - .item_probs(it$fam, pp, tm - e)) / (2 * h)
  rowSums(dP^2 / P)
}

#' Information and standard error functions
#'
#' @inheritParams plot_icc
#' @param type `"test"` (test information and conditional SE, in two
#'   stacked panels sharing the trait axis) or `"item"` (item information
#'   curves, faceted).
#' @export
plot_info <- function(fit, type = c("test", "item"), items = NULL,
                      theta_range = c(-4, 4), dim = NULL, ncol = NULL) {
  type <- match.arg(type)
  d <- .check_uni(fit, dim)
  js <- .select_items(fit, items)
  th <- .theta_seq(theta_range); tm <- .grid_matrix(fit, th, d)
  info <- vapply(js, function(j) .item_info(fit, j, tm, d), numeric(length(th)))
  info <- matrix(info, length(th))
  if (type == "item") {
    df <- data.frame(item = factor(rep(colnames(fit$data)[js], each = length(th)),
                                   levels = colnames(fit$data)[js]),
                     theta = th, info = as.vector(info))
    return(ggplot2::ggplot(df, ggplot2::aes(.data$theta, .data$info)) +
             ggplot2::geom_area(fill = .pal[["blue"]], alpha = 0.15) +
             ggplot2::geom_line(colour = .pal[["blue"]], linewidth = 0.7) +
             ggplot2::facet_wrap(~item, ncol = ncol) +
             ggplot2::labs(x = expression(theta), y = "Information",
                           title = "Item information functions") + theme_irt())
  }
  ti <- rowSums(info)
  df <- rbind(data.frame(panel = "Test information", theta = th, value = ti),
              data.frame(panel = "Standard error", theta = th, value = 1 / sqrt(ti)))
  df$panel <- factor(df$panel, levels = c("Test information", "Standard error"))
  ggplot2::ggplot(df, ggplot2::aes(.data$theta, .data$value)) +
    ggplot2::geom_line(colour = .pal[["blue"]], linewidth = 0.8) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = expression(theta), y = NULL, title = "Test information and conditional SE",
                  subtitle = sprintf("Peak information %.2f at theta = %.2f", max(ti), th[which.max(ti)])) +
    theme_irt()
}

#' Test characteristic curve
#' @inheritParams plot_icc
#' @export
plot_tcc <- function(fit, theta_range = c(-4, 4), dim = NULL) {
  d <- .check_uni(fit, dim)
  th <- .theta_seq(theta_range); tm <- .grid_matrix(fit, th, d)
  ts <- Reduce(`+`, lapply(seq_len(ncol(fit$data)), function(j) {
    it <- fit$spec$items[[j]]
    .item_probs(it$fam, .decode_item(it, fit$par, fit$spec$D), tm) %*% fit$labels[[j]] }))
  ggplot2::ggplot(data.frame(theta = th, score = as.vector(ts)), ggplot2::aes(.data$theta, .data$score)) +
    ggplot2::geom_line(colour = .pal[["blue"]], linewidth = 0.8) +
    ggplot2::labs(x = expression(theta), y = "Expected total score",
                  title = "Test characteristic curve") + theme_irt()
}

#' Wright map (person-item map)
#'
#' Persons (histogram, left) and items (right) on the same logit scale. For
#' polytomous items the item thresholds are shown as smaller points along
#' each item's row.
#'
#' @param fit An `irtfit` object.
#' @param person_method Person measure to plot (`"EAP"`, `"WLE"`, ...).
#' @param dim For multidimensional models, which dimension to map.
#' @param bins Number of histogram bins.
#' @export
plot_wright <- function(fit, person_method = NULL, dim = NULL, bins = 30) {
  d <- .check_uni(fit, dim)
  if (is.null(person_method)) person_method <- .default_theta(fit)
  sc <- irt_scores(fit, person_method)
  pth <- sc[[if (fit$spec$D == 1) "theta" else colnames(fit$Q)[d]]]
  pth <- pth[is.finite(pth)]
  on <- which(fit$Q[, d] != 0)
  items <- colnames(fit$data)[on]
  locs <- lapply(on, function(j) {
    it <- fit$spec$items[[j]]; pp <- .decode_item(it, fit$par, fit$spec$D)
    a <- pp$a[d]
    th <- -pp$d / a
    list(loc = mean(th), th = if (it$K > 2) th else numeric(0))
  })
  loc <- vapply(locs, `[[`, numeric(1), "loc")
  br <- pretty(range(c(pth, loc, unlist(lapply(locs, `[[`, "th")))), n = bins)
  h <- graphics::hist(pth, breaks = br, plot = FALSE)
  hd <- data.frame(ymin = h$breaks[-length(h$breaks)], ymax = h$breaks[-1],
                   xmin = -h$counts / max(h$counts), xmax = 0)
  ord <- order(loc)
  ni <- length(ord)
  xi <- 0.15 + 0.85 * (seq_len(ni) - 1) / max(ni - 1, 1)
  idf <- data.frame(item = items[ord], loc = loc[ord], x = xi)
  tdf <- do.call(rbind, lapply(seq_along(ord), function(k) {
    t <- locs[[ord[k]]]$th; if (!length(t)) return(NULL)
    data.frame(x = xi[k], y = t) }))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = hd, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                               ymin = .data$ymin, ymax = .data$ymax),
                       fill = .pal[["blue"]], colour = "white", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = .ink[["neutral"]])
  if (!is.null(tdf)) p <- p + ggplot2::geom_segment(
    data = do.call(rbind, lapply(split(tdf, tdf$x), function(z) data.frame(x = z$x[1], lo = min(z$y), hi = max(z$y)))),
    ggplot2::aes(x = .data$x, xend = .data$x, y = .data$lo, yend = .data$hi), colour = .ink[["neutral"]], linewidth = 0.6) +
    ggplot2::geom_segment(data = tdf, ggplot2::aes(x = .data$x - 0.025, xend = .data$x + 0.025,
                                                   y = .data$y, yend = .data$y),
                          colour = .pal[["orange"]], linewidth = 0.7)
  p + ggplot2::geom_point(data = idf, ggplot2::aes(.data$x, .data$loc), colour = .pal[["orange"]], size = 2.6) +
    ggplot2::geom_text(data = idf, ggplot2::aes(.data$x, .data$loc, label = .data$item),
                       hjust = -0.25, size = 3, colour = .ink[["primary"]]) +
    ggplot2::scale_x_continuous(breaks = c(-0.5, 0.55), labels = c("Persons", "Items"),
                                limits = c(-1.05, 1.25)) +
    ggplot2::labs(x = NULL, y = "Measure (logits)", title = "Wright map",
                  subtitle = sprintf("%d persons (%s) and %d items on one scale", length(pth),
                                     person_method, ni)) +
    theme_irt() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                                 axis.text.x = ggplot2::element_text(face = "bold"))
}

#' Fit bubble chart
#'
#' Plots infit (or outfit) mean-square against measure, with bubble area
#' proportional to the standard error and the productive range shaded.
#' Works for [item_fit()], [person_fit()] and [facet_report()] output (the
#' last is faceted by facet).
#'
#' @param x Output of [item_fit()], [person_fit()] or [facet_report()], or an
#'   `irtfit` (item fit is computed).
#' @param stat `"infit"` or `"outfit"`.
#' @param range Productive mean-square range to shade.
#' @param label Label elements outside the range (`TRUE`), all (`"all"`) or
#'   none (`FALSE`).
#' @param max_labels Maximum number of misfitting elements to label (the
#'   most extreme are chosen).
#' @export
plot_fit <- function(x, stat = c("infit", "outfit"), range = c(0.5, 1.5), label = TRUE,
                     max_labels = 10) {
  stat <- match.arg(stat)
  if (inherits(x, "irtfit")) x <- item_fit(x)
  if (inherits(x, "mfrmfit")) x <- facet_report(x)
  if (inherits(x, "mfrm_report")) {
    df <- do.call(rbind, lapply(names(x$measures), function(f)
      data.frame(facet = f, x$measures[[f]][, c("element", "measure", "se", "infit", "outfit")])))
    df$facet <- factor(df$facet, levels = names(x$measures))
  } else {
    idcol <- intersect(c("item", "person"), names(x))[1]
    mcol <- intersect(c("measure", "theta"), names(x))[1]
    df <- data.frame(element = x[[idcol]], measure = x[[mcol]],
                     se = if ("se" %in% names(x)) x$se else NA, infit = x$infit, outfit = x$outfit)
  }
  df$value <- df[[stat]]
  df <- df[is.finite(df$value) & is.finite(df$measure), ]
  df$out <- df$value < range[1] | df$value > range[2]
  df$size <- if (all(is.na(df$se))) 1 else ifelse(is.na(df$se), stats::median(df$se, na.rm = TRUE), df$se)
  # facets differ hugely in precision: scale bubbles within each facet
  if ("facet" %in% names(df)) df$size <- stats::ave(df$size, df$facet, FUN = function(v) v / max(v))
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$measure, .data$value)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = range[1], ymax = range[2],
                      fill = .ink[["band"]]) +
    ggplot2::geom_hline(yintercept = 1, colour = .ink[["muted"]], linetype = 2, linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = .data$size, colour = .data$out),
                        alpha = 0.75, stroke = 0.6, shape = 21, fill = NA) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = .pal[["blue"]], `TRUE` = .pal[["red"]]),
                                 labels = c(`FALSE` = "Within range", `TRUE` = "Outside range"),
                                 name = NULL, drop = FALSE) +
    ggplot2::scale_size(range = c(1.8, 6), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
    ggplot2::labs(x = "Measure (logits)", y = paste(tools_title(stat), "mean-square"),
                  title = paste(tools_title(stat), "fit by measure"),
                  subtitle = sprintf("Shaded band: %.1f-%.1f; bubble size proportional to SE",
                                     range[1], range[2])) +
    theme_irt()
  lab <- if (identical(label, "all")) df else if (isTRUE(label)) df[df$out, ] else df[0, ]
  if (!identical(label, "all") && nrow(lab) > max_labels)
    lab <- lab[order(-abs(log(lab$value)))[seq_len(max_labels)], ]
  if (nrow(lab)) p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$element),
                                             size = 3, vjust = -1, colour = .ink[["primary"]])
  if ("facet" %in% names(df)) p <- p + ggplot2::facet_wrap(~facet, scales = "free_x")
  p
}

tools_title <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))

#' Local dependence heatmap (Yen's Q3)
#' @param x Output of [local_dependence()] or an `irtfit`.
#' @export
plot_ld <- function(x) {
  if (inherits(x, "irtfit")) x <- local_dependence(x)
  q <- x$Q3; diag(q) <- NA
  nm <- colnames(q)
  df <- data.frame(item1 = factor(rep(nm, each = length(nm)), levels = nm),
                   item2 = factor(rep(nm, length(nm)), levels = rev(nm)), Q3 = as.vector(q))
  ggplot2::ggplot(df, ggplot2::aes(.data$item1, .data$item2, fill = .data$Q3 - x$mean_Q3)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradient2(low = .div_ramp[["low"]], mid = .div_ramp[["mid"]],
                                  high = .div_ramp[["high"]], midpoint = 0, na.value = "white",
                                  name = "Q3 - mean") +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, title = "Residual correlations (Q3)",
                  subtitle = sprintf("Centred at the average Q3 (%.3f); red = possible local dependence", x$mean_Q3)) +
    theme_irt() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                                 panel.grid = ggplot2::element_blank(), legend.position = "right")
}

#' Residual PCA contrast plot
#' @param x Output of [residual_pca()] or an `irtfit`.
#' @param contrast Which contrast to plot.
#' @export
plot_pca <- function(x, contrast = 1) {
  if (inherits(x, "irtfit")) x <- residual_pca(x)
  df <- x$loadings; df$loading <- df[[paste0("contrast", contrast)]]
  if (all(is.na(df$measure))) df$measure <- seq_len(nrow(df))
  ggplot2::ggplot(df, ggplot2::aes(.data$measure, .data$loading)) +
    ggplot2::geom_hline(yintercept = 0, colour = .ink[["neutral"]]) +
    ggplot2::geom_point(colour = .pal[["blue"]], size = 2.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$item), size = 3, vjust = -0.9, colour = .ink[["primary"]]) +
    ggplot2::labs(x = "Item measure", y = paste("Loading on contrast", contrast),
                  title = "Standardised residual contrast",
                  subtitle = sprintf("Eigenvalue %.2f (%.1f%% of residual variance)",
                                     x$eigenvalues[contrast], 100 * x$proportion[contrast])) +
    theme_irt()
}

#' DIF effect-size plot
#'
#' Lollipop chart of DIF effect sizes with classification thresholds. For
#' [dif_test()] output this shows the MH delta or Delta R-squared; for
#' [dif_irt()] output it shows SIDS (signed expected item-score difference,
#' focal minus reference) with flagged items highlighted, one panel per
#' focal group.
#' @param x Output of [dif_test()] or [dif_irt()].
#' @export
plot_dif <- function(x) {
  if (inherits(x, "irt_dif_irt")) return(.plot_dif_irt(x))
  df <- as.data.frame(x); class(df) <- "data.frame"
  mh <- attr(x, "method") == "mh"
  df$effect <- if (mh) df$delta_MH else df$delta_R2
  df$item <- factor(df$item, levels = df$item[order(df$effect)])
  df$class <- factor(df$class, levels = c("A", "B", "C"))
  cut <- if (mh) c(-1.5, -1, 1, 1.5) else c(0.035, 0.07)
  ggplot2::ggplot(df, ggplot2::aes(.data$effect, .data$item)) +
    ggplot2::geom_vline(xintercept = cut, colour = .ink[["neutral"]], linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, colour = .ink[["muted"]]) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$effect, yend = .data$item), colour = .ink[["neutral"]]) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$class, shape = .data$class), size = 2.8) +
    ggplot2::scale_colour_manual(values = c(A = .pal[["blue"]], B = .pal[["yellow"]], C = .pal[["red"]]),
                                 drop = FALSE, name = "DIF class") +
    ggplot2::scale_shape_manual(values = c(A = 16, B = 17, C = 15), drop = FALSE, name = "DIF class") +
    ggplot2::labs(x = if (mh) "MH delta (ETS)" else expression(Delta * R^2 ~ "(Nagelkerke)"), y = NULL,
                  title = "Differential item functioning",
                  subtitle = if (mh) "Dashed lines: |delta| = 1 (B) and 1.5 (C)" else "Dashed lines: 0.035 (B) and 0.07 (C)") +
    theme_irt()
}

#' FACETS-style vertical ruler (variable map)
#'
#' Shows persons and every facet's elements on one logit ruler, plus a
#' column with the rating scale's expected-score zones, like FACETS'
#' Table 6.
#'
#' @param fit An `mfrmfit` object.
#' @param person_method Person measure used for the person column.
#' @export
plot_facets <- function(fit, person_method = "WLE") {
  rep <- facet_report(fit, person_method)
  pm <- rep$measures[[fit$person]]$measure
  cols <- c(fit$person, fit$facets, "Scale")
  pp <- .mfrm_unpack(fit$ctx, fit$par)
  rng <- range(c(pm, unlist(pp$beta)), finite = TRUE) + c(-0.5, 0.5)
  # person histogram
  br <- pretty(rng, 30)
  h <- graphics::hist(pmin(pmax(pm, min(br)), max(br)), breaks = br, plot = FALSE)
  hd <- data.frame(column = factor(fit$person, levels = cols),
                   ymin = h$breaks[-length(h$breaks)], ymax = h$breaks[-1],
                   xmin = 0, xmax = h$counts / max(h$counts))
  el <- do.call(rbind, lapply(fit$facets, function(f) {
    m <- pp$beta[[f]]; o <- order(m)
    data.frame(column = factor(f, levels = cols), label = fit$levels[[f]][o], y = m[o],
               x = 0.08 + 0.5 * ((seq_along(o) - 1) %% 2))
  }))
  # expected score zones at the average step structure (half-point thresholds)
  tau <- colMeans(pp$tau); k <- 0:(fit$ctx$K - 1); ct <- c(0, cumsum(tau))
  ys <- seq(rng[1] - 1, rng[2] + 1, length.out = 801)
  es <- vapply(ys, function(y) { num <- k * y - ct; p <- exp(num - max(num)); sum(p * k / sum(p)) }, numeric(1))
  cuts <- vapply(seq_len(fit$ctx$K - 1), function(v) ys[which.min(abs(es - (v - 0.5)))], numeric(1))
  zb <- c(rng[1], pmin(pmax(cuts, rng[1]), rng[2]), rng[2])
  zd <- data.frame(column = factor("Scale", levels = cols), ymin = zb[-length(zb)], ymax = zb[-1],
                   cat = factor(fit$categories, levels = fit$categories))
  zd$light <- ifelse(seq_len(nrow(zd)) <= ceiling(nrow(zd) / 2), .ink[["primary"]], "white")
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = zd, ggplot2::aes(xmin = 0, xmax = 1, ymin = .data$ymin, ymax = .data$ymax,
                                               fill = .data$cat), colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(data = zd, ggplot2::aes(x = 0.5, y = (.data$ymin + .data$ymax) / 2, label = .data$cat,
                                               colour = .data$light), fontface = "bold", size = 3.4) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_manual(values = .seq_ramp(fit$ctx$K), guide = "none") +
    ggplot2::geom_rect(data = hd, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                               ymin = .data$ymin, ymax = .data$ymax),
                       fill = .pal[["blue"]], colour = "white", linewidth = 0.3) +
    ggplot2::geom_point(data = el, ggplot2::aes(x = .data$x, y = .data$y), colour = .pal[["orange"]], size = 2.2) +
    ggplot2::geom_text(data = el, ggplot2::aes(x = .data$x + 0.05, y = .data$y, label = .data$label),
                       hjust = 0, size = 3, colour = .ink[["primary"]]) +
    ggplot2::facet_grid(. ~ column, drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1.05), breaks = NULL) +
    ggplot2::coord_cartesian(ylim = rng) +
    ggplot2::labs(x = NULL, y = "Measure (logits)", title = "Variable map",
                  subtitle = "Higher facet measures = harder to score highly on (more difficult / more severe)") +
    theme_irt() + ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0.5),
                                 panel.spacing.x = ggplot2::unit(0.4, "lines"))
}

#' Rating scale category probability curves for a many-facet model
#' @param fit An `mfrmfit` object.
#' @param range Range of the measure relative to the combined facet location.
#' @export
plot_categories <- function(fit, range = c(-5, 5)) {
  pp <- .mfrm_unpack(fit$ctx, fit$par)
  th <- .theta_seq(range); K <- fit$ctx$K; k <- 0:(K - 1)
  groups <- if (is.null(fit$step_facet)) "Rating scale" else fit$levels[[fit$step_facet]]
  df <- do.call(rbind, lapply(seq_len(fit$ctx$S), function(s) {
    ct <- c(0, cumsum(pp$tau[s, ]))
    num <- outer(th, k) - matrix(ct, length(th), K, byrow = TRUE)
    P <- exp(num - apply(num, 1, max)); P <- P / rowSums(P)
    data.frame(group = groups[s], x = th, p = as.vector(P),
               category = factor(rep(fit$categories, each = length(th)), levels = fit$categories))
  }))
  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$p, colour = .data$category)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = .seq_ramp(K), name = "Category") +
    ggplot2::facet_wrap(~group) +
    ggplot2::labs(x = "Measure relative to facet difficulty (logits)", y = "Probability",
                  title = "Category probability curves") + theme_irt()
}

#' Bias / interaction heatmap
#' @param x Output of [mfrm_bias()] with two facets.
#' @param stat `"t"` (default) or `"bias"`.
#' @export
plot_bias <- function(x, stat = c("t", "bias")) {
  stat <- match.arg(stat)
  f <- attr(x, "facets")
  if (length(f) != 2) stop("plot_bias() needs a two-facet bias analysis.", call. = FALSE)
  df <- as.data.frame(x); class(df) <- "data.frame"
  df$value <- df[[stat]]
  df$lab <- ifelse(abs(df$t) >= 2, sprintf("%.1f", df$value), "")
  lim <- max(abs(df$value), if (stat == "t") 2.5 else 0.5)
  ggplot2::ggplot(df, ggplot2::aes(.data[[f[2]]], .data[[f[1]]], fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data$lab), size = 3, colour = .ink[["primary"]]) +
    ggplot2::scale_fill_gradient2(low = .div_ramp[["low"]], mid = .div_ramp[["mid"]],
                                  high = .div_ramp[["high"]], midpoint = 0, limits = c(-lim, lim),
                                  name = if (stat == "t") "Bias t" else "Bias (logits)") +
    ggplot2::labs(title = paste("Bias analysis:", f[1], "x", f[2]),
                  subtitle = "Red = higher scores than expected; labelled cells have |t| >= 2") +
    theme_irt() + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")
}

#' Plot a fitted IRT model
#'
#' @param x An `irtfit` object.
#' @param type One of `"wright"`, `"icc"`, `"expected"`, `"info"`,
#'   `"iteminfo"`, `"tcc"`, `"fit"`, `"ld"`, `"pca"`.
#' @param ... Passed to the underlying `plot_*()` function.
#' @return A ggplot object.
#' @export
plot.irtfit <- function(x, type = c("wright", "icc", "expected", "info", "iteminfo",
                                    "tcc", "fit", "ld", "pca"), ...) {
  type <- match.arg(type)
  switch(type,
         wright = plot_wright(x, ...), icc = plot_icc(x, ...),
         expected = plot_icc(x, type = "expected", ...), info = plot_info(x, ...),
         iteminfo = plot_info(x, type = "item", ...), tcc = plot_tcc(x, ...),
         fit = plot_fit(item_fit(x), ...), ld = plot_ld(x), pca = plot_pca(x, ...))
}

#' Plot a fitted many-facet Rasch model
#'
#' @param x An `mfrmfit` object.
#' @param type `"facets"` (vertical ruler), `"categories"`, `"fit"` or `"bias"`.
#' @param facets For `type = "bias"`, the two facets to cross.
#' @param ... Passed to the underlying `plot_*()` function.
#' @export
plot.mfrmfit <- function(x, type = c("facets", "categories", "fit", "bias"), facets = NULL, ...) {
  type <- match.arg(type)
  switch(type,
         facets = plot_facets(x, ...), categories = plot_categories(x, ...),
         fit = plot_fit(facet_report(x), ...),
         bias = plot_bias(mfrm_bias(x, if (is.null(facets)) x$facets[1:2] else facets), ...))
}

# CAT graphics ----------------------------------------------------------------------

#' Trace of an adaptive test
#'
#' Shows how the trait estimate and its 95% interval evolve item by item,
#' with the true value (for simulations) and any classification cut scores.
#' Points are coloured by whether the response was in the top or bottom
#' category.
#'
#' @param x A `cat_session` or `cat_sim`.
#' @param person For a `cat_sim`, the simulee to show.
#' @export
plot_cat_trace <- function(x, person = 1) {
  if (inherits(x, "cat_sim")) {
    h <- x$history[x$history$person == person, ]
    truth <- x$results$true_theta[x$results$person == person]
    d <- x$design
  } else if (inherits(x, "cat_session")) {
    h <- x$history; truth <- NA; d <- x$design
  } else stop("'x' must be a cat_session or cat_sim.", call. = FALSE)
  if (!nrow(h)) stop("No items administered yet.", call. = FALSE)
  maxs <- vapply(h$item, function(i) max(d$bank$items[[i]]$scores), numeric(1))
  mins <- vapply(h$item, function(i) min(d$bank$items[[i]]$scores), numeric(1))
  h$resp <- factor(ifelse(h$response == maxs, "Highest category", ifelse(h$response == mins, "Lowest category", "Middle")),
                   levels = c("Highest category", "Middle", "Lowest category"))
  p <- ggplot2::ggplot(h, ggplot2::aes(.data$step, .data$theta)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$theta - 1.96 * .data$se, ymax = .data$theta + 1.96 * .data$se),
                         fill = .pal[["blue"]], alpha = 0.15)
  if (!is.null(d$cut)) p <- p + ggplot2::geom_hline(yintercept = d$cut, colour = .pal[["red"]], linetype = 2)
  if (is.finite(truth)) p <- p + ggplot2::geom_hline(yintercept = truth, colour = .ink[["muted"]], linewidth = 0.5,
                                                     linetype = "dashed")
  p + ggplot2::geom_line(colour = .pal[["blue"]], linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$resp), shape = 21, size = 2.6, colour = "white", stroke = 0.6) +
    ggplot2::scale_fill_manual(values = c(`Highest category` = .pal[["aqua"]], Middle = .pal[["yellow"]],
                                          `Lowest category` = .pal[["orange"]]), name = NULL, drop = TRUE) +
    ggplot2::scale_x_continuous(breaks = function(l) unique(round(pretty(l)))) +
    ggplot2::labs(x = "Items administered", y = expression(hat(theta) ~ "with 95% interval"),
                  title = "Adaptive test trace",
                  subtitle = paste0(sprintf("Final estimate %.2f (SE %.2f) after %d items", utils::tail(h$theta, 1),
                                            utils::tail(h$se, 1), nrow(h)),
                                    if (is.finite(truth)) sprintf("; dashed line = true value (%.2f)", truth) else "")) +
    theme_irt()
}

#' Item exposure in a CAT simulation
#'
#' Exposure rate of every bank item plotted against its difficulty, with
#' bubble size proportional to discrimination and an optional target
#' maximum exposure line.
#'
#' @param x A `cat_sim`.
#' @param r_max Optional target maximum exposure rate to draw.
#' @export
plot_cat_exposure <- function(x, r_max = NULL) {
  e <- x$exposure
  e$status <- ifelse(e$exposure == 0, "Never used", if (is.null(r_max)) "Used" else
    ifelse(e$exposure > r_max, "Over target", "Used"))
  cols <- c(Used = .pal[["blue"]], `Over target` = .pal[["red"]], `Never used` = .ink[["neutral"]])
  p <- ggplot2::ggplot(e, ggplot2::aes(.data$b, .data$exposure))
  if (!is.null(r_max)) p <- p + ggplot2::geom_hline(yintercept = r_max, colour = .pal[["red"]], linetype = 2)
  p + ggplot2::geom_point(ggplot2::aes(size = .data$a, colour = .data$status), alpha = 0.7) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_size(range = c(1, 5), name = "Discrimination") +
    ggplot2::labs(x = "Item difficulty (mean threshold)", y = "Exposure rate",
                  title = "Item exposure",
                  subtitle = sprintf("Max %.2f; %.0f%% of the bank never administered",
                                     max(e$exposure), 100 * mean(e$exposure == 0))) +
    theme_irt()
}

#' Conditional precision and test length in a CAT simulation
#'
#' Small multiples of bias, RMSE, mean SE and mean test length across bins of
#' the true trait.
#'
#' @param x A `cat_sim` with known true trait values.
#' @param bins Number of bins.
#' @export
plot_cat_conditional <- function(x, bins = 10) {
  cc <- cat_conditional(x, bins)
  long <- rbind(data.frame(theta = cc$theta, value = cc$bias, stat = "Bias"),
                data.frame(theta = cc$theta, value = cc$rmse, stat = "RMSE"),
                data.frame(theta = cc$theta, value = cc$mean_se, stat = "Mean SE"),
                data.frame(theta = cc$theta, value = cc$mean_length, stat = "Mean test length"))
  long$stat <- factor(long$stat, levels = c("Bias", "RMSE", "Mean SE", "Mean test length"))
  ref <- data.frame(stat = factor("Bias", levels = levels(long$stat)), y = 0)
  ggplot2::ggplot(long, ggplot2::aes(.data$theta, .data$value)) +
    ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = .data$y), colour = .ink[["muted"]]) +
    ggplot2::geom_line(colour = .pal[["blue"]], linewidth = 0.8) +
    ggplot2::geom_point(colour = .pal[["blue"]], size = 2) +
    ggplot2::facet_wrap(~stat, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = expression("True" ~ theta ~ "(bin mean)"), y = NULL,
                  title = "Conditional CAT performance",
                  subtitle = sprintf("%d simulees in %d equal-count bins", nrow(x$results), bins)) +
    theme_irt()
}

#' @export
plot.cat_sim <- function(x, type = c("conditional", "exposure", "trace"), ...) {
  type <- match.arg(type)
  switch(type, conditional = plot_cat_conditional(x, ...), exposure = plot_cat_exposure(x, ...),
         trace = plot_cat_trace(x, ...))
}

#' @export
plot.cat_session <- function(x, ...) plot_cat_trace(x, ...)

.plot_dif_irt <- function(x) {
  df <- as.data.frame(x); class(df) <- "data.frame"
  df$status <- factor(ifelse(df$flag, "Flagged", "Not flagged"), levels = c("Flagged", "Not flagged"))
  ord <- stats::aggregate(sids ~ item, df, mean)
  df$item <- factor(df$item, levels = ord$item[order(ord$sids)])
  ggplot2::ggplot(df, ggplot2::aes(.data$sids, .data$item)) +
    ggplot2::geom_vline(xintercept = 0, colour = .ink[["muted"]]) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$sids, yend = .data$item), colour = .ink[["neutral"]]) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status, shape = .data$status), size = 2.8) +
    ggplot2::scale_colour_manual(values = c(Flagged = .pal[["red"]], `Not flagged` = .pal[["blue"]]), name = NULL, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(Flagged = 15, `Not flagged` = 16), name = NULL, drop = FALSE) +
    ggplot2::facet_wrap(~focal, labeller = ggplot2::labeller(focal = function(v) paste("Focal:", v))) +
    ggplot2::labs(x = "SIDS: expected item score, focal minus reference", y = NULL,
                  title = "Differential item functioning (IRT)",
                  subtitle = "Negative = item harder for the focal group at equal trait levels; flagged by Lord's test") +
    theme_irt()
}

#' @export
plot.irt_dif_irt <- function(x, type = c("summary", "curves"), ...) {
  type <- match.arg(type)
  if (type == "summary") plot_dif(x) else plot_dif_curves(x, ...)
}

#' @export
plot.irt_dif <- function(x, ...) plot_dif(x)
