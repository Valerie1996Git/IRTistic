#' Link two calibrations onto a common scale
#'
#' Finds the linear transformation \eqn{\theta^* = A\theta + B} that places the
#' `from` calibration on the scale of the `to` calibration using common
#' (anchor) items, by the mean-mean, mean-sigma, Haebara or Stocking-Lord
#' method.
#'
#' @param from,to `irtfit` (unidimensional) or `irt_bank` objects. `from` is
#'   transformed onto the scale of `to`.
#' @param anchors Names of the common items. Defaults to all items whose names
#'   appear in both.
#' @param method `"stocking-lord"` (default), `"haebara"`, `"mean-sigma"` or
#'   `"mean-mean"`.
#' @param fix_slope Fix \eqn{A = 1}? Defaults to `TRUE` when every anchor is a
#'   Rasch-family item in both calibrations.
#' @return An object of class `irt_link` with the constants `A` and `B`, the
#'   transformed bank (`bank`), the method and the anchors.
#' @examples
#' set.seed(1)
#' x <- sim_irt(800, 15, "2pl")
#' fx <- irt_fit(x[1:400, ], "2pl", se = "none")
#' fy <- irt_fit(x[401:800, ], "2pl", se = "none")
#' irt_link(fx, fy)
#' @export
irt_link <- function(from, to, anchors = NULL, method = c("stocking-lord", "haebara", "mean-sigma", "mean-mean"),
                     fix_slope = NULL) {
  method <- match.arg(method)
  bx <- if (inherits(from, "irtfit")) item_bank(from) else from
  by <- if (inherits(to, "irtfit")) item_bank(to) else to
  if (is.null(anchors)) anchors <- intersect(names(bx$items), names(by$items))
  anchors <- intersect(anchors, intersect(names(bx$items), names(by$items)))
  if (length(anchors) < 2) stop("At least two common items are needed for linking.", call. = FALSE)
  ax <- bx$items[anchors]; ay <- by$items[anchors]
  rasch <- function(its) all(vapply(its, function(it) it$model %in% c("rasch", "pcm", "rsm"), logical(1)))
  if (is.null(fix_slope)) fix_slope <- rasch(ax) && rasch(ay)
  bvals <- function(its) unlist(lapply(its, function(it) -it$d / it$a))
  avals <- function(its) vapply(its, `[[`, numeric(1), "a")
  Bx <- bvals(ax); By <- bvals(ay)
  if (method == "mean-mean") {
    A <- if (fix_slope) 1 else mean(avals(ax)) / mean(avals(ay))
  } else {
    A <- if (fix_slope) 1 else stats::sd(By) / stats::sd(Bx)
  }
  B <- mean(By) - A * mean(Bx)
  if (method %in% c("stocking-lord", "haebara")) {
    q <- .quad(by$mu, by$sigma, 61)
    Py <- lapply(ay, function(it) .bp(it, q$theta))
    Ey <- Reduce(`+`, lapply(ay, function(it) .be(it, q$theta)))
    crit <- function(p) {
      AA <- if (fix_slope) 1 else exp(p[2])
      tb <- .bank_transform(list(items = ax, mu = 0, sigma = 1), AA, p[1])$items
      if (method == "stocking-lord") {
        Ex <- Reduce(`+`, lapply(tb, function(it) .be(it, q$theta)))
        sum(q$w * (Ey - Ex)^2)
      } else {
        sum(vapply(seq_along(tb), function(k) sum(q$w * rowSums((Py[[k]] - .bp(tb[[k]], q$theta))^2)),
                   numeric(1)))
      }
    }
    start <- if (fix_slope) B else c(B, log(max(A, 0.05)))
    o <- if (fix_slope) stats::optimize(crit, c(B - 5, B + 5)) else stats::optim(start, crit, method = "BFGS")
    if (fix_slope) { B <- o$minimum } else { B <- o$par[1]; A <- exp(o$par[2]) }
  }
  out <- list(A = A, B = B, bank = .bank_transform(bx, A, B), method = method,
              anchors = anchors, fix_slope = fix_slope)
  class(out) <- "irt_link"
  out
}

#' @export
print.irt_link <- function(x, ...) {
  cat(sprintf("<irt_link> %s linking on %d anchor item(s)%s\n", x$method, length(x$anchors),
              if (x$fix_slope) " (slope fixed at 1)" else ""))
  cat(sprintf("  theta(to) = %.4f * theta(from) + %.4f\n", x$A, x$B))
  invisible(x)
}

#' IRT true-score equating
#'
#' Equates number-correct (summed) scores on form X to the scale of form Y
#' through the test characteristic curves, after the two forms have been
#' placed on a common scale with [irt_link()].
#'
#' @param x Form X: an `irtfit`, `irt_bank`, or an `irt_link` whose `bank` is
#'   form X already on Y's scale.
#' @param y Form Y: an `irtfit` or `irt_bank`.
#' @return A data frame with each form-X score, the corresponding trait value
#'   and the equivalent form-Y true score.
#' @examples
#' set.seed(2)
#' d <- sim_irt(1000, 20, "2pl")
#' fx <- irt_fit(d[1:500, 1:15], "2pl", se = "none")     # form X: items 1-15
#' fy <- irt_fit(d[501:1000, 6:20], "2pl", se = "none")  # form Y: items 6-20
#' lk <- irt_link(fx, fy)                                  # anchors: items 6-15
#' head(equate_true_score(lk, fy))
#' @export
equate_true_score <- function(x, y) {
  bx <- if (inherits(x, "irt_link")) x$bank else if (inherits(x, "irtfit")) item_bank(x) else x
  by <- if (inherits(y, "irtfit")) item_bank(y) else y
  tcc <- function(b, th) sum(vapply(b$items, function(it) .be(it, th), numeric(1)))
  lo_x <- sum(vapply(bx$items, function(it) min(it$scores) + it$g * (max(it$scores) - min(it$scores)), numeric(1)))
  minx <- sum(vapply(bx$items, function(it) min(it$scores), numeric(1)))
  maxx <- sum(vapply(bx$items, function(it) max(it$scores), numeric(1)))
  scores <- minx:maxx
  hi_x <- sum(vapply(bx$items, function(it) min(it$scores) +
                      (if (is.null(it$u)) 1 else it$u) * (max(it$scores) - min(it$scores)), numeric(1)))
  th <- vapply(scores, function(s) {
    if (s <= lo_x + 1e-8) return(-Inf)
    if (s >= hi_x - 1e-8) return(Inf)
    stats::uniroot(function(t) tcc(bx, t) - s, c(-15, 15), tol = 1e-8)$root
  }, numeric(1))
  ymin <- sum(vapply(by$items, function(it) min(it$scores) + it$g * (max(it$scores) - min(it$scores)), numeric(1)))
  ymax <- sum(vapply(by$items, function(it) max(it$scores), numeric(1)))
  eq <- vapply(th, function(t) if (t == -Inf) ymin else if (t == Inf) ymax else tcc(by, t), numeric(1))
  # Kolen & Brennan's convention below the guessing floor: linear interpolation
  low <- scores <= lo_x
  if (any(low) && lo_x > minx) eq[low] <- minx + (scores[low] - minx) * (ymin - minx) / (lo_x - minx)
  data.frame(score_x = scores, theta = th, score_y = eq)
}
