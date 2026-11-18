make_dif_data <- function() {
  set.seed(4)
  d <- sim_irt(1000, 12, "rasch")
  g <- rep(c("ref", "foc"), 500)
  th <- attr(d, "truth")$theta[, 1]
  foc <- g == "foc"
  d$I03[foc] <- rbinom(500, 1, plogis(th[foc] - 1.2))
  d$I08[foc] <- rbinom(500, 1, plogis(th[foc] + 1.0))
  list(d = d, g = g)
}

test_that("MH statistics match stats::mantelhaen.test", {
  z <- make_dif_data()
  res <- dif_test(z$d, z$g, method = "mh")
  X <- as.matrix(z$d); tot <- rowSums(X)
  mt <- stats::mantelhaen.test(factor(z$g, levels = c("ref", "foc")),
                               factor(X[, 3], levels = c(1, 0)), factor(tot))
  expect_equal(res$statistic[3], unname(mt$statistic), tolerance = 1e-8)
  expect_equal(res$alpha_MH[3], 1 / unname(mt$estimate), tolerance = 1e-8)
})

test_that("MH purification removes DIF items from the matching score", {
  z <- make_dif_data()
  res <- dif_test(z$d, z$g, method = "mh", purify = TRUE)
  expect_equal(res$item[!res$anchor], c("I03", "I08"))
  pur <- attr(res, "purification")
  expect_true(pur$converged)
  # purified statistic for a clean item equals MH on the anchor-only total
  X <- as.matrix(z$d); tot <- rowSums(X[, res$anchor])
  mt <- stats::mantelhaen.test(factor(z$g, levels = c("ref", "foc")),
                               factor(X[, 5], levels = c(1, 0)), factor(tot))
  expect_equal(res$statistic[5], unname(mt$statistic), tolerance = 1e-8)
  # studied DIF item is matched on anchors plus itself
  tot3 <- tot + X[, 3]
  mt3 <- stats::mantelhaen.test(factor(z$g, levels = c("ref", "foc")),
                                factor(X[, 3], levels = c(1, 0)), factor(tot3))
  expect_equal(res$statistic[3], unname(mt3$statistic), tolerance = 1e-8)
  # purification reduces contamination of clean items
  raw <- dif_test(z$d, z$g, method = "mh")
  clean <- !res$item %in% c("I03", "I08")
  expect_lt(mean(abs(res$delta_MH[clean])), mean(abs(raw$delta_MH[clean])))
})

test_that("LR purification works with total-score matching", {
  z <- make_dif_data()
  res <- dif_test(z$d, z$g, method = "lr", match = "total", purify = TRUE)
  expect_true(all(c("I03", "I08") %in% res$item[!res$anchor]))
  expect_error(dif_test(irt_fit(z$d, "rasch"), z$g, method = "lr", purify = TRUE), "match = 'total'")
})

test_that("helpful errors for common data problems", {
  d <- sim_irt(200, 5, "rasch")
  expect_error(irt_fit(cbind(d, gi33 = 0), "rasch"), "only one observed value")
  expect_error(irt_fit(cbind(d, id = 1:200), "rasch"), "looks like an ID")
  d2 <- d; d2$I01 <- ifelse(d2$I01 == 1, "yes", "no")
  expect_error(irt_fit(d2, "rasch"), "contain text")
  expect_error(dif_test(d, c(NA, rep("a", 199)), method = "mh"), "missing values")
})
