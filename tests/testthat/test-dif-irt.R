test_that("Lord's test flags a planted DIF item and keeps type I error low", {
  skip_on_cran()
  set.seed(61)
  n <- 1600
  base <- sim_irt(10, 10, "2pl"); tr <- attr(base, "truth")
  tr$a[5] <- 1.3; tr$b[5, 1] <- 0
  g <- rep(c("ref", "foc"), n / 2)
  th <- rnorm(n) - 0.3 * (g == "foc")
  d <- sim_irt(n, a = tr$a, b = tr$b, model = "2pl", theta = matrix(th))
  f <- g == "foc"
  d$I05[f] <- rbinom(sum(f), 1, plogis(1.3 * (th[f] - 0.8)))
  res <- dif_irt(d, g, "2pl", reference = "ref")
  expect_true(res$flag[res$item == "I05"])
  expect_lte(sum(res$flag), 2)
  expect_lt(res$sids[res$item == "I05"], 0)    # harder for focal
  expect_false(res$anchor[res$item == "I05"])  # purified out of the anchor set
  expect_s3_class(plot_dif_curves(res), "ggplot")
  expect_s3_class(plot_dif(res), "ggplot")
})
