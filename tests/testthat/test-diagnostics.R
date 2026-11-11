set.seed(21)
dat <- sim_irt(500, 8, "rasch")
dat$I08 <- rbinom(500, 1, 0.5)   # pure noise item
fit <- irt_fit(dat, "rasch")

test_that("person scoring methods agree in ordering", {
  for (m in c("EAP", "MAP", "ML", "WLE")) {
    s <- irt_scores(fit, m)
    expect_equal(nrow(s), 500)
    # in a Rasch model every method is a monotone function of the raw score
    expect_lt(max(tapply(s$theta, s$raw, function(v) diff(range(v)))), 1e-6)
    expect_true(all(diff(tapply(s$theta, s$raw, mean)) > 0))
  }
  expect_true("extreme" %in% names(irt_scores(fit, "ML")))
})

test_that("item fit flags the noise item", {
  f <- item_fit(fit)
  expect_equal(f$item[which.max(f$infit)], "I08")
  expect_gt(f$infit_z[8], 2)
})

test_that("reliability, Q3 and residual PCA return sensible output", {
  r <- reliability(fit)
  expect_true(all(r$reliability > 0 & r$reliability < 1, na.rm = TRUE))
  ld <- local_dependence(fit)
  expect_equal(dim(ld$Q3), c(8, 8))
  expect_lt(abs(ld$mean_Q3 + 1 / 7), 0.08)
  pc <- residual_pca(fit)
  expect_length(pc$eigenvalues, 8)
})

test_that("DIF detects an item with a planted shift", {
  set.seed(22)
  d <- sim_irt(1200, 8, "rasch")
  g <- rep(c("ref", "foc"), 600)
  th <- attr(d, "truth")$theta[, 1]
  d$I03[g == "foc"] <- rbinom(600, 1, plogis(th[g == "foc"] - attr(d, "truth")$b[3, 1] - 1.2))
  f <- irt_fit(d, "rasch", se = "none")
  mh <- dif_test(f, g, "mh")
  lr <- dif_test(f, g, "lr")
  expect_equal(mh$item[which.max(abs(mh$delta_MH))], "I03")
  expect_equal(lr$item[which.max(lr$statistic)], "I03")
  expect_true(mh$class[3] %in% c("B", "C"))
})
