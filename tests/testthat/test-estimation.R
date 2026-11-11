test_that("all unidimensional models recover generating parameters", {
  skip_on_cran()
  set.seed(11)
  for (m in c("rasch", "2pl", "grm", "pcm", "gpcm", "rsm")) {
    dat <- sim_irt(800, nitems = 6, model = m, ncat = 4)
    tr <- attr(dat, "truth")
    fit <- irt_fit(dat, model = m)
    expect_true(fit$converged, info = m)
    cf <- coef(fit)
    b <- cf$estimate[grepl("^b[0-9]*$", cf$term) & !(m %in% c("pcm", "gpcm") & cf$term == "b")]
    true_b <- if (m %in% c("rasch", "2pl")) tr$b[, 1] else if (m == "rsm") rowMeans(tr$b) else as.vector(t(tr$b))
    expect_gt(cor(b, true_b), 0.97)
    expect_true(all(is.finite(cf$se)), info = m)
  }
})

test_that("3PL runs with the default guessing prior", {
  skip_on_cran()
  set.seed(12)
  dat <- sim_irt(1000, 8, "3pl")
  fit <- irt_fit(dat, "3pl")
  g <- coef(fit)$estimate[coef(fit)$term == "g"]
  expect_true(all(g > 0 & g < 0.5))
})

test_that("Rasch item difficulties agree with eRm (CML)", {
  skip_if_not_installed("eRm")
  set.seed(13)
  dat <- sim_irt(500, 8, "rasch")
  fit <- irt_fit(dat, "rasch")
  b <- coef(fit)$estimate
  e <- eRm::RM(dat)
  expect_lt(max(abs((b - mean(b)) - (-e$betapar))), 0.02)
  # infit/outfit with ML person measures match eRm
  ifit <- eRm::itemfit(eRm::person.parameter(e))
  ours <- item_fit(fit, "ML")
  expect_equal(ours$infit, unname(ifit$i.infitMSQ), tolerance = 0.02)
  expect_equal(ours$outfit, unname(ifit$i.outfitMSQ), tolerance = 0.02)
})

test_that("multidimensional model recovers the latent correlation", {
  set.seed(14)
  Q <- cbind(rep(1:0, each = 5), rep(0:1, each = 5))
  dat <- sim_irt(1000, 10, "2pl", Q = Q, Sigma = matrix(c(1, .6, .6, 1), 2))
  fit <- irt_fit(dat, "2pl", Q = Q, se = "none")
  expect_equal(fit$Sigma[1, 2], 0.6, tolerance = 0.12)
})

test_that("missing data and arbitrary category codes are handled", {
  set.seed(15)
  dat <- as.matrix(sim_irt(400, 6, "grm", ncat = 3)) + 1
  dat[sample(length(dat), 200)] <- NA
  dat <- as.data.frame(dat)
  fit <- irt_fit(dat, "grm")
  expect_s3_class(fit, "irtfit")
  expect_equal(fit$labels[[1]], c(1, 2, 3))
  expect_error(irt_fit(dat, "2pl"), "polytomous")
})

test_that("anova compares nested models", {
  set.seed(16)
  dat <- sim_irt(500, 6, "2pl")
  a <- anova(irt_fit(dat, "rasch", se = "none"), irt_fit(dat, "2pl", se = "none"))
  expect_equal(nrow(a), 2)
  expect_true(a$LR[2] > 0)
})
