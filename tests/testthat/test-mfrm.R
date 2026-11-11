set.seed(31)
r <- sim_mfrm(120, 5, 4, ncat = 4, design = 0.8)
fit <- mfrm_fit(r, "person", c("rater", "criterion"), "score")

test_that("MFRM recovers facet measures", {
  tr <- attr(r, "truth")
  cf <- coef(fit)
  expect_gt(cor(cf$estimate[cf$facet == "rater"], tr$rater), 0.95)
  expect_gt(cor(cf$estimate[cf$facet == "criterion"], tr$criterion), 0.95)
  expect_true(all(is.finite(cf$se)))
})

test_that("one-facet MFRM equals the IRT rating scale model", {
  r1 <- r[r$rater == "R01", c("person", "criterion", "score")]
  f1 <- mfrm_fit(r1, "person", "criterion", "score")
  wide <- stats::reshape(r1, idvar = "person", timevar = "criterion", direction = "wide")
  f2 <- irt_fit(wide[, -1], "rsm", se = "none")
  expect_equal(f1$loglik, f2$loglik, tolerance = 1e-3)
})

test_that("reports, bias and unexpected responses work", {
  rep <- facet_report(fit)
  expect_named(rep$measures, c("person", "rater", "criterion"))
  expect_true(all(c("fair_avg", "infit", "outfit_z") %in% names(rep$measures$rater)))
  expect_equal(nrow(rep$categories), 4)
  b <- mfrm_bias(fit, c("rater", "criterion"))
  expect_equal(nrow(b), 20)
  expect_s3_class(unexpected(fit, 2), "data.frame")
  pcm <- mfrm_fit(r, "person", c("rater", "criterion"), "score", model = "pcm", step_facet = "criterion")
  expect_equal(nrow(facet_report(pcm)$categories), 16)
})
