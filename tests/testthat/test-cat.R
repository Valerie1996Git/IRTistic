set.seed(71)
bank <- make_bank(b = rnorm(150), a = rlnorm(150, 0.1, 0.3), model = "2pl",
                  content = rep(c("A", "B", "C"), 50))

test_that("live CAT session runs and stops on SE", {
  d <- cat_design(bank, max_items = 30, se_target = 0.4)
  s <- cat_start(d)
  while (!cat_done(s)) s <- cat_answer(s, cat_next(s), rbinom(1, 1, 0.5))
  expect_true(s$se <= 0.4 || length(s$idx) == 30)
  expect_equal(attr(cat_done(s), "reason") %in% c("se_target", "max_items"), TRUE)
  expect_error(cat_answer(s, s$items[1], 1), "already administered")
  expect_true(is.na(cat_next(s)))
})

test_that("all selection rules and estimators recover theta", {
  skip_on_cran()
  th <- rnorm(60)
  for (sel in c("MFI", "MPWI", "KL", "random")) for (est in c("EAP", "WLE")) {
    sim <- cat_simulate(cat_design(bank, select = sel, estimator = est, max_items = 25), theta = th, seed = 1)
    expect_gt(summary(sim)$recovery[["correlation"]], 0.85)
  }
})

test_that("content balancing, exposure control and classification work", {
  d <- cat_design(bank, max_items = 20, exposure = "randomesque",
                  content_target = c(A = 0.5, B = 0.25, C = 0.25))
  sim <- cat_simulate(d, theta = rnorm(40), seed = 2)
  props <- prop.table(table(d$content[match(sim$history$item, names(bank$items))]))
  expect_equal(as.vector(props), c(0.5, 0.25, 0.25), tolerance = 0.01)
  dc <- cat_design(bank, cut = 0, min_items = 5, max_items = 40)
  sc <- cat_simulate(dc, theta = rnorm(60), seed = 3)
  expect_gt(summary(sc)$classification_accuracy, 0.8)
  expect_s3_class(plot(sim, "exposure"), "ggplot")
  expect_s3_class(plot(sc, "trace", person = 2), "ggplot")
  expect_s3_class(plot(sc, "conditional", bins = 4), "ggplot")
})

test_that("post-hoc simulation uses real responses", {
  dat <- sim_irt(150, 30, "grm", ncat = 3)
  fit <- irt_fit(dat, "grm", se = "none")
  ph <- cat_simulate(cat_design(fit, max_items = 12), responses = as.matrix(dat))
  expect_gt(cor(ph$results$theta, irt_scores(fit)$theta), 0.85)
})
