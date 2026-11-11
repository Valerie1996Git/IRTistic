test_that("linking recovers a known transformation", {
  set.seed(51)
  base <- sim_irt(10, 20, "2pl"); tr <- attr(base, "truth")
  dy <- sim_irt(1500, a = tr$a, b = tr$b, model = "2pl", theta = matrix(rnorm(1500)))
  dx <- sim_irt(1500, a = tr$a, b = tr$b, model = "2pl", theta = matrix(rnorm(1500, 0.5, 1.2)))
  fy <- irt_fit(dy, "2pl", se = "none"); fx <- irt_fit(dx, "2pl", se = "none")
  for (m in c("stocking-lord", "haebara", "mean-mean", "mean-sigma")) {
    lk <- irt_link(fx, fy, method = m)
    expect_equal(lk$A, 1.2, tolerance = 0.1, info = m)
    expect_equal(lk$B, 0.5, tolerance = 0.15, info = m)
  }
})

test_that("Rasch linking fixes the slope", {
  set.seed(52)
  b <- make_bank(b = rnorm(10), model = "rasch")
  lk <- irt_link(.bank_transform(b, 1, 0.7), b)
  expect_equal(lk$A, 1)
  expect_equal(lk$B, -0.7, tolerance = 1e-4)
})

test_that("score tables and true-score equating behave", {
  bank <- make_bank(b = seq(-2, 2, length.out = 10), a = 1.3)
  st <- score_table(bank)
  expect_equal(nrow(st), 11)
  expect_equal(sum(st$expected_prop), 1, tolerance = 1e-8)
  expect_true(all(diff(st$theta) > 0))
  eq <- equate_true_score(bank, bank)
  expect_equal(eq$score_y, eq$score_x, tolerance = 1e-6)
})
