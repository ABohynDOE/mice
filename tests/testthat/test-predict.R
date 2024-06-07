test_that("predict.mira: nhanes lm", {
  data(nhanes)
  imp <- mice::mice(nhanes, maxit = 2, m = 2, seed = 1, print = FALSE, use.matcher = TRUE)
  fit_mira <- with(imp, lm(chl ~ age + bmi))
  tmp <- stats::predict(fit_mira)
  expect_true(inherits(tmp, "data.frame"))
})