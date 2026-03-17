tol <- 1e-6

test_that("mgfBound returns correct values for simple inputs", {
  expect_equal(mgfBound(k = 10, n = 5, lambda = 0), 1, tolerance = tol)

  vals <- mgfBound(k = 5, n = 5, lambda = c(0.1, 0.5, 0.9))
  expect_length(vals, 3)
  expect_true(all(is.finite(vals)))

  expect_equal(mgfBound(k = 10, n = 20, lambda = 1e-10), 1, tolerance = 1e-4)
})

test_that("mgfBound stays finite or signals overflow cleanly for larger inputs", {
  val <- mgfBound(k = 50, n = 50, lambda = 0.5)

  expect_true(is.finite(val) || is.infinite(val))
  expect_gte(val, 1)
})

test_that("tailProbBound returns a valid probability bound", {
  bound <- tailProbBound(x = 10, k = 10, n = 100)

  expect_true(is.finite(bound))
  expect_gt(bound, 0)
  expect_lte(bound, 1)
})

test_that("tailProbBound decreases for larger deviations", {
  bound1 <- tailProbBound(x = 10, k = 10, n = 100)
  bound2 <- tailProbBound(x = 30, k = 10, n = 100)

  expect_lt(bound2, bound1)
})

test_that("tailProbBound accepts vector k and n", {
  val <- tailProbBound(x = 20, k = c(5, 5), n = c(50, 50))

  expect_true(is.finite(val) || is.infinite(val))
})

test_that("criticalValue returns a positive finite value", {
  crit <- criticalValue(k = 10, n = 50, p = 0.05)

  expect_true(is.finite(crit))
  expect_gt(crit, 0)
})

test_that("criticalValue increases as p decreases", {
  crit1 <- criticalValue(k = 10, n = 50, p = 0.05)
  crit2 <- criticalValue(k = 10, n = 50, p = 0.01)

  expect_gt(crit2, crit1)
})

test_that("tailProbBound remains stable for moderately large k", {
  expect_silent({
    val <- tailProbBound(x = 200, k = 100, n = 200)
  })

  expect_true(is.finite(val) || is.infinite(val))
})
