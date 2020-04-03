test_that("excursions::gaussint", {

  independent_quadrant <- function(d) {
    excursions::gaussint(mu = rep(0, d),
                         Q.chol = diag(1, nrow = d),
                         a = rep(0, d),
                         b = rep(Inf, d),
                         seed = 1L)
  }

  expect_equal(independent_quadrant(d = 1)$P, 0.5)
  expect_equal(independent_quadrant(d = 2)$P, 0.5^2)
  expect_equal(independent_quadrant(d = 10)$P, 0.5^10)
})
