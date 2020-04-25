
# internal utils ####
test_that("Triangular solves", {
  expect_equal(
    tri_solve(matrix(c(1, 0, 1, 1), 2, 2), lower_tri = FALSE),
    matrix(c(1, 0, -1, 1), 2, 2)
  )
  expect_equal(
    tri_solve(matrix(c(1, 1, 0, 1), 2, 2), lower_tri = TRUE),
    matrix(c(1, -1, 0, 1), 2, 2)
  )
  expect_equal(
    tri_solve(matrix(c(1, 0, 1, 1), 2, 2), c(1, 1), lower_tri = FALSE),
    c(0, 1)
  )
  expect_equal(
    tri_solve(matrix(c(1, 1, 0, 1), 2, 2), c(1, 1), lower_tri = TRUE),
    c(1, 0)
  )
  expect_equal(
    tri_solve(matrix(c(1, 1, 0, 1), 2, 2), c(1, 1), lower_tri = TRUE),
    c(1, 0)
  )
})

test_that("Quantile transformations", {
  df <- 2
  shape1 <- 3
  shape2 <- 4
  ncp <- 5
  expect_equal(
    qchisq_pnorm(0, df),
    qchisq(0.5, df)
  )
  expect_equal(
    qnorm_pchisq(qchisq(0.5, df), df),
    qnorm(0.5)
  )
  expect_equal(
    qbeta_pchisq(
      qchisq(0.5, df, ncp),
      df,
      shape1, shape2,
      ncp
    ),
    qbeta(0.5, shape1, shape2, ncp)
  )
  expect_equal(
    qchisq_pbeta(
      qbeta(
        0.5,
        shape1, shape2,
        ncp
      ),
      shape1, shape2,
      df,
      ncp
    ),
    qchisq(0.5, df, ncp)
  )
})
