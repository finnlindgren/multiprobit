# latent/wishart transformation self consistent ####
test_that("self consistent latent/wishart transformation", {
  LV <- matrix(c(1,2,3, 0,4,5, 0,0,6), 3, 3)
  LV <- t(chol(LV %*% t(LV)))
  df <- nrow(LV) + 0.5
  x <- 1:6

  x_to_LC <- latent_to_wishart(x = x, LV = LV, df = df)

  LC_check <- LV %*% x_to_LC$L
  expect_equal(x_to_LC$LC, LC_check)

  LC_to_x <- latent_from_wishart(LC = x_to_LC$LC, LV = LV, df = df)

  expect_equal(LC_to_x$x, x)

  L_check <- solve(LV, x_to_LC$LC)
  expect_equal(LC_to_x$L, L_check)

  expect_equal(x_to_LC$L, LC_to_x$L)
})


# latent/nwishart transformation self consistent ####
test_that("self consistent latent/nwishart transformation", {
  LV <- matrix(c(1,2,3, 0,4,5, 0,0,6), 3, 3)
  LV <- t(chol(LV %*% t(LV)))
  df <- nrow(LV) + 0.5
  x <- 1:3

  x_to_LC <- latent_to_nwishart(x = x, LV = LV, df = df)

  LC_check <- x_to_LC$s * (LV %*% x_to_LC$L)
  expect_equal(x_to_LC$LC, LC_check)

  LC_to_x <- latent_from_nwishart(LC = x_to_LC$LC, LV = LV, df = df)

  expect_equal(LC_to_x$x, x)

  L_check <- solve(LV, x_to_LC$LC / LC_to_x$s)
  expect_equal(LC_to_x$L, L_check)

  expect_equal(x_to_LC$s, LC_to_x$s)
  expect_equal(x_to_LC$L, LC_to_x$L)
})
