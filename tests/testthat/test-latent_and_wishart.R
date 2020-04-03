# latent/wishart transformation self consistent ####
test_that("self consistent latent/wishart transformation", {
  V_chol <- matrix(c(1,2,3, 0,4,5, 0,0,6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:6

  x_to_S <- latent_to_wishart(x = x, V_chol = V_chol, df = df, lower_chol = TRUE)

  S_check <- V_chol %*% x_to_S$B_chol
  expect_equal(x_to_S$S_chol, S_check)

  S_to_x <- latent_from_wishart(S_chol = x_to_S$S_chol, V_chol = V_chol, df = df,
                                lower_chol = TRUE)

  expect_equal(S_to_x$x, x)
  expect_equal(x_to_S$B_chol, S_to_x$B_chol)

  B_check <- solve(V_chol, x_to_S$S_chol)
  expect_equal(S_to_x$B_chol, B_check)

  x_to_S_u <- latent_to_wishart(x = x,
                                V_chol = t(V_chol),
                                df = df,
                                lower_chol = FALSE)
  S_to_x_u <- latent_from_wishart(S_chol = t(x_to_S$S_chol),
                                  V_chol = t(V_chol),
                                  df = df,
                                  lower_chol = FALSE)

  expect_equal(x_to_S_u$S_chol, t(x_to_S$S_chol))
  expect_equal(x_to_S_u$B_chol, t(x_to_S$B_chol))
  expect_equal(S_to_x_u$x, S_to_x$x)
  expect_equal(S_to_x_u$B_chol, t(S_to_x$B_chol))
})


# latent/nwishart transformation self consistent ####
test_that("self consistent latent/nwishart transformation", {
  V_chol <- matrix(c(1,2,3, 0,4,5, 0,0,6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:3

  x_to_C <- latent_to_nwishart(x = x, V_chol = V_chol, df = df,
                               lower_chol = TRUE)

  C_check <- x_to_C$s * (V_chol %*% x_to_C$B_chol)
  expect_equal(x_to_C$C_chol, C_check)

  C_to_x <- latent_from_nwishart(C_chol = x_to_C$C_chol,
                                 V_chol = V_chol,
                                 df = df,
                                 lower_chol = TRUE)

  expect_equal(C_to_x$x, x)
  expect_equal(x_to_C$s, C_to_x$s)
  expect_equal(x_to_C$B_chol, C_to_x$B_chol)

  B_check <- solve(V_chol, x_to_C$C_chol / C_to_x$s)
  expect_equal(C_to_x$B_chol, B_check)

  x_to_C_u <- latent_to_nwishart(x = x,
                                 V_chol = t(V_chol),
                                 df = df,
                                 lower_chol = FALSE)
  C_to_x_u <- latent_from_nwishart(C_chol = t(x_to_C$C_chol),
                                   V_chol = t(V_chol),
                                   df = df,
                                   lower_chol = FALSE)

  expect_equal(x_to_C_u$C_chol, t(x_to_C$C_chol))
  expect_equal(C_to_x_u$s, C_to_x$s)
  expect_equal(x_to_C_u$B_chol, t(x_to_C$B_chol))
  expect_equal(C_to_x_u$x, C_to_x$x)
  expect_equal(C_to_x_u$s, C_to_x$s)
  expect_equal(C_to_x_u$B_chol, t(C_to_x$B_chol))
})
