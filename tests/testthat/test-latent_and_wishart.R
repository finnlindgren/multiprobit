# latent/wishart transformation self consistent ####
test_that("self consistent latent/wishart transformation", {
  V_chol <- matrix(c(1,2,3, 0,4,5, 0,0,6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:6

  x_to_W <- latent_to_wishart(x = x,
                              V_chol = V_chol,
                              df = df,
                              lower_chol = TRUE)

  W_check <- V_chol %*% x_to_W$B_chol
  expect_equal(x_to_W$W_chol, W_check)

  W_to_x <- latent_from_wishart(W_chol = x_to_W$W_chol,
                                V_chol = V_chol,
                                df = df,
                                lower_chol = TRUE)

  expect_equal(W_to_x$x, x)
  expect_equal(x_to_W$B_chol, W_to_x$B_chol)

  B_check <- solve(V_chol, x_to_W$W_chol)
  expect_equal(W_to_x$B_chol, B_check)

  x_to_W_u <- latent_to_wishart(x = x,
                                V_chol = t(V_chol),
                                df = df,
                                lower_chol = FALSE)
  W_to_x_u <- latent_from_wishart(W_chol = t(x_to_W$W_chol),
                                  V_chol = t(V_chol),
                                  df = df,
                                  lower_chol = FALSE)

  expect_equal(x_to_W_u$W_chol, t(x_to_W$W_chol))
  expect_equal(x_to_W_u$B_chol, t(x_to_W$B_chol))
  expect_equal(W_to_x_u$x, W_to_x$x)
  expect_equal(W_to_x_u$B_chol, t(W_to_x$B_chol))
})


# latent/nwishart transformation self consistent ####
test_that("self consistent latent/nwishart transformation", {
  V_chol <- matrix(c(1,2,3, 0,4,5, 0,0,6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:3

  x_to_W <- latent_to_nwishart(x = x, V_chol = V_chol, df = df,
                               lower_chol = TRUE)

  W_check <- x_to_W$s * (V_chol %*% x_to_W$B_chol)
  expect_equal(x_to_W$W_chol, W_check)

  W_to_x <- latent_from_nwishart(W_chol = x_to_W$W_chol,
                                 V_chol = V_chol,
                                 df = df,
                                 lower_chol = TRUE)

  expect_equal(W_to_x$x, x)
  expect_equal(x_to_W$s, W_to_x$s)
  expect_equal(x_to_W$B_chol, W_to_x$B_chol)

  B_check <- solve(V_chol, x_to_W$W_chol / W_to_x$s)
  expect_equal(W_to_x$B_chol, B_check)

  x_to_W_u <- latent_to_nwishart(x = x,
                                 V_chol = t(V_chol),
                                 df = df,
                                 lower_chol = FALSE)
  W_to_x_u <- latent_from_nwishart(W_chol = t(x_to_W$W_chol),
                                   V_chol = t(V_chol),
                                   df = df,
                                   lower_chol = FALSE)

  expect_equal(x_to_W_u$W_chol, t(x_to_W$W_chol))
  expect_equal(W_to_x_u$s, W_to_x$s)
  expect_equal(x_to_W_u$B_chol, t(x_to_W$B_chol))
  expect_equal(W_to_x_u$x, W_to_x$x)
  expect_equal(W_to_x_u$s, W_to_x$s)
  expect_equal(W_to_x_u$B_chol, t(W_to_x$B_chol))
})



# latent/iwishart transformation self consistent ####
test_that("self consistent latent/iwishart transformation", {
  V_chol <- matrix(c(1,2,3, 0,4,5, 0,0,6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:6

  x_to_W <- latent_to_iwishart(x = x,
                               V_chol = V_chol,
                               df = df,
                               lower_chol = TRUE)

  W_check <- V_chol %*% x_to_W$B_chol
  expect_equal(x_to_W$W_chol, W_check)

  W_to_x <- latent_from_iwishart(W_chol = x_to_W$W_chol,
                                 V_chol = V_chol,
                                 df = df,
                                 lower_chol = TRUE)

  expect_equal(W_to_x$x, x)
  expect_equal(x_to_W$B_chol, W_to_x$B_chol)

  B_check <- solve(V_chol, x_to_W$W_chol)
  expect_equal(W_to_x$B_chol, B_check)

  x_to_W_u <- latent_to_iwishart(x = x,
                                 V_chol = t(V_chol),
                                 df = df,
                                 lower_chol = FALSE)
  W_to_x_u <- latent_from_iwishart(W_chol = t(x_to_W$W_chol),
                                   V_chol = t(V_chol),
                                   df = df,
                                   lower_chol = FALSE)

  expect_equal(x_to_W_u$W_chol, t(x_to_W$W_chol))
  expect_equal(x_to_W_u$B_chol, t(x_to_W$B_chol))
  expect_equal(W_to_x_u$x, W_to_x$x)
  expect_equal(W_to_x_u$B_chol, t(W_to_x$B_chol))
})



# latent/niwishart transformation self consistent ####
test_that("self consistent latent/niwishart transformation", {
  V_chol <- matrix(c(1,2,3, 0,4,5, 0,0,6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:3

  x_to_W <- latent_to_niwishart(x = x, V_chol = V_chol, df = df,
                               lower_chol = TRUE)

  W_check <- x_to_W$s * (V_chol %*% x_to_W$B_chol)
  expect_equal(x_to_W$W_chol, W_check)

  W_to_x <- latent_from_niwishart(W_chol = x_to_W$W_chol,
                                 V_chol = V_chol,
                                 df = df,
                                 lower_chol = TRUE)

  expect_equal(W_to_x$x, x)
  expect_equal(x_to_W$s, W_to_x$s)
  expect_equal(x_to_W$B_chol, W_to_x$B_chol)

  B_check <- solve(V_chol, x_to_W$W_chol / W_to_x$s)
  expect_equal(W_to_x$B_chol, B_check)

  x_to_W_u <- latent_to_niwishart(x = x,
                                 V_chol = t(V_chol),
                                 df = df,
                                 lower_chol = FALSE)
  W_to_x_u <- latent_from_niwishart(W_chol = t(x_to_W$W_chol),
                                   V_chol = t(V_chol),
                                   df = df,
                                   lower_chol = FALSE)

  expect_equal(x_to_W_u$W_chol, t(x_to_W$W_chol))
  expect_equal(W_to_x_u$s, W_to_x$s)
  expect_equal(x_to_W_u$B_chol, t(x_to_W$B_chol))
  expect_equal(W_to_x_u$x, W_to_x$x)
  expect_equal(W_to_x_u$s, W_to_x$s)
  expect_equal(W_to_x_u$B_chol, t(W_to_x$B_chol))
})
