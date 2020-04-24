# latent/wishart transformation self consistent ####
test_that("self consistent latent/wishart transformation", {
  V_chol <- matrix(c(1, 2, 3, 0, 4, 5, 0, 0, 6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:6

  x_to_W <- latent_to_wishart(
    x = x,
    V_chol = V_chol,
    df = df,
    lower_chol = TRUE
  )

  W_check <- V_chol %*% x_to_W$B_chol
  expect_equal(x_to_W$W_chol, W_check)

  W_to_x <- latent_from_wishart(
    W_chol = x_to_W$W_chol,
    V_chol = V_chol,
    df = df,
    lower_chol = TRUE
  )

  expect_equal(W_to_x$x, x)
  expect_equal(x_to_W$B_chol, W_to_x$B_chol)

  B_check <- solve(V_chol, x_to_W$W_chol)
  expect_equal(W_to_x$B_chol, B_check)

  x_to_W_u <- latent_to_wishart(
    x = x,
    V_chol = t(V_chol),
    df = df,
    lower_chol = FALSE
  )
  W_to_x_u <- latent_from_wishart(
    W_chol = t(x_to_W$W_chol),
    V_chol = t(V_chol),
    df = df,
    lower_chol = FALSE
  )

  expect_equal(x_to_W_u$W_chol, t(x_to_W$W_chol))
  expect_equal(x_to_W_u$B_chol, t(x_to_W$B_chol))
  expect_equal(W_to_x_u$x, W_to_x$x)
  expect_equal(W_to_x_u$B_chol, t(W_to_x$B_chol))
})


# latent/nwishart transformation self consistent ####
test_that("self consistent latent/nwishart transformation", {
  V_chol <- matrix(c(1, 2, 3, 0, 4, 5, 0, 0, 6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:3

  x_to_W <- latent_to_nwishart(
    x = x, V_chol = V_chol, df = df,
    lower_chol = TRUE
  )

  W_check <- x_to_W$s * (V_chol %*% x_to_W$B_chol)
  expect_equal(x_to_W$W_chol, W_check)

  W_to_x <- latent_from_nwishart(
    W_chol = x_to_W$W_chol,
    V_chol = V_chol,
    df = df,
    lower_chol = TRUE
  )

  expect_equal(W_to_x$x, x)
  expect_equal(x_to_W$s, W_to_x$s)
  expect_equal(x_to_W$B_chol, W_to_x$B_chol)

  B_check <- solve(V_chol, x_to_W$W_chol / W_to_x$s)
  expect_equal(W_to_x$B_chol, B_check)

  x_to_W_u <- latent_to_nwishart(
    x = x,
    V_chol = t(V_chol),
    df = df,
    lower_chol = FALSE
  )
  W_to_x_u <- latent_from_nwishart(
    W_chol = t(x_to_W$W_chol),
    V_chol = t(V_chol),
    df = df,
    lower_chol = FALSE
  )

  expect_equal(x_to_W_u$W_chol, t(x_to_W$W_chol))
  expect_equal(W_to_x_u$s, W_to_x$s)
  expect_equal(x_to_W_u$B_chol, t(x_to_W$B_chol))
  expect_equal(W_to_x_u$x, W_to_x$x)
  expect_equal(W_to_x_u$s, W_to_x$s)
  expect_equal(W_to_x_u$B_chol, t(W_to_x$B_chol))
})



# latent/iwishart transformation self consistent ####
test_that("self consistent latent/iwishart transformation", {
  V_chol <- matrix(c(1, 2, 3, 0, 4, 5, 0, 0, 6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:6

  x_to_W <- latent_to_iwishart(
    x = x,
    V_chol = V_chol,
    df = df,
    lower_chol = TRUE
  )

  W_check <- V_chol %*% x_to_W$B_chol
  expect_equal(x_to_W$W_chol, W_check)

  W_to_x <- latent_from_iwishart(
    W_chol = x_to_W$W_chol,
    V_chol = V_chol,
    df = df,
    lower_chol = TRUE
  )

  expect_equal(W_to_x$x, x)
  expect_equal(x_to_W$B_chol, W_to_x$B_chol)

  B_check <- solve(V_chol, x_to_W$W_chol)
  expect_equal(W_to_x$B_chol, B_check)

  x_to_W_u <- latent_to_iwishart(
    x = x,
    V_chol = t(V_chol),
    df = df,
    lower_chol = FALSE
  )
  W_to_x_u <- latent_from_iwishart(
    W_chol = t(x_to_W$W_chol),
    V_chol = t(V_chol),
    df = df,
    lower_chol = FALSE
  )

  expect_equal(x_to_W_u$W_chol, t(x_to_W$W_chol))
  expect_equal(x_to_W_u$B_chol, t(x_to_W$B_chol))
  expect_equal(W_to_x_u$x, W_to_x$x)
  expect_equal(W_to_x_u$B_chol, t(W_to_x$B_chol))
})



# latent/niwishart transformation self consistent ####
test_that("self consistent latent/niwishart transformation", {
  V_chol <- matrix(c(1, 2, 3, 0, 4, 5, 0, 0, 6), 3, 3)
  V_chol <- t(chol(V_chol %*% t(V_chol)))
  df <- nrow(V_chol) + 0.5
  x <- 1:3

  x_to_W <- latent_to_niwishart(
    x = x, V_chol = V_chol, df = df,
    lower_chol = TRUE
  )

  W_check <- x_to_W$s * (V_chol %*% x_to_W$B_chol)
  expect_equal(x_to_W$W_chol, W_check)

  W_to_x <- latent_from_niwishart(
    W_chol = x_to_W$W_chol,
    V_chol = V_chol,
    df = df,
    lower_chol = TRUE
  )

  expect_equal(W_to_x$x, x)
  expect_equal(x_to_W$s, W_to_x$s)
  expect_equal(x_to_W$B_chol, W_to_x$B_chol)

  B_check <- solve(V_chol, x_to_W$W_chol / W_to_x$s)
  expect_equal(W_to_x$B_chol, B_check)

  x_to_W_u <- latent_to_niwishart(
    x = x,
    V_chol = t(V_chol),
    df = df,
    lower_chol = FALSE
  )
  W_to_x_u <- latent_from_niwishart(
    W_chol = t(x_to_W$W_chol),
    V_chol = t(V_chol),
    df = df,
    lower_chol = FALSE
  )

  expect_equal(x_to_W_u$W_chol, t(x_to_W$W_chol))
  expect_equal(W_to_x_u$s, W_to_x$s)
  expect_equal(x_to_W_u$B_chol, t(x_to_W$B_chol))
  expect_equal(W_to_x_u$x, W_to_x$x)
  expect_equal(W_to_x_u$s, W_to_x$s)
  expect_equal(W_to_x_u$B_chol, t(W_to_x$B_chol))
})


# Wishart/ density
test_that("Wishart density", {
  if (requireNamespace("CholWishart")) {
    mydiag <- function(d, x = 1.0) {
      Matrix::sparseMatrix(
        i = seq_len(d),
        j = seq_len(d),
        x = x,
        dims = c(d, d)
      )
    }

    V_chol <- mydiag(3)
    df <- 5
    expect_equal(
      dwishart(
        W = mydiag(3), x = NULL, W_chol = NULL, V_chol, df,
        lower_chol = FALSE, log = FALSE
      ),
      0.0001879003
    )
    expect_equal(
      dwishart(
        W = mydiag(3), x = NULL, W_chol = NULL, V_chol, df,
        lower_chol = FALSE, log = TRUE
      ),
      log(0.0001879003),
      tolerance = .Machine$double.eps^0.5 / 0.0001879003
    )
    expect_equal(
      dwishart(
        W = mydiag(3), x = NULL, W_chol = NULL, V_chol, df,
        lower_chol = TRUE, log = FALSE
      ),
      0.0001879003
    )
    x <- latent_from_wishart(mydiag(3), V_chol, df, lower_chol = TRUE)
    expect_equal(
      dwishart(
        W = NULL, x = x$x, W_chol = NULL, V_chol, df,
        lower_chol = TRUE, log = FALSE
      ),
      0.0001879003
    )

    # Input errors
    # Missing df
    expect_error(
      dwishart(
        W = mydiag(3), x = NULL, W_chol = NULL, V_chol,
        lower_chol = FALSE, log = FALSE
      )
    )
    # Missing V_chol
    expect_error(
      dwishart(
        W = mydiag(3), x = NULL, W_chol = NULL, df = df,
        lower_chol = FALSE, log = FALSE
      )
    )
    # Missing input
    expect_error(
      dwishart(
        W = NULL, x = NULL, W_chol = NULL, V_chol, df,
        lower_chol = FALSE, log = FALSE
      )
    )
    # Overspecified input
    expect_error(
      dwishart(
        W = mydiag(3), x = rep(0, 6), W_chol = NULL, V_chol, df,
        lower_chol = FALSE, log = FALSE
      )
    )
    expect_error(
      dwishart(
        W = mydiag(3), x = NULL, W_chol = mydiag(3), V_chol, df,
        lower_chol = FALSE, log = FALSE
      )
    )
    expect_error(
      dwishart(
        W = NULL, x = rep(0, 6), W_chol = mydiag(3), V_chol, df,
        lower_chol = FALSE, log = FALSE
      )
    )
  }
})

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
