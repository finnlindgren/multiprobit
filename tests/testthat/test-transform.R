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




# self consistent wm_model transformations ####
test_that("self consistent latent/wishart transformation", {
  V_chol <- matrix(c(1, 2, 3, 0, 4, 5, 0, 0, 6), 3, 3)
  V <- V_chol %*% t(V_chol)
  df <- nrow(V_chol) + 0.5

  for (type in c("wishart", "nwishart", "iwishart", "niwishart")) {
    for (lower_chol in c(FALSE, TRUE)) {
      for (use_chol in c(FALSE, TRUE)) {
        if (use_chol) {
          V_chol <- chol(V)
          if (lower_chol) {
            V_chol <- t(V_chol)
          }
          model <- wm_model(
            type,
            V_chol = V_chol,
            df = df,
            lower_chol = lower_chol
          )
        } else {
          model <- wm_model(type, V = V, df = df, lower_chol = lower_chol)
        }
        x <- seq_len(model$N_latent) / model$N_latent

        W_from_x <- wm_matrix(model, latent = x)
        x_from_W <- wm_latent(model, W = W_from_x)
        expect_equal(x_from_W, x)

        W_chol_from_x <- wm_chol(model, latent = x, lower_chol = lower_chol)
        W_chol_from_x_null <- wm_chol(model, latent = x)
        expect_equal(W_chol_from_x_null, W_chol_from_x)

        W_chol_from_x_t <- wm_chol(model, latent = x, lower_chol = !lower_chol)
        expect_equal(t(W_chol_from_x_t), W_chol_from_x)

        x_from_W_chol <- wm_latent(model,
          W_chol = W_chol_from_x,
          lower_chol = lower_chol
        )
        x_from_W_chol_t <- wm_latent(model,
          W_chol = t(W_chol_from_x),
          lower_chol = !lower_chol
        )
        expect_equal(x_from_W_chol, x)
        expect_equal(x_from_W_chol_t, x)
      }
    }
  }
})






# Wishart density
test_that("Wishart density", {
  if (requireNamespace("CholWishart", quietly = TRUE)) {
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



wm_moment_testing <- function() {
  model <- wm_model("nwishart", V = multiprobit:::sparse_identity(4), df = 5)
  mean_latent <- seq_len(model$N_latent) / model$N_latent
  cov_latent <- multiprobit:::sparse_identity(model$N_latent) / 4

  Mpoint <- wm_matrix(model, latent = mean_latent)
  Spoint <- matrix(1 / model$df^0.5, model$d, model$d) *
    (1 - diag(model$d)) * cov_latent[1, 1]
  M1 <-
    wm_moments_linear(
      model,
      mean_latent = mean_latent,
      cov_latent = cov_latent,
      order = c(1, 1)
    )$mean
  M2 <-
    wm_moments_linear(
      model,
      mean_latent = mean_latent,
      cov_latent = cov_latent,
      order = c(2, 2)
    )$mean
  S1 <-
    wm_moments_linear(
      model,
      mean_latent = mean_latent,
      cov_latent = cov_latent,
      order = c(1, 1)
    )$sd
  S2 <-
    wm_moments_linear(
      model,
      mean_latent = mean_latent,
      cov_latent = cov_latent,
      order = c(2, 2)
    )$sd
  Mmc <- 0
  Smc <- 0
  N_mc <- 10000
  for (loop in seq_len(N_mc)) {
    mat <- wm_matrix(model,
      latent = rnorm(model$N_latent,
        mean = mean_latent,
        sd = Matrix::diag(cov_latent)^0.5
      )
    )
    Mmc <- Mmc + mat
    Smc <- Smc + (mat - Mpoint * 0)^2
  }
  # sum (mat - E)^2 = sum (mat - M + M - E)^2 = sum (mat - M)^2 + Nmc (M - E)^2
  Mmc <- Mmc / N_mc
  Smc <- sqrt(Smc / N_mc - (Mpoint * 0 - Mmc)^2)
  Mmc_se <- Smc / sqrt(N_mc)
  Smc_se <- sqrt(2 * Smc^4) / sqrt(N_mc)

  rbind(
    mean_Mmc =
      c(
        Mpoint = sum((Mpoint - Mmc)^2),
        M1 = sum((M1 - Mmc)^2),
        M2 = sum((M2 - Mmc)^2),
        Mmc = sum(Mmc_se^2)
      ) / sum(Mmc^2),
    mean_M2 =
      c(
        Mpoint = sum((Mpoint - M2)^2),
        M1 = sum((M1 - M2)^2),
        M2 = sum((M2 - M2)^2),
        Mmc = sum((M2 - Mmc)^2) + sum(Mmc_se^2)
      ) / sum(Mmc^2),
    sd =
      c(
        Spoint = sum((Spoint - Smc)^2),
        S1 = sum((S1 - Smc)^2),
        S2 = sum((S2 - Smc)^2),
        Mmc = sum(Smc_se^2)
      ) / sum(Smc^2)
  )
}
