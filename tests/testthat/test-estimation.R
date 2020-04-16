test_that("multiprobit", {
  N <- 6
  d <- 2
  J <- 1

  set.seed(1L)
  X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
  B <- matrix(0.5, J, d)
  Y <- matrix(rnorm(N * d, mean = as.vector(X %*% B)) > 0, N, d)
  df <- d + 1
  prec_beta <- 0.1

  timing <- list()
  opt <- list()
  for (method in c("alternating", "joint", "stepwise")) {
    opt[[method]] <- multiprobit(response = Y,
                                 X = X,
                                 df = df,
                                 prec_beta = prec_beta,
                                 max_iter = 1,
                                 gaussint_options = list(max.threads = 1),
                                 optim_options = list(),
                                 method = method)
    expect_equal(class(opt[[method]]), "list")
    expect_equal(class(opt[[method]][["result"]])[1], "data.frame")
    expect_equal(class(opt[[method]][["counts"]])[1], "data.frame")
  }
})
