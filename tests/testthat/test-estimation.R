test_that("multiprobit", {
  N <- 6
  d <- 2
  J <- 2

  set.seed(1L)
  X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
  B <- matrix(0.5, J, d)
  Y <- matrix(rnorm(N * d, mean = as.vector(X %*% B)) > 0, N, d)
  df <- d + 1
  prec_beta <- 0.1

  model <- mp_model(response = Y, X = X, df = df, prec_beta = prec_beta)
  options <-
    mp_options(
      max_iter = 1,
      gaussint = list(max.threads = 1)
    )
  expect_true(mp_options_check(mp_options(
    mp_options_default(),
    options
  )))

  timing <- list()
  opt <- list()
  for (strategy in c("alternating", "joint", "stepwise")) {
    options <- mp_options(options, strategy = strategy)
    expect_true(mp_options_check(mp_options(
      mp_options_default(),
      options
    )))
    opt[[strategy]] <- multiprobit(
      model = model,
      options = options
    )
    expect_equal(class(opt[[strategy]]), "mp_estimate")
    expect_equal(class(opt[[strategy]][["result"]])[1], "data.frame")
    expect_equal(class(opt[[strategy]][["counts"]])[1], "data.frame")
  }
})
