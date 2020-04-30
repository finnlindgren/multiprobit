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

  for (hessian in c("block", "full")) {
    for (direct_X in c(TRUE, FALSE)) {
      if (direct_X) {
        model <- mp_model(response = Y, X = X, df = df, prec_beta = prec_beta)
      } else {
        model <- mp_model(
          response = Y, formula = ~ -1 + x1,
          data = data.frame(x1 = X[, 1]),
          df = df, prec_beta = prec_beta
        )
      }

      expect_equal(inherits(model, "mp_model"), TRUE)
      s <- summary(model)
      expect_equal(inherits(s, "mp_model_summary"), TRUE)
      expect_equal(names(s), c("beta", "u", "Sigma"))

      options <-
        mp_options(
          max_iter = 1,
          gaussint = list(max.threads = 1),
          hessian = hessian
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

        expect_equal(inherits(opt[[strategy]], "mp_estimate"), TRUE)
        expect_equal(class(opt[[strategy]]), "mp_estimate")
        expect_equal(class(opt[[strategy]][["result"]])[1], "data.frame")
        expect_equal(class(opt[[strategy]][["counts"]])[1], "data.frame")

        s <- summary(opt[[strategy]])
        expect_equal(inherits(s, "mp_estimate_summary"), TRUE)
        expect_equal(names(s), c("beta", "u", "Sigma"))
      }
    }
  }
})
