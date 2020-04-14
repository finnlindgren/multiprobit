test_that("excursions::gaussint", {

  independent_quadrant <- function(d) {
    excursions::gaussint(mu = rep(0, d),
                         Q.chol = diag(1, nrow = d),
                         a = rep(0, d),
                         b = rep(Inf, d),
                         seed = 1L)
  }

  for (d in c(1, 2, 10)) {
    expect_equal(independent_quadrant(d)$P, 0.5^d)
  }
})

test_that("mpp (multivariate_probit_probability)", {

  myident <- function(d) {
    Matrix::sparseMatrix(i = seq_len(d),
                         j = seq_len(d),
                         x = 1.0,
                         dims = c(d, d))
  }

  independent_quadrant <- function(d, y, mu, ...) {
    Q_chol <- myident(d)
    mpp(y, mu,
        Q_chol = Q_chol,
        ...,
        seed = 1L)
  }
  independent_quadrant_sigma <- function(d, y, mu, ...) {
    Sigma_chol <- myident(d)
    mpp(y, mu,
        Sigma_chol = Sigma_chol,
        ...,
        seed = 1L)
  }

  for (d in c(1, 2, 10)) {
    # vector & vector
    expect_equal(independent_quadrant(
      d = d,
      y = rep(1, d),
      mu = rep(0, d)
    )$P, 0.5^d)
    # matrix & vector
    expect_equal(independent_quadrant(
      d = d,
      y = matrix(1, 1, d),
      mu = rep(0, d)
    )$P, 0.5^d)
    # vector & matrix
    expect_equal(independent_quadrant(
      d = d,
      y = rep(1, d),
      mu = matrix(0, 1, d)
    )$P, 0.5^d)
    # matrix & matrix
    expect_equal(independent_quadrant(
      d = d,
      y = matrix(1, 1, d),
      mu = matrix(0, 1, d)
    )$P, 0.5^d)
    # Equal dimensions?
    expect_error(independent_quadrant(
      d = d,
      y = matrix(1, 1, d),
      mu = matrix(0, 1, d+1)
    ))
    # Equal rows?
    expect_error(independent_quadrant(
      d = d,
      y = matrix(1, 2, d),
      mu = matrix(0, 3, d)
    ))
    # y,mu,chol matchind dimensions?
    expect_error(independent_quadrant(
      d = d+1,
      y = matrix(1, 2, d),
      mu = matrix(0, 2, d)
    ))
    # Sigma_chol
    expect_equal(independent_quadrant_sigma(
      d = d,
      y = rep(1, d),
      mu = rep(0, d)
    )$P, 0.5^d)
    expect_error(independent_quadrant_sigma(
      d = d,
      y = rep(1, d),
      mu = rep(0, d),
      Q_chol = myident(d)
    ))
    # log
    expect_equal(independent_quadrant(
      d = d,
      y = rep(1, d),
      mu = rep(0, d),
      log = TRUE
    )$P, d * log(0.5))
    # lower_chol
    expect_equal(independent_quadrant(
      d = d,
      y = rep(1, d),
      mu = rep(0, d),
      log = TRUE,
      lower_chol = TRUE
    )$P, d * log(0.5))
    expect_equal(independent_quadrant_sigma(
      d = d,
      y = rep(1, d),
      mu = rep(0, d),
      log = TRUE,
      lower_chol = TRUE
    )$P, d * log(0.5))
  }
})


test_that("mpp (sum of all combinations)", {
  d <- 5
  prob <- mpp(
    y = (-1)^floor((.row(dim = c(2^d, d))-1) / 2^(.col(dim = c(2^d, d))-1)),
    mu = rep(1, d))
  expect_equal(sum(prob$P), 1.0)
})


test_that("mpp_gradient_mu", {

  myident <- function(d) {
    Matrix::sparseMatrix(i = seq_len(d),
                         j = seq_len(d),
                         x = 1.0,
                         dims = c(d, d))
  }

  test_grad <- function(symmetric) {
    d <- 3
    Q_chol <- myident(d)
    y <- cbind(rep(c(0, 1), times = c(15, 5)),
               rep(c(0, 1), times = c(10, 10)),
               rep(c(0, 1), times = c(5, 15)))
    mu <- qnorm(c(0.25, 0.5, 0.75))
    D <- c(0, 0, 0)
    for (i in seq_len(nrow(y))) {
      D <- D + mpp_gradient_mu(
        y = y[i, ],
        mu = mu,
        Q_chol = Q_chol,
        log = TRUE,
        symmetric = symmetric)
    }
    D
  }
  expect_equal(test_grad(FALSE),
               c(0, 0, 0),
               tolerance = .Machine$double.eps^(1/3))
  expect_equal(test_grad(TRUE),
               c(0, 0, 0),
               tolerance = .Machine$double.eps^(1/3))
})



test_that("mpp_hessian_mu", {

  myident <- function(d) {
    Matrix::sparseMatrix(i = seq_len(d),
                         j = seq_len(d),
                         x = 1.0,
                         dims = c(d, d))
  }

  test_hess <- function(diagonal) {
    d <- 3
    Q_chol <- myident(d)
    y <- cbind(rep(c(0, 1), times = c(15, 5)),
               rep(c(0, 1), times = c(10, 10)),
               rep(c(0, 1), times = c(5, 15)))
    mu <- qnorm(c(0.25, 0.5, 0.75))
    if (diagonal) {
      H <- rep(0, 3)
    } else {
      H <- matrix(0, 3, 3)
    }
    for (i in seq_len(nrow(y))) {
      H <- H + mpp_hessian_mu(
        y = y[i, ],
        mu = mu,
        Q_chol = Q_chol,
        log = TRUE,
        diagonal = diagonal)
    }
    H
  }
  H <- test_hess(FALSE)
  expect_equal(H[lower.tri(H, diag = TRUE)],
               c(-10.77145, 1.84297e-5, -2.065015e-4,
                 -12.7324, 1.84297e-5, -10.77145),
               tolerance = .Machine$double.eps^(1/3))
  expect_equal(test_hess(TRUE), diag(H))
})


test_that("mpp_gradient_u", {

  myident <- function(d) {
    Matrix::sparseMatrix(i = seq_len(d),
                         j = seq_len(d),
                         x = 1.0,
                         dims = c(d, d))
  }

  test_grad <- function(symmetric) {
    d <- 3
    V_chol <- myident(d)
    df <- 10
    y <- cbind(rep(c(0, 1), times = c(15, 5)),
               rep(c(0, 1), times = c(10, 10)),
               rep(c(0, 1), times = c(5, 15)))
    mu <- qnorm(c(0.25, 0.5, 0.75))
    u <- rep(0, 3)
    D <- rep(0, 3)
    for (i in seq_len(nrow(y))) {
      D <- D + mpp_gradient_u(
        y = y[i, ],
        mu = mu,
        u = u,
        V_chol = V_chol,
        df = df,
        log = TRUE,
        symmetric = symmetric)
    }
    D
  }
  expect_equal(test_grad(TRUE),
               c( 2.319965, 1.270981, 2.387166),
               tolerance = .Machine$double.eps^(1/3))
  expect_equal(test_grad(FALSE),
               c( 2.319965, 1.270981, 2.387166),
               tolerance = .Machine$double.eps^(1/3))
})
