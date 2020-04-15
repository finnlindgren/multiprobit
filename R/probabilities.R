#' @title Multivariate Probit Event Probabilities
#' @description Evaluate one or more probabilities for outcomes of a
#' multivariate probit model, with given location and covariance scale
#' parameteters.
#' @param y A matrix with n-by-d elements, where each row is a multivariate
#'   observation, see Details. A vector is interpreted as a single row matrix.
#' @param mu A matrix with n-by-d elements, where each row is an expectation
#'   vector parameter, see Details. A vector is interpreted as a single row
#'   matrix.
#' @param Sigma_chol The Cholesky factor of the covariance matrix parameter,
#'   Default: NULL
#' @param Q_chol The Cholesky factor of the precision matrix matrix parameter,
#'   Default: NULL
#' @param log logical indicating if the log-probability should be returned,
#'   Default: FALSE
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   Default: FALSE
#' @param seed The random seed for \code{excursions::gaussint}, Default: NULL
#' @param ... Further parameters passed on to \code{excursions::gaussint}
#' @return A list with components
#' \describe{
#' \item{P}{A vector of probabilities}
#' \item{E}{A vector with the estimated approximation error for each
#'   probability}
#' }
#' @details Computes the probability
#'   \deqn{P(y_1 > 0, ..., y_d > 0|\mu,\Sigma)}
#' when \eqn{y} is a \eqn{d}-dimensional indicator vector with elements
#' \eqn{y_i=I(z_i \ge 0)}, and \eqn{z} is a \eqn{d}-dimensional
#' Gaussian vector with distribution \eqn{N(\mu,\Sigma)}. Only the inequality
#' for \eqn{y_i} is used, so alternative data representations such as
#' \eqn{-1/+1} will also work as expected.
#'
#' The \eqn{\Sigma} paramter can either be specified though its Cholesky factor
#' \code{Sigma_chol} or through the Cholesky factor of the precision (inverse
#' of \eqn{\Sigma}) \code{Q_chol}.
#' The logical parameter \code{lower_chol} determines if a lower or upper
#' triangular Cholesky factor was supplied.
#'
#' The \code{seed} parameter is passed on to \code{excursions::gaussint}, and
#' allows consistent approximation error when calculating numerical derivatives.
#' @examples
#' if (interactive()) {
#'   mpp(
#'     y = c(1, 0),
#'     mu = c(1, 2),
#'     Sigma_chol = chol(matrix(c(1, 0.5, 0.5, 1), 2, 2))
#'   )
#' }
#' @export
#' @rdname mpp
#'
#' @import excursions
#' @import Matrix

mpp <- function(y, mu,
                Sigma_chol = NULL,
                Q_chol = NULL,
                log = FALSE,
                lower_chol = FALSE,
                seed = NULL,
                ...) {
  # TODO:: generalise permitted y; vector vs matrix etc
  if (is.vector(y)) {
    y <- matrix(y, 1, length(y))
  } else {
    y <- as.matrix(y)
  }
  if (is.vector(mu)) {
    mu <- matrix(mu, 1, length(mu))
  } else {
    mu <- as.matrix(mu)
  }

  d <- NCOL(y)
  if (d != ncol(mu)) {
    stop("mu must have the same number of columns as y")
  }
  if ((NROW(y) > 1) &&
    (NROW(mu) > 1) &&
    (NROW(y) != NROW(mu))) {
    stop("y and mu must either have equal #rows or a single row")
  }
  n <- max(NROW(y), NROW(mu))

  reverse <- FALSE
  if (is.null(Sigma_chol)) {
    if (is.null(Q_chol)) {
      # Use sparseMatrix instead of Diagonal to workaround bug in excursions:
      Q_chol <- Matrix::sparseMatrix(
        i = seq_len(d),
        j = seq_len(d),
        x = 1, ,
        dims = c(d, d)
      )
    } else if (lower_chol) {
      Q_chol <- t(Q_chol)
    }
  } else {
    if (!is.null(Q_chol)) {
      stop("Q_chol must be NULL if Sigma_chol is provided.")
    }
    Q_chol <- tri_solve(Sigma_chol, lower = lower_chol)
    reverse <- TRUE
    perm <- rev(seq_len(d))
    if (lower_chol) {
      Q_chol <- Q_chol[perm, perm]
    } else {
      Q_chol <- t(Q_chol)[perm, perm]
    }
  }

  prob <- list(P = numeric(n), E = numeric(n))
  for (obs in seq_len(n)) {
    y_ <- y[1 + (obs - 1) %% NROW(y), ]
    mu_ <- mu[1 + (obs - 1) %% NROW(mu), ]
    a <- ifelse(y_ > 0, 0, -Inf)
    b <- ifelse(y_ > 0, Inf, 0)
    if (reverse) {
      prob_ <-
        excursions::gaussint(
          mu = mu_[perm],
          Q.chol = Q_chol,
          a = a[perm],
          b = b[perm],
          seed = seed,
          ...
        )
    } else {
      prob_ <-
        excursions::gaussint(
          mu = mu_,
          Q.chol = Q_chol,
          a = a,
          b = b,
          seed = seed,
          ...
        )
    }
    prob$P[obs] <- prob_$P
    prob$E[obs] <- prob_$E
  }
  if (log) {
    prob$E <- prob$E / prob$P
    prob$P <- log(prob$P)
  }
  prob
}




#' @title Multivariate Probit Event Probability Derivatives
#' @description FUNCTION_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @param mu PARAM_DESCRIPTION
#' @param ... Further parameters passed on to
#'   \code{\link{mpp}}
#' @param h Step size for finite differences, Default: 1e-06 (for gradients)
#'   or 1e-04 (for hessian)
#' @param symmetric For gradients, whether to use symmetric finite differences,
#'   Default: FALSE
#' @return OUTPUT_DESCRIPTION gradient
#' @details DETAILS gradient
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1 gradient
#' }
#' }
#' @export
#' @rdname mpp_derivatives
#' @seealso \code{\link{mpp}}

mpp_gradient_mu <- function(y, mu, ..., seed = 1L,
                            h = 1e-6,
                            symmetric = FALSE) {
  d <- length(mu)
  mu_ <- matrix(mu, d, d, byrow = TRUE)
  H <- diag(x = h, nrow = d, ncol = d)
  if (symmetric) {
    prob <- mpp(
      y = y,
      mu = rbind(mu_ - H, mu_ + H),
      ...,
      seed = seed
    )
    D <- (prob$P[d + seq_len(d)] - prob$P[seq_len(d)]) / (2 * h)
  } else {
    prob <- mpp(
      y = y,
      mu = rbind(mu, mu_ + H),
      ...,
      seed = seed
    )
    D <- (prob$P[1 + seq_len(d)] - prob$P[1]) / h
  }
  D
}


#' @param diagonal Logical; if \code{TRUE}, only the diagonal of the hessian
#'   is evaluated, Default: FALSE
#' @return OUTPUT_DESCRIPTION hessian
#' @details DETAILS hessain
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1 hessian
#' }
#' }
#' @export
#' @rdname mpp_derivatives

mpp_hessian_mu <- function(y, mu, ..., seed = 1L,
                           h = 1e-4,
                           diagonal = FALSE) {
  H <- function(h, d) {
    Matrix::sparseMatrix(
      i = seq_len(d), j = seq_len(d),
      x = h, dims = c(d, d)
    )
  }
  Hi <- function(h, d) {
    Matrix::sparseMatrix(
      i = rep(seq_len(d), times = d) +
        rep((seq_len(d) - 1) * d, each = d),
      j = rep(seq_len(d), times = d),
      x = h / 2, dims = c(d^2, d)
    )
  }
  Hj <- function(h, d) {
    Matrix::sparseMatrix(
      i = rep(seq_len(d), times = d) +
        rep((seq_len(d) - 1) * d, each = d),
      j = rep(seq_len(d), each = d),
      x = h / 2, dims = c(d^2, d)
    )
  }

  d <- length(mu)
  if (diagonal) {
    mu_ <- matrix(mu, d, d, byrow = TRUE)
    prob <- mpp(
      y = y,
      mu = rbind(mu, mu_ + H(h, d), mu_ - H(h, d)),
      ...,
      seed = seed
    )
    D <- (prob$P[1 + seq_len(d)]
    + prob$P[1 + d + seq_len(d)]
      - 2 * prob$P[1]) / h^2
  } else {
    mu_ <- matrix(mu, d^2, d, byrow = TRUE)
    prob <- mpp(
      y = y,
      mu = rbind(
        mu,
        mu_ + Hi(h, d) + Hj(h, d),
        mu_ - Hi(h, d) - Hj(h, d),
        mu_ + Hi(h, d) - Hj(h, d) # ,
        # mu_ - Hi(h, d) + Hj(h, d)
      ),
      ...,
      seed = seed
    )
    p_0 <- prob$P[1]
    p_pp <- matrix(prob$P[1 + seq_len(d^2)], d, d)
    p_mm <- matrix(prob$P[1 + d^2 + seq_len(d^2)], d, d)
    p_pm <- matrix(prob$P[1 + 2 * d^2 + seq_len(d^2)], d, d)
    D <- (p_pp + p_mm - p_pm - t(p_pm)) / h^2
  }
  D
}




#' @param u A vector of latent variables identifying the Normalised Wishart
#' matrix, length \eqn{d(d-1)/2}
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used,
#'   Default: \code{FALSE})
#' @param seed The random seed for \code{excursions::gaussint}, Default: NULL
#' @param log Whether to compute gradient of the log-probability,
#'   Default: \code{FALSE}
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname mpp_derivatives

mpp_gradient_u <- function(y,
                           mu,
                           u,
                           V_chol,
                           df,
                           lower_chol = FALSE,
                           seed = 1L,
                           h = 1e-6,
                           symmetric = FALSE,
                           log = FALSE,
                           ...) {
  dof <- length(u)
  d <- (1 + sqrt(1 + 8 * dof)) / 2
  g <- numeric(dof)
  if (!symmetric) {
    C0 <- latent_to_nwishart(x = u,
                             V_chol = V_chol,
                             df = df,
                             lower_chol = FALSE)
    prob0 <- mpp(y = y, mu = mu, Sigma_chol = C0$W_chol, log = log,
                 seed = seed, ...)
    for (loop in seq_len(dof)) {
      H <- rep(c(0, h, 0), times = c(loop - 1, 1, dof - loop))
      C <- latent_to_nwishart(x = u + H,
                              V_chol = V_chol,
                              df = df,
                              lower_chol = FALSE)
      prob <- mpp(y = y, mu = mu, Sigma_chol = C$W_chol, log = log,
                  seed = seed, ...)
      g[loop] <- (prob$P - prob0$P) / h
    }
  } else {
    for (loop in seq_len(dof)) {
      H <- rep(c(0, h, 0), times = c(loop - 1, 1, dof - loop))
      C_p <- latent_to_nwishart(x = u + H,
                                V_chol = V_chol,
                                df = df,
                                lower_chol = FALSE)
      C_m <- latent_to_nwishart(x = u - H,
                                V_chol = V_chol,
                                df = df,
                                lower_chol = FALSE)
      prob_p <- mpp(y = y, mu = mu, Sigma_chol = C_p$W_chol, log = log,
                    seed = seed, ...)
      prob_m <- mpp(y = y, mu = mu, Sigma_chol = C_m$W_chol, log = log,
                    seed = seed, ...)
      g[loop] <- (prob_p$P - prob_m$P) / (2 * h)
    }
  }
  g
}

