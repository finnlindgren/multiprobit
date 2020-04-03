#' @import stats

#' @title Transform latent variables to Wishart
#' @description Transforms latent iid \eqn{N(0,1)} to Wishart matrices
#' @param x A vector of latent variables identifying the Wishart matrix,
#' length \eqn{d(d+1)/2}
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{S_chol}{The Cholesky factor of the Wishart matrix}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When \code{lower_chol} is \code{FALSE}, \code{S_chol} is
#' \code{B_chol \%*\% V_chol}, otherwise \code{V_chol \%*\% B_chol}.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname latent_to_wishart

latent_to_wishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  R <- diag(sqrt(qchisq(pnorm(x[seq_len(d)]), df = df - seq_len(d) + 1)))
  R[upper.tri(R)] <- x[d + seq_len(d * (d-1) / 2)]
  if (lower_chol) {
    S_chol <- V_chol %*% t(R)
    list(S_chol = S_chol, B_chol = t(R))
  } else {
    S_chol <- R %*% V_chol
    list(S_chol = S_chol, B_chol = R)
  }
}

#' @title Transform Wishart to latent variables
#' @description Transforms a Wishart matrix to latent iid \eqn{N(0,1)} variables
#' @param S_chol The Cholesky factor of the Wishart matrix
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{x}{A vector of latent variables identifying the Wishart matrix,
#'   length \eqn{d(d+1)/2}}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname latent_from_wishart

latent_from_wishart <- function(S_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  if (!lower_chol) {
    S_chol <- t(S_chol)
    V_chol <- t(V_chol)
  }
  B_chol <- solve(V_chol, S_chol)
  x <- numeric(d * (d+1) / 2)
  x[seq_len(d)] <- qnorm(pchisq(diag(B_chol)^2, df = df - seq_len(d) + 1))
  x[d + seq_len(d * (d-1) / 2)] <- t(B_chol)[upper.tri(B_chol)]
  if (lower_chol) {
    list(x = x, B_chol = B_chol)
  } else {
    list(x = x, B_chol = t(B_chol))
  }
}

#' @title Transform latent variables to Normalised Wishart
#' @description Transforms latent iid \eqn{N(0,1)} to Normalised Wishart
#'   matrices
#' @param x A vector of latent variables identifying the Normalised Wishart
#' matrix, length \eqn{d(d-1)/2}
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{C_chol}{The Cholesky factor of the Normalised Wishart matrix}
#'   \item{s}{Vector of scaling values}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When \code{lower_chol} is \code{FALSE},
#' \code{C_chol = (B_chol \%*\% V_chol) \%*\% diag(s)}.
#'
#' When \code{lower_chol} is \code{TRUE},
#' \code{C_chol = s * (V_chol \%*\% B_chol)}.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname latent_to_nwishart

latent_to_nwishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  if (!lower_chol) {
    V_chol <- t(V_chol)
  }
  s_vec <- c(1 / V_chol[1, 1], numeric(d - 1))
  B_chol <- C_chol <- matrix(0, d, d)
  B_chol[1, 1] <- C_chol[1, 1] <- 1
  prev_index <- cumsum(seq_len(d) - 1)
  for (k in seq_len(d - 1)) {
    B_chol[k + 1, seq_len(k)] <- u <- x[prev_index[k] + seq_len(k)]
    lambda_vec <- (V_chol[k + 1, seq_len(k)] / V_chol[k + 1, k + 1]) %*%
      B_chol[seq_len(k), seq_len(k), drop = FALSE]
    uu <- lambda_vec + u
    uu_sq <- sum(uu^2)
    uu_norm <- sqrt(uu_sq)
    lambda <- sum(lambda_vec^2)
    z <-
      qbeta(
        pchisq(uu_sq, df = k, ncp = lambda),
        shape1 = k / 2,
        shape2 = (df - k) / 2,
        ncp = lambda
      )
    tmp <- sqrt(z) / uu_norm
    s_vec[k + 1] <- tmp / V_chol[k + 1, k + 1]
    B_chol[k + 1, k + 1] <- sqrt(1 - z) / tmp
    C_chol[k + 1, seq_len(k)] <- tmp * uu
    C_chol[k + 1, k + 1] <- sqrt(1 - z)
  }
  if (lower_chol) {
    list(C_chol = C_chol, s = s_vec, B_chol = B_chol)
  } else {
    list(C_chol = t(C_chol), s = s_vec, B_chol = t(B_chol))
  }
}

#' @title Transform Normalised Wishart to latent variables
#' @description Transform a Normalised Wishart matrix to latent iid
#' \eqn{N(0,1)} variables
#' @param C_chol The lower Cholesky factor of the Normalised Wishart matrix
#' @param V_chol The lower Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{x}{A vector of latent variables identifying the Normalised Wishart
#'   matrix, length \eqn{d(d-1)/2}}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When \code{lower_chol} is \code{FALSE},
#' \code{C_chol = (B_chol \%*\% V_chol) \%*\% diag(s)}.
#'
#' When \code{lower_chol} is \code{TRUE},
#' \code{C_chol = s * (V_chol \%*\% B_chol)}.
#'
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname latent_from_nwishart

latent_from_nwishart <- function(C_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(C_chol)
  if (!lower_chol) {
    C_chol <- t(C_chol)
    V_chol <- t(V_chol)
  }
  # Make sure the rows are normalised:
  C_chol <- rowSums(C_chol^2)^0.5 * C_chol
  s_vec <- c(1 / V_chol[1, 1], numeric(d - 1))
  x <- numeric(d * (d - 1) / 2)
  B_chol <- matrix(0, d, d)
  B_chol[1, 1] <- 1
  prev_index <- cumsum(seq_len(d) - 1)
  for (k in seq_len(d - 1)) {
    z <- 1 - C_chol[k + 1, k + 1]^2
    lambda_vec <-
      (V_chol[k + 1, seq_len(k), drop = FALSE] / V_chol[k + 1, k + 1]) %*%
      B_chol[seq_len(k), seq_len(k), drop = FALSE]
    lambda <- sum(lambda_vec^2)
    uu_norm <-
      sqrt(qchisq(
        pbeta(z, shape1 = k / 2, shape2 = (df - k) / 2, ncp = lambda),
        df = k,
        ncp = lambda
      ))
    u <- (uu_norm / sqrt(z)) * C_chol[k + 1, seq_len(k), drop = FALSE] - lambda_vec

    # Store the latent variables
    x[prev_index[k] + seq_len(k)] <- u
    B_chol[k + 1, seq_len(k)] <- u

    tmp <- sqrt(z) / uu_norm
    s_vec[k + 1] <- tmp / V_chol[k + 1, k + 1]
    B_chol[k + 1, k + 1] <- C_chol[k + 1, k + 1] / tmp
  }
  if (lower_chol) {
    list(x = x, s = s_vec, B_chol = B_chol)
  } else {
    list(x = x, s = s_vec, B_chol = t(B_chol))
  }
}
