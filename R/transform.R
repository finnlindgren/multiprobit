#' @title Internal utilities
#' @details \code{tri_solve} solves triangular systems with back/forwardsolve
#' @keywords internal
#' @rdname internal_utils

tri_solve <- function(A, b, lower_tri = FALSE) {
  if (missing(b)) {
    I <- diag(1.0, nrow(A))
    if (!lower_tri) {
      backsolve(A, I)
    } else {
      forwardsolve(A, I)
    }
  } else {
    if (!lower_tri) {
      backsolve(A, b)
    } else {
      forwardsolve(A, b)
    }
  }
}


#' @details \code{qchisq_pnorm} evaluates \code{qchisq(pnorm(...))} with
#' attempt at numerical stability.
#' @keywords internal
#' @importFrom stats qchisq pnorm
#' @rdname internal_utils

qchisq_pnorm <- function(x, df) {
  qchisq(pnorm(x, log.p = TRUE),
    df = df, log.p = TRUE
  )
}

#' @keywords internal
#' @importFrom stats pchisq qnorm
#' @rdname internal_utils

qnorm_pchisq <- function(x, df) {
  qnorm(pchisq(x, df = df, log.p = TRUE), log.p = TRUE)
}

#' @keywords internal
#' @importFrom stats qbeta pchisq
#' @rdname internal_utils

qbeta_pchisq <- function(x,
                         df,
                         shape1, shape2,
                         ncp) {
  qbeta(pchisq(q = x, df = df, ncp = ncp, log.p = TRUE),
    shape1 = shape1, shape2 = shape2, ncp = ncp, log.p = TRUE
  )
}

#' @keywords internal
#' @importFrom stats pbeta qchisq
#' @rdname internal_utils

qchisq_pbeta <- function(x,
                         shape1, shape2,
                         df,
                         ncp) {
  qchisq(
    pbeta(x, shape1 = shape1, shape2 = shape2, ncp = ncp, log.p = TRUE),
    df = df, ncp = ncp, log.p = TRUE
  )
}


# Wishart ####

#' @title Transform latent variables to Wishart
#' @description Transform latent iid \eqn{N(0,1)} to Wishart matrices
#' @param x A vector of latent variables identifying the Wishart matrix,
#' length \eqn{d(d+1)/2}
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{W_chol}{The Cholesky factor of the Wishart matrix}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When \code{lower_chol} is \code{FALSE}, \code{W_chol} is
#' \code{B_chol \%*\% V_chol}, otherwise \code{V_chol \%*\% B_chol}.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname latent_to_wishart

latent_to_wishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  R <- diag(sqrt(qchisq_pnorm(x[seq_len(d)], df = df - seq_len(d) + 1)), d, d)
  R[upper.tri(R)] <- x[d + seq_len(d * (d - 1) / 2)]
  if (lower_chol) {
    W_chol <- V_chol %*% t(R)
    list(W_chol = W_chol, B_chol = t(R))
  } else {
    W_chol <- R %*% V_chol
    list(W_chol = W_chol, B_chol = R)
  }
}

#' @title Transform Wishart to latent variables
#' @description Transforms a Wishart matrix to latent iid \eqn{N(0,1)} variables
#' @param W_chol The Cholesky factor of the Wishart matrix
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
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname latent_to_wishart

latent_from_wishart <- function(W_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  if (!lower_chol) {
    # Convert input to lower Cholesky format
    W_chol <- t(W_chol)
    V_chol <- t(V_chol)
  }
  B_chol <- tri_solve(V_chol, W_chol, lower_tri = TRUE)
  x <- numeric(d * (d + 1) / 2)
  x[seq_len(d)] <- qnorm_pchisq(diag(B_chol)^2, df = df - seq_len(d) + 1)
  x[d + seq_len(d * (d - 1) / 2)] <- t(B_chol)[upper.tri(B_chol)]
  if (lower_chol) {
    list(x = x, B_chol = B_chol)
  } else {
    list(x = x, B_chol = t(B_chol))
  }
}

# Normalised Wishart ####

#' @title Transform latent variables to Normalised Wishart
#' @description Transform latent iid \eqn{N(0,1)} to Normalised Wishart
#'   matrices
#' @param x A vector of latent variables identifying the Normalised Wishart
#' matrix, length \eqn{d(d-1)/2}
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{W_chol}{The Cholesky factor of the Normalised Wishart matrix}
#'   \item{s}{Vector of scaling values}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When \code{lower_chol} is \code{FALSE},
#' \code{W_chol = (B_chol \%*\% V_chol) \%*\% diag(s)}.
#'
#' When \code{lower_chol} is \code{TRUE},
#' \code{W_chol = s * (V_chol \%*\% B_chol)}.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname latent_to_nwishart

latent_to_nwishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  if (!lower_chol) {
    V_chol <- t(V_chol)
  }
  s_vec <- c(1 / V_chol[1, 1], numeric(d - 1))
  B_chol <- W_chol <- matrix(0, d, d)
  B_chol[1, 1] <- W_chol[1, 1] <- 1
  prev_index <- cumsum(seq_len(d) - 1)
  for (k in seq_len(d - 1)) {
    B_chol[k + 1, seq_len(k)] <- u <- x[prev_index[k] + seq_len(k)]
    lambda_vec <- (V_chol[k + 1, seq_len(k)] / V_chol[k + 1, k + 1]) %*%
      B_chol[seq_len(k), seq_len(k), drop = FALSE]
    uu <- lambda_vec + u
    uu_sq <- sum(uu^2)
    uu_norm <- sqrt(uu_sq)
    lambda <- sum(lambda_vec^2)
    if (uu_norm < .Machine$double.eps^0.5) {
      # Close to the origin; handle specially
      uu_norm <- .Machine$double.eps^0.5
      uu_sq <- .Machine$double.eps
    }
    z <-
      qbeta_pchisq(
        uu_sq,
        df = k,
        shape1 = k / 2, shape2 = (df - k) / 2,
        ncp = lambda
      )
    tmp <- sqrt(z) / uu_norm
    s_vec[k + 1] <- tmp / V_chol[k + 1, k + 1]
    B_chol[k + 1, k + 1] <- sqrt(1 - z) / tmp
    W_chol[k + 1, seq_len(k)] <- tmp * uu
    W_chol[k + 1, k + 1] <- sqrt(1 - z)
  }
  if (lower_chol) {
    list(W_chol = W_chol, s = s_vec, B_chol = B_chol)
  } else {
    list(W_chol = t(W_chol), s = s_vec, B_chol = t(B_chol))
  }
}

#' @title Transform Normalised Wishart to latent variables
#' @description Transform a Normalised Wishart matrix to latent iid
#' \eqn{N(0,1)} variables
#' @param W_chol The Cholesky factor of the Normalised Wishart matrix
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{x}{A vector of latent variables identifying the Normalised Wishart
#'   matrix, length \eqn{d(d-1)/2}}
#'   \item{s}{Vector of scaling values}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When \code{lower_chol} is \code{FALSE},
#' \code{W_chol = (B_chol \%*\% V_chol) \%*\% diag(s)}.
#'
#' When \code{lower_chol} is \code{TRUE},
#' \code{W_chol = s * (V_chol \%*\% B_chol)}.
#'
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname latent_to_nwishart

latent_from_nwishart <- function(W_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(W_chol)
  if (!lower_chol) {
    W_chol <- t(W_chol)
    V_chol <- t(V_chol)
  }
  # Make sure the rows are normalised:
  W_chol <- (1 / rowSums(W_chol^2)^0.5) * W_chol
  s_vec <- c(1 / V_chol[1, 1], numeric(d - 1))
  x <- numeric(d * (d - 1) / 2)
  B_chol <- matrix(0, d, d)
  B_chol[1, 1] <- 1
  prev_index <- cumsum(seq_len(d) - 1)
  for (k in seq_len(d - 1)) {
    z <- 1 - W_chol[k + 1, k + 1]^2
    lambda_vec <-
      (V_chol[k + 1, seq_len(k), drop = FALSE] / V_chol[k + 1, k + 1]) %*%
      B_chol[seq_len(k), seq_len(k), drop = FALSE]
    lambda <- sum(lambda_vec^2)
    uu_norm <-
      sqrt(qchisq_pbeta(z,
        shape1 = k / 2, shape2 = (df - k) / 2,
        df = k,
        ncp = lambda
      ))
    u <- ((uu_norm / sqrt(z)) *
      W_chol[k + 1, seq_len(k), drop = FALSE] - lambda_vec)

    # Store the latent variables
    x[prev_index[k] + seq_len(k)] <- u
    B_chol[k + 1, seq_len(k)] <- u

    tmp <- sqrt(z) / uu_norm
    s_vec[k + 1] <- tmp / V_chol[k + 1, k + 1]
    B_chol[k + 1, k + 1] <- W_chol[k + 1, k + 1] / tmp
  }
  if (lower_chol) {
    list(x = x, s = s_vec, B_chol = B_chol)
  } else {
    list(x = x, s = s_vec, B_chol = t(B_chol))
  }
}






# Inverse Wishart ####

#' @title Transform latent variables to Inverse Wishart
#' @description Transform latent iid \eqn{N(0,1)} to Inverse Wishart matrices
#' @param x A vector of latent variables identifying the Inverse Wishart matrix,
#' length \eqn{d(d+1)/2}
#' @param V_chol Transposed inverse of the Cholesky factor of the Inverse
#' Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{W_chol}{The Cholesky factor of the Inverse Wishart matrix}
#'   \item{B_chol}{The inner Cholesky factor of the inverse Barlett
#'     decomposition}
#' When \code{lower_chol} is \code{FALSE}, \code{W_chol} is
#' \code{B_chol \%*\% V_chol}, otherwise \code{V_chol \%*\% B_chol}.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname latent_to_iwishart

latent_to_iwishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  rd <- rev(seq_len(d))
  # Compute Wishart matrix, in reverse ordering
  W <- latent_to_wishart(
    x = x,
    V_chol = t(tri_solve(V_chol, lower_tri = lower_chol))[rd, rd],
    df = df,
    lower_chol = lower_chol
  )
  # Invert, and go back to original ordering
  W$W_chol <- t(tri_solve(W$W_chol, lower_tri = lower_chol))[rd, rd]
  W$B_chol <- t(tri_solve(W$B_chol, lower_tri = lower_chol))[rd, rd]
  W
}




#' @title Transform Inverse Wishart to latent variables
#' @description Transforms an Inverse Wishart matrix to latent iid
#' \eqn{N(0,1)} variables
#' @param W_chol The Cholesky factor of the Inverse Wishart matrix
#' @return A list with components
#'   \item{x}{A vector of latent variables identifying the Wishart matrix,
#'   length \eqn{d(d+1)/2}}
#'   \item{B_chol}{The inner Cholesky factor of the inverse Barlett
#'     decomposition}
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname latent_to_iwishart

latent_from_iwishart <- function(W_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  rd <- rev(seq_len(d))
  # Reverse order and invert:
  # Pre-normalise V_chol for numerical stability; renormalisation afterwards
  # means we don't need to store this rescaling
  W <- latent_from_wishart(
    t(tri_solve(W_chol, lower_tri = lower_chol))[rd, rd],
    t(tri_solve(V_chol, lower_tri = lower_chol))[rd, rd],
    df = df,
    lower_chol = lower_chol
  )
  list(
    x = W$x,
    B_chol = t(tri_solve(W$B_chol, lower_tri = lower_chol))[rd, rd]
  )
}

# Normalised Inverse Wishart ####

normalise_ti_V_chol_rd <- function(V_chol, rd, lower_chol) {
  ti_V_chol_rd <- t(tri_solve(V_chol, lower_tri = lower_chol))[rd, rd]
  if (lower_chol) {
    ti_V_chol_rd <- (1 / rowSums(ti_V_chol_rd^2)^0.5) * ti_V_chol_rd
  } else {
    ti_V_chol_rd <- ti_V_chol_rd %*% diag(1 / colSums(ti_V_chol_rd^2)^0.5)
  }
  ti_V_chol_rd
}

#' @title Transform latent variables to Normalised Inverse Wishart
#' @description Transform latent iid \eqn{N(0,1)} to Normalised Inverse Wishart
#'   matrices
#' @param x A vector of latent variables identifying the Normalised Inverse
#' Wishart matrix, length \eqn{d(d-1)/2}
#' @param V_chol The Cholesky factor of the Inverse Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol \code{TRUE} if lower triangular Cholesky factors are used
#'   (default = \code{FALSE})
#' @return A list with components
#'   \item{W_chol}{The Cholesky factor of the Normalised Inverse Wishart matrix}
#'   \item{s}{Vector of scaling values}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When \code{lower_chol} is \code{FALSE},
#' \code{W_chol = (B_chol \%*\% V_chol) \%*\% diag(s)}.
#'
#' When \code{lower_chol} is \code{TRUE},
#' \code{W_chol = s * (V_chol \%*\% B_chol)}.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname latent_to_niwishart

latent_to_niwishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  rd <- rev(seq_len(d))
  # Compute NWishart matrix, in reverse ordering
  W <- latent_to_nwishart(
    x = x,
    V_chol = normalise_ti_V_chol_rd(
      V_chol = V_chol,
      rd = rd,
      lower_chol = lower_chol
    ),
    df = df,
    lower_chol = lower_chol
  )
  # Invert, and go back to original ordering
  W$B_chol <- t(tri_solve(W$B_chol, lower_tri = lower_chol))[rd, rd]
  # Renormalise
  if (lower_chol) {
    # Need to construct W_chol here, since V_chol was renormalised
    W$W_chol <- V_chol %*% W$B_chol
    scale <- 1 / rowSums(W$W_chol^2)^0.5
    W$W_chol <- scale * W$W_chol
  } else {
    # Need to construct W_chol here, since V_chol was renormalised
    W$W_chol <- W$B_chol %*% V_chol
    scale <- 1 / colSums(W$W_chol^2)^0.5
    W$W_chol <- W$W_chol %*% diag(scale)
  }
  W$s <- scale
  W
}

#' @title Transform Normalised Inverse Wishart to latent variables
#' @description Transform a Normalised Inverse Wishart matrix to latent iid
#' \eqn{N(0,1)} variables
#' @param W_chol The Cholesky factor of the Normalised Inverse Wishart matrix
#' @return A list with components
#'   \item{x}{A vector of latent variables identifying the Normalised Inverse
#'   Wishart matrix, length \eqn{d(d-1)/2}}
#'   \item{s}{Vector of scaling values}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#'
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname latent_to_niwishart

latent_from_niwishart <- function(W_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(W_chol)
  rd <- rev(seq_len(d))
  # Reverse order and invert
  W <- latent_from_nwishart(
    t(tri_solve(W_chol, lower_tri = lower_chol))[rd, rd],
    V_chol = normalise_ti_V_chol_rd(
      V_chol = V_chol,
      rd = rd,
      lower_chol = lower_chol
    ),
    df = df,
    lower_chol = lower_chol
  )
  W$B_chol <- t(tri_solve(W$B_chol, lower_tri = lower_chol))[rd, rd]
  if (lower_chol) {
    W$s <- 1 / rowSums((V_chol %*% W$B_chol)^2)^0.5
  } else {
    W$s <- 1 / colSums((W$B_chol %*% V_chol)^2)^0.5
  }
  W
}



# Wishart density ####

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param W PARAM_DESCRIPTION, Default: NULL
#' @param x PARAM_DESCRIPTION, Default: NULL
#' @param W_chol PARAM_DESCRIPTION, Default: NULL
#' @param V_chol PARAM_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param lower_chol PARAM_DESCRIPTION, Default: FALSE
#' @param log PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname dwishart

dwishart <- function(W = NULL, x = NULL, W_chol = NULL, V_chol, df,
                     lower_chol = FALSE, log = FALSE) {
  if (!requireNamespace("CholWishart")) {
    stop("Missing package 'CholWishart'.")
  }
  if (is.null(W) && is.null(x) && is.null(W_chol)) {
    stop("One of W, x and W_chol must be non-NULL.")
  }
  if ((!is.null(W)) + (!is.null(x)) + (!is.null(W_chol)) > 1) {
    stop("Only one of W, x and W_chol may be non-NULL.")
  }
  if (missing(V_chol)) {
    stop("V_chol must be supplied")
  }
  if (missing(df)) {
    stop("df must be supplied")
  }
  if (!is.null(W)) {
    W_chol <- chol(W)
    if (lower_chol) {
      W_chol <- t(W_chol)
    }
  } else if (!is.null(x)) {
    W_chol <- latent_to_wishart(x,
      V_chol = V_chol, df = df,
      lower_chol = lower_chol
    )$W_chol
  }
  d <- nrow(W_chol)
  if (lower_chol) {
    LVi_LW <- tri_solve(V_chol, W_chol, lower_tri = TRUE)
  } else {
    LVi_LW <- tri_solve(t(V_chol), t(W_chol), lower_tri = TRUE)
  }
  log_p <-
    as.vector(
      -df * d / 2 * log(2) - CholWishart::lmvgamma(df / 2, d) -
        df * sum(log(diag(V_chol))) +
        (df - d - 1) * sum(log(diag(W_chol))) - 0.5 * sum(LVi_LW * LVi_LW)
    )
  if (log) {
    log_p
  } else {
    exp(log_p)
  }
}
