

# Wishart ####

#' @title Transform latent variables to Wishart
#' @description Transform latent iid \eqn{N(0,1)} to Wishart matrices
#' @param x A vector of latent variables identifying the Wishart matrix,
#' length \eqn{d(d+1)/2}
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol `TRUE` if lower triangular Cholesky factors are used
#'   (default = `FALSE`)
#' @return A list with components
#'   \item{W_chol}{The Cholesky factor of the Wishart matrix}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When `lower_chol` is `FALSE`, `W_chol` is
#' `B_chol \%*\% V_chol`, otherwise `V_chol \%*\% B_chol`.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @keywords internal
#' @note This is an internal non-exported function
#' @rdname latent_to_wishart

latent_to_wishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  R <- Matrix::diag(sqrt(qchisq_pnorm(
    x[seq_len(d)],
    df = df - seq_len(d) + 1
  )), d, d)
  R[upper.tri(R)] <- x[d + seq_len(d * (d - 1) / 2)]
  if (lower_chol) {
    W_chol <- V_chol %*% Matrix::t(R)
    list(W_chol = W_chol, B_chol = Matrix::t(R))
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
#' @param lower_chol `TRUE` if lower triangular Cholesky factors are used
#'   (default = `FALSE`)
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
#' @keywords internal
#' @note This is an internal non-exported function
#' @rdname latent_to_wishart

latent_from_wishart <- function(W_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  if (!lower_chol) {
    # Convert input to lower Cholesky format
    W_chol <- Matrix::t(W_chol)
    V_chol <- Matrix::t(V_chol)
  }
  B_chol <- tri_solve(V_chol, W_chol, lower_tri = TRUE)
  x <- numeric(d * (d + 1) / 2)
  x[seq_len(d)] <- qnorm_pchisq(Matrix::diag(B_chol)^2,
    df = df - seq_len(d) + 1
  )
  x[d + seq_len(d * (d - 1) / 2)] <- Matrix::t(B_chol)[upper.tri(B_chol)]
  if (lower_chol) {
    list(x = x, B_chol = B_chol)
  } else {
    list(x = x, B_chol = Matrix::t(B_chol))
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
#' @param lower_chol `TRUE` if lower triangular Cholesky factors are used
#'   (default = `FALSE`)
#' @return A list with components
#'   \item{W_chol}{The Cholesky factor of the Normalised Wishart matrix}
#'   \item{s}{Vector of scaling values}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When `lower_chol` is `FALSE`,
#' `W_chol = (B_chol \%*\% V_chol) \%*\% diag(s)`.
#'
#' When `lower_chol` is `TRUE`,
#' `W_chol = s * (V_chol \%*\% B_chol)`.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @keywords internal
#' @note This is an internal non-exported function
#' @rdname latent_to_nwishart

latent_to_nwishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  if (!lower_chol) {
    V_chol <- Matrix::t(V_chol)
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
    list(W_chol = Matrix::t(W_chol), s = s_vec, B_chol = Matrix::t(B_chol))
  }
}

#' @title Transform Normalised Wishart to latent variables
#' @description Transform a Normalised Wishart matrix to latent iid
#' \eqn{N(0,1)} variables
#' @param W_chol The Cholesky factor of the Normalised Wishart matrix
#' @param V_chol The Cholesky factor of the Wishart \eqn{V} parameter
#' @param df The Wishart degrees of freedom
#' @param lower_chol `TRUE` if lower triangular Cholesky factors are used
#'   (default = `FALSE`)
#' @return A list with components
#'   \item{x}{A vector of latent variables identifying the Normalised Wishart
#'   matrix, length \eqn{d(d-1)/2}}
#'   \item{s}{Vector of scaling values}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When `lower_chol` is `FALSE`,
#' `W_chol = (B_chol \%*\% V_chol) \%*\% diag(s)`.
#'
#' When `lower_chol` is `TRUE`,
#' `W_chol = s * (V_chol \%*\% B_chol)`.
#'
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @keywords internal
#' @note This is an internal non-exported function
#' @rdname latent_to_nwishart

latent_from_nwishart <- function(W_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(W_chol)
  if (!lower_chol) {
    W_chol <- Matrix::t(W_chol)
    V_chol <- Matrix::t(V_chol)
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
    list(x = x, s = s_vec, B_chol = Matrix::t(B_chol))
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
#' @param lower_chol `TRUE` if lower triangular Cholesky factors are used
#'   (default = `FALSE`)
#' @return A list with components
#'   \item{W_chol}{The Cholesky factor of the Inverse Wishart matrix}
#'   \item{B_chol}{The inner Cholesky factor of the inverse Barlett
#'     decomposition}
#' When `lower_chol` is `FALSE`, `W_chol` is
#' `B_chol \%*\% V_chol`, otherwise `V_chol \%*\% B_chol`.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @keywords internal
#' @note This is an internal non-exported function
#' @rdname latent_to_iwishart

latent_to_iwishart <- function(x, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  rd <- rev(seq_len(d))
  # Compute Wishart matrix, in reverse ordering
  W <- latent_to_wishart(
    x = x,
    V_chol = Matrix::t(tri_solve(V_chol, lower_tri = lower_chol))[rd, rd],
    df = df,
    lower_chol = lower_chol
  )
  # Invert, and go back to original ordering
  W$W_chol <- Matrix::t(tri_solve(W$W_chol, lower_tri = lower_chol))[rd, rd]
  W$B_chol <- Matrix::t(tri_solve(W$B_chol, lower_tri = lower_chol))[rd, rd]
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
#' @keywords internal
#' @note This is an internal non-exported function
#' @rdname latent_to_iwishart

latent_from_iwishart <- function(W_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(V_chol)
  rd <- rev(seq_len(d))
  # Reverse order and invert:
  # Pre-normalise V_chol for numerical stability; renormalisation afterwards
  # means we don't need to store this rescaling
  W <- latent_from_wishart(
    Matrix::t(tri_solve(W_chol, lower_tri = lower_chol))[rd, rd],
    Matrix::t(tri_solve(V_chol, lower_tri = lower_chol))[rd, rd],
    df = df,
    lower_chol = lower_chol
  )
  list(
    x = W$x,
    B_chol = Matrix::t(tri_solve(W$B_chol, lower_tri = lower_chol))[rd, rd]
  )
}

# Normalised Inverse Wishart ####

normalise_ti_V_chol_rd <- function(V_chol, rd, lower_chol) {
  ti_V_chol_rd <- t(tri_solve(V_chol, lower_tri = lower_chol))[rd, rd]
  if (lower_chol) {
    ti_V_chol_rd <- (1 / rowSums(ti_V_chol_rd^2)^0.5) * ti_V_chol_rd
  } else {
    ti_V_chol_rd <-
      ti_V_chol_rd %*% Matrix::diag(1 / colSums(ti_V_chol_rd^2)^0.5)
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
#' @param lower_chol `TRUE` if lower triangular Cholesky factors are used
#'   (default = `FALSE`)
#' @return A list with components
#'   \item{W_chol}{The Cholesky factor of the Normalised Inverse Wishart matrix}
#'   \item{s}{Vector of scaling values}
#'   \item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
#' When `lower_chol` is `FALSE`,
#' `W_chol = (B_chol \%*\% V_chol) \%*\% diag(s)`.
#'
#' When `lower_chol` is `TRUE`,
#' `W_chol = s * (V_chol \%*\% B_chol)`.
#' @author Finn Lindgren
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @keywords internal
#' @note This is an internal non-exported function
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
  W$B_chol <- Matrix::t(tri_solve(W$B_chol, lower_tri = lower_chol))[rd, rd]
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
    W$W_chol <- W$W_chol %*% Matrix::diag(scale)
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
#' @keywords internal
#' @note This is an internal non-exported function
#' @rdname latent_to_niwishart

latent_from_niwishart <- function(W_chol, V_chol, df, lower_chol = FALSE) {
  d <- nrow(W_chol)
  rd <- rev(seq_len(d))
  # Reverse order and invert
  W <- latent_from_nwishart(
    Matrix::t(tri_solve(W_chol, lower_tri = lower_chol))[rd, rd],
    V_chol = normalise_ti_V_chol_rd(
      V_chol = V_chol,
      rd = rd,
      lower_chol = lower_chol
    ),
    df = df,
    lower_chol = lower_chol
  )
  W$B_chol <- Matrix::t(tri_solve(W$B_chol, lower_tri = lower_chol))[rd, rd]
  if (lower_chol) {
    W$s <- 1 / rowSums((V_chol %*% W$B_chol)^2)^0.5
  } else {
    W$s <- 1 / colSums((W$B_chol %*% V_chol)^2)^0.5
  }
  W
}



# Wishart density ####

#' @title Wishart Density
#' @description Calculate the density for a Wishart model
#' @param W A Wishart matrix, Default: NULL
#' @param x A latent parameter vector for Wishart model, Default: NULL
#' @param W_chol The Cholesky factor of a Wishart matrix, Default: NULL
#' @param V_chol The Cholesky factor of the V-parameter for a Wishart model
#' @param df The Wishart degrees of freedom parametre
#' @param lower_chol `TRUE` if lower triangular Cholesky factors are used,
#' Default: FALSE
#' @param log If `TRUE`, return the log-density, Default: FALSE
#' @return The density or log-density (if `log == TRUE`) for `W`, or the `W`
#' matrix constructed from a Cholesky factor `W_chol`
#' @details This function requires the `CholWishart::lmvgamma` function from
#' the `CholWishart` package.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @keywords internal
#' @note This is an internal non-exported function
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
      W_chol <- Matrix::t(W_chol)
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
    LVi_LW <- tri_solve(Matrix::t(V_chol), Matrix::t(W_chol), lower_tri = TRUE)
  }
  log_p <-
    as.vector(
      -df * d / 2 * log(2) - CholWishart::lmvgamma(df / 2, d) -
        df * sum(log(Matrix::diag(V_chol))) +
        (df - d - 1) * sum(log(Matrix::diag(W_chol))) -
        0.5 * sum(LVi_LW * LVi_LW)
    )
  if (log) {
    log_p
  } else {
    exp(log_p)
  }
}




# Generalised Wishart model object ####


#' @title Wishart-related Models
#' @description Defining a wrapper object class `wm_model` that can represent
#' Wishart, Normalised Wishart, Inverse Wishart, and Normalised Inverse Wishart
#' @param type Either `'wishart'`,  `'nwishart'`,  `'iwishart'`,  `'niwishart'`
#' @param V The matrix parameter for the distribution
#' @param df The degrees-of-freedom parameter for the distribution
#' @param V_chol The Cholesky factor of `V`. Type must match the `lower_chol`
#' parameter
#' @param lower_chol logical; For `wm_model`, whether the internal
#' representation should use lower triangular Cholesky factors. For other methods, determines what
#' Cholesky type input is or output should be returned as, with `NULL`
#' inheriting the internal `model` setting.
#' Default for `wm_model` is `FALSE`, for other methods default is `NULL`.
#' @return [wm_model()] returns a `wm_model` object that encapsulates the
#' parameters of one of the four Wishart model types, as defined by `type`.
#' @export
#' @rdname wishart_model

wm_model <- function(type,
                     V = NULL, df = NULL,
                     V_chol = NULL, lower_chol = FALSE) {
  type <- match.arg(type, c("wishart", "nwishart", "iwishart", "niwishart"))
  if (!is.null(V)) {
    V <- as.matrix(V)
    stopifnot(is.null(V_chol))
    V_chol <- chol(V)
    if (lower_chol) {
      V_chol <- Matrix::t(V_chol)
    }
  } else {
    V_chol <- as.matrix(V_chol)
    stopifnot(!is.null(V_chol))
    if (lower_chol) {
      V <- V_chol %*% Matrix::t(V_chol)
    } else {
      V <- Matrix::t(V_chol) %*% V_chol
    }
  }
  d <- nrow(V)
  stopifnot(
    "df most not be NULL" = !is.null(df),
    "df should be > dimension - 1" = df > d - 1
  )
  if (type %in% c("wishart", "iwishart")) {
    N_latent <- d * (d + 1) / 2
  } else {
    N_latent <- d * (d - 1) / 2
  }

  model <- list(
    type = type,
    d = d, N_latent = N_latent, V = V, df = df,
    V_chol = V_chol, lower_chol = lower_chol
  )
  class(model) <-
    c(paste0("wm_model_", type), "wm_model")
  model
}


#' @param model A `wm_model` object
#' @param W A symmetric matrix valid for the value space of the `model` type
#' @param W_chol The Cholesky factor of a matrix valid for the value space of
#' the `model` type. Only one of `W` and `W_chol` may be given.
#' @param ... Further parameters passed on to other methods
#' @return [wm_latent()] returns the latent variables for the representation
#' of a (W/NW/IW/NIW) matrix, given either the matrix itself in `W`, or its
#' Cholesky factor in `W_chol`.
#' @export
#' @rdname wishart_model

wm_latent <-
  function(model, W = NULL, W_chol = NULL, lower_chol = NULL, ...) {
    stopifnot(inherits(model, "wm_model"))
    if (is.null(lower_chol)) {
      lower_chol <- model$lower_chol
    }
    if (!is.null(W)) {
      W_chol <- chol(W)
      if (model$lower_chol) {
        W_chol <- Matrix::t(W_chol)
      }
    } else {
      stopifnot(!is.null(W_chol))
      if (lower_chol != model$lower_chol) {
        W_chol <- Matrix::t(W_chol)
      }
    }
    switch(model$type,
      "wishart" = latent_from_wishart(
        W_chol = W_chol, df = model$df, V_chol = model$V_chol,
        lower_chol = model$lower_chol
      ),
      "nwishart" = latent_from_nwishart(
        W_chol = W_chol, df = model$df, V_chol = model$V_chol,
        lower_chol = model$lower_chol
      ),
      "iwishart" = latent_from_iwishart(
        W_chol = W_chol, df = model$df, V_chol = model$V_chol,
        lower_chol = model$lower_chol
      ),
      "niwishart" = latent_from_niwishart(
        W_chol = W_chol, df = model$df, V_chol = model$V_chol,
        lower_chol = model$lower_chol
      )
    )$x
  }

#' @param latent A numeric vector of length `model$N_latent` for the latent
#' representation of a model outcome
#' @return [wm_chol()] returns the Cholesky factor of a (W/NW/IW/NIW) matrix.
#' @export
#' @rdname wishart_model

wm_chol <-
  function(model, latent, lower_chol = NULL, ...) {
    stopifnot(inherits(model, "wm_model"))
    if (is.null(lower_chol)) {
      lower_chol <- model$lower_chol
    }
    W_chol <-
      switch(model$type,
        "wishart" = latent_to_wishart(
          x = latent, df = model$df, V_chol = model$V_chol,
          lower_chol = model$lower_chol
        ),
        "nwishart" = latent_to_nwishart(
          x = latent, df = model$df, V_chol = model$V_chol,
          lower_chol = model$lower_chol
        ),
        "iwishart" = latent_to_iwishart(
          x = latent, df = model$df, V_chol = model$V_chol,
          lower_chol = model$lower_chol
        ),
        "niwishart" = latent_to_niwishart(
          x = latent, df = model$df, V_chol = model$V_chol,
          lower_chol = model$lower_chol
        )
      )$W_chol
    if (lower_chol != model$lower_chol) {
      W_chol <- Matrix::t(W_chol)
    }
    W_chol
  }

#' @return [wm_matrix()] returns a (W/NW/IW/NIW) matrix.
#' @export
#' @rdname wishart_model

wm_matrix <-
  function(model, latent, ...) {
    stopifnot(inherits(model, "wm_model"))
    W_chol <- wm_chol(model, latent = latent, lower_chol = FALSE)
    Matrix::t(W_chol) %*% W_chol
  }


#' @param symmetric logical; If `TRUE`, use symmetric finite differences to
#' compute derivatives
#' @param h positive delta for finite differences
#' @return [wm_matrix_jacobian()] returns the Jacobian for the
#' (column-)vectorised (see [cvec()]) matrix with respect to the latent
#' variables.
#'
#' @export
#' @rdname wishart_model

wm_matrix_jacobian <-
  function(model, latent, symmetric = TRUE, h = 1e-4, ...) {
    stopifnot(inherits(model, "wm_model"))
    jacobian <- Matrix::Matrix(0, model$d * model$d, model$N_latent)
    if (symmetric) {
      for (loop in seq_len(model$N_latent)) {
        H <- rep(c(0, h, 0), times = c(loop - 1, 1, model$N_latent - loop))
        W_p <- wm_matrix(model, latent + H)
        W_m <- wm_matrix(model, latent - H)
        jacobian[, loop] <- cvec(W_p - W_m) / (2 * h)
      }
    } else {
      W <- wm_matrix(model, latent)
      for (loop in seq_len(model$N_latent)) {
        H <- rep(c(0, h, 0), times = c(loop - 1, 1, model$N_latent - loop))
        W_p <- wm_matrix(model, latent + H)
        jacobian[, loop] <- cvec(W_p - W) / h
      }
    }
    jacobian
  }


#' @return [wm_chol_jacobian()] returns the Jacobian for the
#' column-vectorised (see [cvec()]) Cholesky matrix with respect to the
#' latent variables, with type determined by the `lower_chol` setting.
#' @export
#' @rdname wishart_model

wm_chol_jacobian <-
  function(model, latent, symmetric = TRUE, h = 1e-4, lower_chol = NULL, ...) {
    stopifnot(inherits(model, "wm_model"))
    jacobian <- Matrix::Matrix(0, model$d * model$d, model$N_latent)
    if (symmetric) {
      for (loop in seq_len(model$N_latent)) {
        H <- rep(c(0, h, 0), times = c(loop - 1, 1, model$N_latent - loop))
        W_p <- wm_chol(model, latent + H, lower_chol = lower_chol)
        W_m <- wm_chol(model, latent - H)
        jacobian[, loop] <- cvec(W_p - W_m) / (2 * h)
      }
    } else {
      W <- wm_chol(model, latent, lower_chol = lower_chol)
      for (loop in seq_len(model$N_latent)) {
        H <- rep(c(0, h, 0), times = c(loop - 1, 1, model$N_latent - loop))
        W_p <- wm_chol(model, latent + H, lower_chol = lower_chol)
        jacobian[, loop] <- cvec(W_p - W) / h
      }
    }
    jacobian
  }



#' @param mean_latent Expectation vector for the latent variables
#' @param cov_latent Covariance matrix for the latent variables
#' @param order A vector of two integers defining the Taylor expansion orders
#' used for linearised moment calculations for expectation and variance,
#' respectively. Can be either \eqn{1} or \eqn{2}. Default: `c(2, 1)`
#' @return [wm_moments_linear()] returns linearised approximations of propagated
#' mean and standard deviation of the Wishart matrix given mean and covariance
#' of the latent variables.
#' @details For `wm_moments_linear`, the further `...` parameters are passed
#' on to `wm_matrix_jacobian`.
#'
#' @export
#' @rdname wishart_model

wm_moments_linear <- function(model,
                              mean_latent = rep(0, model$N_latent),
                              cov_latent = sparse_identity(model$N_latent),
                              order = c(2, 1),
                              h = 1e-4,
                              ...) {
  stopifnot(inherits(model, "wm_model"))
  m <- wm_matrix(model, latent = mean_latent)
  jacobian <- wm_matrix_jacobian(model, latent = mean_latent, h = h, ...)
  S <-
    cvec(
      Matrix::rowSums(
        jacobian * (jacobian %*% cov_latent)
      ),
      d = model$d,
      sparse = FALSE
    )
  if (any(order > 1)) {
    for (k in seq_len(model$N_latent)) {
      H <- rep(c(0, h, 0), times = c(k - 1, 1, model$N_latent - k))
      jacobian_k <-
        (wm_matrix_jacobian(model, latent = mean_latent + H, h = h, ...) -
          wm_matrix_jacobian(model, latent = mean_latent - H, h = h, ...)) /
        (2 * h)
      if (order[1] > 1) {
        m <- m + 0.5 * cvec(jacobian_k %*% cov_latent[, k],
          d = model$d,
          sparse = FALSE
        )
      }
      if (order[2] > 1) {
        S <- S + 0.5 * cvec(Matrix::rowSums((jacobian_k %*% cov_latent)^2),
          d = model$d,
          sparse = FALSE
        )
      }
    }
  }
  list(mean = m, sd = sqrt(S))
}





#' @param log If `TRUE`, return the log-density, Default: FALSE
#' @return `wm_density()` returns the density or log-density (if `log == TRUE`)
#' for `W`, or the `W` matrix constructed from a Cholesky factor `W_chol`
#' @details The `wm_density()` method requires the `CholWishart::lmvgamma`
#' function from the `CholWishart` package.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname wishart_model

wm_density <- function(model,
                       latent = NULL,
                       W = NULL, W_chol = NULL,
                       lower_chol = NULL, log = FALSE) {
  stopifnot(inherits(model, "wm_model"))
  stopifnot("Wishart density only implemented for plain Wishart models" =
              (model[["type"]] == "wishart"))
  if (is.null(lower_chol)) {
    lower_chol <- model$lower_chol
  }
  if (!is.null(W_chol) && (lower_chol != model$lower_chol)) {
    W_chol <- t(W_chol)
  }
  dwishart(W = W, x = latent, W_chol = W_chol,
           V_chol = model$V_chol, V_chol$df,
           lower_chol = model$lower_chol, log = log)
}
