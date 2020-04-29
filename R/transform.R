

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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @details This function requires the `CholWishart::lmvgamma` function from
#' the `CholWishart` package.
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


#' @title Wishart model wrapper
#' @param type DOC
#' @param V DOC
#' @param df DOC
#' @param V_chol DOC
#' @param lower_chol DOC
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
    N_u <- d * (d + 1) / 2
  } else {
    N_u <- d * (d - 1) / 2
  }

  model <- list(
    type = type,
    d = d, N_u = N_u, V = V, df = df,
    V_chol = V_chol, lower_chol = lower_chol
  )
  class(model) <-
    c(paste0("wm_model_", type), "wm_model")
  model
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param model PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname wishart_model

wm_latent <- function(model, ...) {
  UseMethod("wm_latent")
}

#' @export
#' @rdname wishart_model

wm_chol <- function(model, ...) {
  UseMethod("wm_chol")
}

#' @export
#' @rdname wishart_model

wm_matrix <- function(model, ...) {
  UseMethod("wm_matrix")
}

#' @export
#' @rdname wishart_model

wm_matrix_jacobian <- function(model, ...) {
  UseMethod("wm_matrix_jacobian")
}

#' @param W DOC
#' @param W_chol DOC
#' @param lower_chol DOC
#' @export
#' @rdname wishart_model

wm_latent.wm_model <-
  function(model, W = NULL, W_chol = NULL, lower_chol = NULL, ...) {
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

#' @param x DOC
#' @export
#' @rdname wishart_model

wm_chol.wm_model <-
  function(model, x, lower_chol = NULL, ...) {
    if (is.null(lower_chol)) {
      lower_chol <- model$lower_chol
    }
    W_chol <-
      switch(model$type,
        "wishart" = latent_to_wishart(
          x = x, df = model$df, V_chol = model$V_chol,
          lower_chol = model$lower_chol
        ),
        "nwishart" = latent_to_nwishart(
          x = x, df = model$df, V_chol = model$V_chol,
          lower_chol = model$lower_chol
        ),
        "iwishart" = latent_to_iwishart(
          x = x, df = model$df, V_chol = model$V_chol,
          lower_chol = model$lower_chol
        ),
        "niwishart" = latent_to_niwishart(
          x = x, df = model$df, V_chol = model$V_chol,
          lower_chol = model$lower_chol
        )
      )$W_chol
    if (lower_chol != model$lower_chol) {
      W_chol <- Matrix::t(W_chol)
    }
    W_chol
  }

#' @export
#' @rdname wishart_model

wm_matrix.wm_model <-
  function(model, x, ...) {
    W_chol <- wm_chol(model, x = x, lower_chol = FALSE)
    Matrix::t(W_chol) %*% W_chol
  }


#' @param symmetric logical; If `TRUE`, use symmetric finite differences to
#' compute derivatives
#' @param h positive delta for finite differences
#' @return `wm_matrix_jacobian` returns the Jacobian for the (column-)vectorised
#' (see [cvec()]) matrix with respect to the latent variables.
#'
#' @export
#' @rdname wishart_model

wm_matrix_jacobian.wm_model <-
  function(model, x, symmetric = TRUE, h = 1e-4, ...) {
    jacobian <- Matrix::Matrix(0, model$d * model$d, model$N_u)
    if (symmetric) {
      for (loop in seq_len(model$N_u)) {
        H <- rep(c(0, h, 0), times = c(loop - 1, 1, model$N_u - loop))
        W_p <- wm_matrix(model, x + H)
        W_m <- wm_matrix(model, x - H)
        jacobian[, loop] <- cvec(W_p - W_m) / (2 * h)
      }
    } else {
      W <- wm_matrix(model, x)
      for (loop in seq_len(model$N_u)) {
        H <- rep(c(0, h, 0), times = c(loop - 1, 1, model$N_u - loop))
        W_p <- wm_matrix(model, x + H)
        jacobian[, loop] <- cvec(W_p - W) / h
      }
    }
    jacobian
  }



#' @param symmetric logical; If `TRUE`, use symmetric finite differences to
#' compute derivatives
#' @param h positive delta for finite differences
#' @return `wm_chol_jacobian` returns the Jacobian for the column-vectorised
#' (see [cvec()]) Cholesky matrix (lower Cholesky if `lower_chol` is `TRUE`,
#' upper if `lower_chol` is `FALSE`, and for the model internal version if
#' `lower_chol` is `NULL`) with respect to the latent
#'  variables.
#'
#' @export
#' @rdname wishart_model

wm_chol_jacobian.wm_model <-
  function(model, x, symmetric = TRUE, h = 1e-4, lower_chol = NULL, ...) {
    jacobian <- Matrix::Matrix(0, model$d * model$d, model$N_u)
    if (symmetric) {
      for (loop in seq_len(model$N_u)) {
        H <- rep(c(0, h, 0), times = c(loop - 1, 1, model$N_u - loop))
        W_p <- wm_chol(model, x + H, lower_chol = lower_chol)
        W_m <- wm_chol(model, x - H)
        jacobian[, loop] <- cvec(W_p - W_m) / (2 * h)
      }
    } else {
      W <- wm_chol(model, x, lower_chol = lower_chol)
      for (loop in seq_len(model$N_u)) {
        H <- rep(c(0, h, 0), times = c(loop - 1, 1, model$N_u - loop))
        W_p <- wm_chol(model, x + H, lower_chol = lower_chol)
        jacobian[, loop] <- cvec(W_p - W) / h
      }
    }
    jacobian
  }
