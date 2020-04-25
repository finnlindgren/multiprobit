#' @title Internal multiprobit utilities
#' @details `tri_solve` solves triangular systems with back/forwardsolve
#' @keywords internal
#' @rdname multiprobit_utils

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


#' @details `qchisq_pnorm` evaluates `qchisq(pnorm(...))` with
#' attempt at numerical stability.
#' @keywords internal
#' @importFrom stats qchisq pnorm
#' @rdname multiprobit_utils

qchisq_pnorm <- function(x, df) {
  qchisq(pnorm(x, log.p = TRUE),
         df = df, log.p = TRUE
  )
}

#' @keywords internal
#' @importFrom stats pchisq qnorm
#' @rdname multiprobit_utils

qnorm_pchisq <- function(x, df) {
  qnorm(pchisq(x, df = df, log.p = TRUE), log.p = TRUE)
}

#' @keywords internal
#' @importFrom stats qbeta pchisq
#' @rdname multiprobit_utils

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
#' @rdname multiprobit_utils

qchisq_pbeta <- function(x,
                         shape1, shape2,
                         df,
                         ncp) {
  qchisq(
    pbeta(x, shape1 = shape1, shape2 = shape2, ncp = ncp, log.p = TRUE),
    df = df, ncp = ncp, log.p = TRUE
  )
}


#' @rdname multiprobit_utils

cvec <- function(A, d = NULL) {
  if (is.matrix(A)) {
    as.vector(A)
  } else {
    matrix(A, length(A) / d, d, byrow = FALSE)
  }
}

#' @rdname multiprobit_utils

rvec <- function(A, d = NULL) {
  if (is.matrix(A)) {
    as.vector(t(A))
  } else {
    matrix(A, length(A) / d, d, byrow = TRUE)
  }
}


#' @importFrom Matrix sparseMatrix
#' @rdname multiprobit_utils

sparse_identity <- function(d) {
  Matrix::sparseMatrix(
    i = seq_len(d),
    j = seq_len(d),
    x = 1.0,
    dims = c(d, d)
  )
}
