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


#' @title Vectorisation utitilies
#' @description Convert between matrix and vector representations
#' @param A A matrix or a vector
#' @param d The number of matrix columns
#' @param sparse logical; If `TRUE`, use `Matrix::Matrix` to construct matrix
#' output, that detects sparsity. Default: `FALSE`
#' @return
#'   If `A` is a matrix, `cvec` returns the columnwise vectorisation of `A`.
#'   and `rvec` returns the rowwise vectorisation of `A`.
#'   If `A` is a vector, `cvec` and `rvec` return a matrix with `d` columns,
#'   filled columnwise for `cvec` and rowwise for `rvec`.
#' @examples
#' if (interactive()) {
#'   Ac <- cvec(1:6, 2) # As matrix(1:6, 3, 2, byrow = FALSE)
#'   Ar <- rvec(1:6, 2) # As matrix(1:6, 3, 2, byrow = TRUE)
#'   cvec(Ac) # As as.vector(Ac)
#'   rvec(Ar) # As as.vector(t(Ar))
#' }
#' @export
#' @rdname vectorisation

cvec <- function(A, d = NULL, sparse = FALSE) {
  if (is.matrix(A)) {
    as.vector(A)
  } else if (sparse) {
    Matrix::Matrix(as.vector(A), length(A) / d, d, byrow = FALSE)
  } else {
    matrix(A, length(A) / d, d, byrow = FALSE)
  }
}

#' @export
#' @rdname vectorisation

rvec <- function(A, d = NULL, sparse = FALSE) {
  if (is.matrix(A)) {
    as.vector(t(A))
  } else if (sparse) {
    Matrix::Matrix(as.vector(A), length(A) / d, d, byrow = TRUE)
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
