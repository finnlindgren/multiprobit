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



#' @title Inverse Reverse Cholesky
#' @description Compute inverse of a Cholesky matrix,
#' storing the result in reverse order, see Details.
#' @param Sigma_chol Holesky factor of a covariance matrix
#' @param lower_chol `TRUE` if lower triangular Cholesky factors are used
#'   Default: FALSE
#' @return A Cholesky matrix of the same triangle orientation as the input,
#'   with rows and columns in reverse order.
#' @details For `lower_chol == FALSE`,
#' the input `Sigma_chol` is the matrix \eqn{R} in the Cholesky factorisation
#' \eqn{\Sigma=R^T R} with inverse \eqn{Q=R^{-1}R^{-T}}. Since \eqn{R^{-1}} has
#' the opposite upper/lower triangular property to \eqn{R^T}, this \eqn{Q}
#' factorisation isn't a regular Cholesky factorisation. Let \eqn{P} be the
#' permutation matrix that reverses element order. Then
#' \eqn{PQP=PR^{-1}PPR^{-T}P}, and \eqn{PR^{-T}P} is the Cholesky factor of
#' \eqn{PQP} with the same upper/lower triangular property as \eqn{R}.
#'
#' For `lower_chol == TRUE`, the input is \eqn{L} in \eqn{\Sigma=LL^T} and the
#' output is \eqn{PL^{-T}P}, with \eqn{PQP=PL^{-T}PPL^{-1}P}.
#' @examples
#' if(interactive()){
#'   inverse_chol_reverse(matrix(c(1, 0, 2, 3), 2, 2))
#' }
#' @keywords internal
#' @rdname inverse_chol_reverse

inverse_chol_reverse <- function(Sigma_chol, lower_chol = FALSE) {
  d <- nrow(Sigma_chol)
  Q_chol <- tri_solve(Sigma_chol, lower_tri = lower_chol)
  perm <- rev(seq_len(d))
  Matrix::t(Q_chol)[perm, perm]
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
