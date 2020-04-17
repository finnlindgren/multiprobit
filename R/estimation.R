cvec <- function(A, d = NULL) {
  if (is.matrix(A)) {
    as.vector(A)
  } else {
    matrix(A, length(A) / d, d, byrow = FALSE)
  }
}
rvec <- function(A, d = NULL) {
  if (is.matrix(A)) {
    as.vector(t(A))
  } else {
    matrix(A, length(A) / d, d, byrow = TRUE)
  }
}
myident <- function(d) {
  Matrix::sparseMatrix(
    i = seq_len(d),
    j = seq_len(d),
    x = 1.0,
    dims = c(d, d)
  )
}


mp_loglike <- function(Y, mu, Sigma_chol,
                       lower_chol = FALSE, ...) {
  sum(mpp(Y,
    mu = mu, Sigma_chol = Sigma_chol,
    lower_chol = lower_chol,
    ..., log = TRUE
  )$P)
}
mp_loglike_gradient_beta <- function(Y, X, mu, Sigma_chol, lower_chol = FALSE,
                                     ...) {
  Q_chol <- inverse_chol_reverse(Sigma_chol, lower_chol = lower_chol)
  d <- ncol(Y)
  perm <- rev(seq_len(d))
  L <- 0
  for (i in seq_len(nrow(Y))) {
    L <- L + kronecker(
      t(X[i, , drop = FALSE]),
      mpp_gradient_mu(Y[i, perm, drop = FALSE],
        mu = mu[i, perm, drop = FALSE],
        Q_chol = Q_chol,
        lower_chol = lower_chol,
        ..., log = TRUE
      )[perm]
    )
  }
  L
}
mp_loglike_gradient_u <- function(Y, u, mu, V_chol, df, ...) {
  mpp_gradient_u(Y,
    u = u, mu = mu,
    V_chol = V_chol, df = df, ..., log = TRUE
  )
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Y PARAM_DESCRIPTION
#' @param X PARAM_DESCRIPTION
#' @param V_chol PARAM_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param prec_beta PARAM_DESCRIPTION
#' @param beta PARAM_DESCRIPTION
#' @param u PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @param lower_chol PARAM_DESCRIPTION
#' @param what PARAM_DESCRIPTION, Options: "loglike", "grad", with additional
#'   options "grad_beta", "grad_u" for \code{mp_logposterior} only
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_logposterior

mp_logposterior <- function(Y, X, V_chol, df, prec_beta,
                            beta, u, ...,
                            lower_chol = FALSE,
                            what = c("loglike", "grad", "grad_beta", "grad_u")) {
  what <- match.arg(what)
  d <- ncol(Y)
  dd <- d * (d - 1) / 2
  J <- ncol(X)

  stopifnot(length(beta) + length(u) == (J * d + dd))

  mu <- X %*% rvec(beta, d)
  if (what != "grad_u") {
    Sigma_chol <- latent_to_nwishart(
      x = u, V_chol = V_chol, df = df,
      lower_chol = lower_chol
    )$W_chol
  }
  if (what == "loglike") {
    L <- -sum(u)^2 / 2 - prec_beta / 2 * sum(beta^2) +
      mp_loglike(
        Y = Y, mu = mu, Sigma_chol = Sigma_chol,
        lower_chol = lower_chol,
        ...
      )
  } else {
    if (what %in% c("grad", "grad_beta")) {
      dL_dbeta <- -prec_beta * beta +
        mp_loglike_gradient_beta(
          Y = Y, X = X, mu = mu,
          Sigma_chol = Sigma_chol,
          lower_chol = lower_chol,
          ...
        )
    }
    if (what %in% c("grad", "grad_u")) {
      dL_du <- -u + mp_loglike_gradient_u(
        Y = Y, u = u, mu = mu,
        V_chol = V_chol, df = df,
        lower_chol = lower_chol,
        ...
      )
    }
    switch(what,
      "grad" = L <- c(dL_dbeta, dL_du),
      "grad_beta" = L <- dL_dbeta,
      "grad_u" = L <- dL_du
    )
  }
  L
}

#' @param latent PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_logposterior

mp_logposterior_joint <- function(latent, Y, X, V_chol, df, prec_beta, ...,
                                  lower_chol = FALSE,
                                  what = c("loglike", "grad")) {
  what <- match.arg(what)

  d <- ncol(Y)
  dd <- d * (d - 1) / 2
  J <- ncol(X)

  stopifnot(length(latent) == (J * d + dd))

  index_beta <- seq_len(J * d)
  index_u <- J * d + seq_len(dd)

  mp_logposterior(Y, X, V_chol, df, prec_beta,
    beta = latent[index_beta],
    u = latent[index_u],
    lower_chol = lower_chol,
    what = what,
    ...
  )
}

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_logposterior

mp_logposterior_fixed_beta <- function(latent, Y, X, V_chol, df, prec_beta,
                                       beta, ...,
                                       lower_chol = FALSE,
                                       what = c("loglike", "grad")) {
  what <- match.arg(what)

  d <- ncol(Y)
  dd <- d * (d - 1) / 2

  stopifnot(length(latent) == (dd))

  mp_logposterior(Y, X, V_chol, df, prec_beta,
    beta = beta,
    u = latent,
    lower_chol = lower_chol,
    what = ifelse(what == "loglike", "loglike", "grad_u"),
    ...
  )
}

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_logposterior

mp_logposterior_fixed_u <- function(latent, Y, X, V_chol, df, prec_beta,
                                    u, ...,
                                    lower_chol = FALSE,
                                    what = c("loglike", "grad")) {
  what <- match.arg(what)

  d <- ncol(Y)
  J <- ncol(X)

  stopifnot(length(latent) == (J * d))

  mp_logposterior(Y, X, V_chol, df, prec_beta,
    beta = latent,
    u = u,
    lower_chol = lower_chol,
    what = ifelse(what == "loglike", "loglike", "grad_beta")
  )
}



#' @title Estimate a multivariate probit model
#' @description Estimate a multivariate probit model from multivariate binary
#'   data in a Bayesian generalised linear model framework
#' @param response A matrix with n-by-d elements, where each row is a multivariate
#'   observation, see Details. A vector is interpreted as a single row matrix.
#' @param X An optimally precomputed n-by-J model matrix, where J is the number
#'   regression coeficcients for each of the d dimensions.
#' @param formula A formula interpretable by \code{model.matrix}.
#' @param data A \code{data.frame} containing the veriables needed by the
#'   formula.
#' @param df Degrees of freedom for the normalised Wishart prior for the
#'   correlation matrix. See Details.
#' @param prec_beta Prior precision for the regression coefficients
#' @param gaussint_options Optional list of options for
#'   \code{excursions::gaussint}.
#'   Specific relevant option:
#'   \code{num_threads} The maximum number of allowed threads for parallel
#' computing by \code{excursions::gaussint}, Default: 0, meaning no limit.
#' @param optim_options Optional list of control options for \code{optim()}
#' @param method Optimisation method. Options: "alternating", "joint",
#'   "stepwise"
#' @param max_iter The maximum number of steps for
#'   \code{method == "alternating"}
#' @return OUTPUT_DESCRIPTION
#' @details Estimates a multivariate probit model for the response variable,
#' here denoted \eqn{Y}. The model is built from a linear predictor
#' \deqn{M = X B}
#' where \eqn{X} is a n-by-J matrix of \eqn{J} predictors, and \eqn{B} is
#' a J-by-d matrix of regression coefficients.
#' Each row of \eqn{M} is the linear
#' predictor for one multivariate observation. The response variables \eqn{Y}
#' are linked to \eqn{M} by first defining latent Gaussian variables
#' \deqn{Z=M+E} where each row of \eqn{E} is a multivariate Normal vector,
#' \eqn{E ~ N(0,\Sigma)}. Then, \deqn{Y_{i,k}=I(Z_{i,k} > 0).}
#' Conditionally on \eqn{B}, each row of \eqn{Y} has a multinomial distribution
#' on the set of all \eqn{0/1} combinations, with each probability equal to a
#' hyperquadrant probability of a the multivariate Normal distribution
#' \eqn{N(\mu,\Sigma)}, where \eqn{\mu} is the corresponding row of \eqn{M}.
#' Only the inequality \eqn{Y_{i,k} > 0} for the response variables is used,
#' so alternative data representations such as \eqn{-1/+1} will also work as
#' expected.
#'
#' The degrees of freedom for the normalised Wishart prior are linked to the
#' LKJ prior by \code{df = ...}, which makes the two models equivalent.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @importFrom stats model.matrix
#' @importFrom stats optim
#' @rdname multiprobit

multiprobit <- function(response, X = NULL,
                        formula = NULL, data = NULL,
                        df, prec_beta = 0,
                        gaussint_options = NULL,
                        optim_options = NULL,
                        max_iter = 10,
                        method = c("alternating", "joint", "stepwise")) {
  Y <- response > 0
  if (!is.null(X)) {
    if (!is.null(formula) || !is.null(data)) {
      stop("When X is provided, do not also provide formula & data")
    }
    X <- as.matrix(X)
  } else {
    if (is.null(formula) || is.null(data)) {
      stop("Either X or formula & data must be provided")
    }
    X <- model.matrix(formula, data = data)
  }

  stopifnot(nrow(X) == nrow(Y))
  V_chol <- myident(ncol(Y))
  model <- list(
    Y = Y, X = X, V_chol = V_chol, df = df, prec_beta = prec_beta)

  method <- match.arg(method)
  if (method == "alternating") {
    opt <- optim_alternating(model,
                             max_iter = max_iter,
                             optim_options = optim_options,
                             gaussint_options = gaussint_options)
  } else if (method == "joint") {
    opt <- optim_joint(model,
                       optim_options = optim_options,
                       gaussint_options = gaussint_options)
  } else if (method == "stepwise") {
    opt <- optim_stepwise(model,
                          optim_options = optim_options,
                          gaussint_options = gaussint_options)
  } else {
    stop("Unknown method; This cannot happen.")
  }

  opt
}











fn <- function(x, model, ..., opt_type = c("joint", "beta", "u")) {
  opt_type <- match.arg(opt_type)
  what <- "loglike"
  switch(opt_type,
         "joint" = do.call(mp_logposterior_joint,
                           c(list(latent = x, what = what, ...),
                             model)),
         "beta" = do.call(mp_logposterior_fixed_u,
                          c(list(latent = x, what = what, ...),
                            model)),
         "u" = do.call(mp_logposterior_fixed_beta,
                       c(list(latent = x, what = what, ...),
                         model))
  )
}
gr <- function(x, model, ..., opt_type = c("joint", "beta", "u")) {
  opt_type <- match.arg(opt_type)
  what <- "grad"
  switch(opt_type,
         "joint" = do.call(mp_logposterior_joint,
                           c(list(latent = x, what = what, ...),
                             model)),
         "beta" = do.call(mp_logposterior_fixed_u,
                          c(list(latent = x, what = what, ...),
                            model)),
         "u" = do.call(mp_logposterior_fixed_beta,
                       c(list(latent = x, what = what, ...),
                         model))
  )
}


optim_alternating <- function(model, max_iter,
                              optim_options = NULL,
                              gaussint_options = NULL) {
  optim_options <- as.list(optim_options)
  optim_options[["fnscale"]] <- -1
  if (is.null(optim_options[["method"]])) {
    method <- "BFGS"
  } else {
    method <- optim_options[["method"]]
    optim_options[["method"]] <- NULL
  }
  gaussint_options <- as.list(gaussint_options)
  if (is.null(gaussint_options[["seed"]])) {
    gaussint_options[["seed"]] <- 1L
  }

  N <- nrow(model$Y)
  d <- ncol(model$Y)
  J <- ncol(model$X)
  N_beta <- J * d
  N_u <- d * (d - 1) / 2
  N_latent <- N_beta + N_u
  opt_beta = list()
  opt_u = list()
  result <- data.frame(index = rep(seq_len(N_latent),
                                   times = max_iter),
                       latent = 0,
                       iteration = rep(seq_len(max_iter), each = N_latent),
                       type = "alternating")
  counts <- data.frame(iteration = seq_len(max_iter),
                       fn = numeric(max_iter),
                       gr_beta = numeric(max_iter),
                       gr_u = numeric(max_iter))
  beta <- rep(0, N_beta)
  u <- rep(0, N_u)
  for (loop in seq_len(max_iter)) {
    if (loop > 1) {
      counts[loop, c("fn", "gr_beta", "gr_u")] <-
        counts[loop - 1, c("fn", "gr_beta", "gr_u")]
    }
    opt_beta[[loop]] <-
      optim(par = beta,
            fn = fn, gr = gr,
            opt_type = "beta",
            model = model,
            u = u,
            gaussint_options = gaussint_options,
            method = method,
            control = optim_options)
    beta <- opt_beta[[loop]]$par
    counts$fn[loop] <- counts$fn[loop] + opt_beta[[loop]]$counts["function"]
    counts$gr_beta[loop] <- counts$gr_beta[loop] + opt_beta[[loop]]$counts["gradient"]
    opt_u[[loop]] <-
      optim(par = u,
            fn = fn, gr = gr,
            opt_type = "u",
            model = model,
            beta = beta,
            gaussint_options = gaussint_options,
            method = method,
            control = optim_options)
    u <- opt_u[[loop]]$par
    counts$fn[loop] <- counts$fn[loop] + opt_u[[loop]]$counts["function"]
    counts$gr_u[loop] <- counts$gr_u[loop] + opt_u[[loop]]$counts["gradient"]

    result$latent[(loop - 1) * N_latent + seq_len(N_latent)] <- c(beta, u)
  }
  list(result = result, counts = counts, opt_beta = opt_beta, opt_u = opt_u)
}

optim_joint <- function(model, optim_options = NULL, gaussint_options = NULL) {
  optim_options <- as.list(optim_options)
  optim_options[["fnscale"]] <- -1
  if (is.null(optim_options[["method"]])) {
    method <- "BFGS"
  } else {
    method <- optim_options[["method"]]
    optim_options[["method"]] <- NULL
  }
  gaussint_options <- as.list(gaussint_options)
  if (is.null(gaussint_options[["seed"]])) {
    gaussint_options[["seed"]] <- 1L
  }

  N <- nrow(model$Y)
  d <- ncol(model$Y)
  J <- ncol(model$X)
  N <- nrow(model$Y)
  d <- ncol(model$Y)
  J <- ncol(model$X)
  N_beta <- J * d
  N_u <- d * (d - 1) / 2
  N_latent <- N_beta + N_u
  beta <- rep(0, N_beta)
  u <- rep(0, N_u)
  opt_joint <-
    optim(par = c(beta, u),
          fn = fn, gr = gr,
          opt_type = "joint",
          model = model,
          gaussint_options = gaussint_options,
          method = method,
          control = optim_options)
  result <- data.frame(index = seq_len(N_latent),
                       latent = opt_joint$par,
                       iteration = 1,
                       type = "joint")
  counts <- data.frame(iteration = 1,
                       fn = opt_joint$counts["function"],
                       gr_beta = opt_joint$counts["gradient"],
                       gr_u = opt_joint$counts["gradient"])
  list(result = result, counts = counts, opt_joint = opt_joint)
}


optim_stepwise <- function(model, optim_options = NULL, gaussint_options = NULL) {
  optim_options <- as.list(optim_options)
  optim_options[["fnscale"]] <- -1
  if (is.null(optim_options[["method"]])) {
    method <- "BFGS"
  } else {
    method <- optim_options[["method"]]
    optim_options[["method"]] <- NULL
  }
  gaussint_options <- as.list(gaussint_options)
  if (is.null(gaussint_options[["seed"]])) {
    gaussint_options[["seed"]] <- 1L
  }

  N <- nrow(model$Y)
  d <- ncol(model$Y)
  J <- ncol(model$X)
  N_beta <- J * d
  N_u <- d * (d - 1) / 2
  N_latent <- N_beta + N_u
  opt_beta = list()
  opt_u = list()
  beta <- rep(0, N_beta)
  u <- rep(0, N_u)
  opt_beta <-
    optim(par = beta,
          fn = fn, gr = gr,
          opt_type = "beta",
          model = model,
          u = u,
          gaussint_options = gaussint_options,
          method = method,
          control = optim_options)
  beta <- opt_beta$par
  opt_joint <-
    optim(par = c(beta, u),
          fn = fn, gr = gr,
          opt_type = "joint",
          model = model,
          gaussint_options = gaussint_options,
          method = method,
          control = optim_options)
  result <- data.frame(index = seq_len(N_latent),
                       latent = opt_joint$par,
                       iteration = 1,
                       type = "stepwise")
  counts <- data.frame(iteration = 1,
                       fn = opt_beta$counts["function"] +
                         opt_joint$counts["function"],
                       gr_beta = opt_beta$counts["gradient"] +
                         opt_joint$counts["gradient"],
                       gr_u = opt_joint$counts["gradient"])
  list(result = result, counts = counts,
       opt_beta = opt_beta, opt_joint = opt_joint)
}
