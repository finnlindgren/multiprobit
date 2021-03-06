# Loglikelihood, logposterior, and derivatives ####

mp_loglike <- function(Y, mu, Sigma_chol, lower_chol, ...) {
  sum(mpp(Y,
    mu = mu, Sigma_chol = Sigma_chol,
    lower_chol = lower_chol,
    ..., log = TRUE
  )$P)
}
mp_loglike_gradient_beta <- function(Y, X, mu, Sigma_chol, lower_chol,
                                     ...) {
  Q_chol <- inverse_chol_reverse(Sigma_chol, lower_chol = lower_chol)
  d <- ncol(Y)
  perm <- rev(seq_len(d))
  L <- 0
  for (i in seq_len(nrow(Y))) {
    L <- L + kronecker(
      Matrix::t(X[i, , drop = FALSE]),
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
mp_loglike_gradient_u <- function(Y, u, mu, Sigma_model, ...) {
  mpp_gradient_u(Y, u = u, mu = mu, Sigma_model, ..., log = TRUE)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Y PARAM_DESCRIPTION
#' @param X PARAM_DESCRIPTION
#' @param Sigma_model PARAM_DESCRIPTION
#' @param prec_beta PARAM_DESCRIPTION
#' @param beta PARAM_DESCRIPTION
#' @param u PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @param what PARAM_DESCRIPTION, Options: "loglike", "grad", with additional
#'   options "grad_beta", "grad_u" for `mp_logposterior` only
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

mp_logposterior <-
  function(Y, X, Sigma_model, prec_beta,
           beta, u, ...,
           what = c("loglike", "grad", "grad_beta", "grad_u")) {
    what <- match.arg(what)
    d <- ncol(Y)
    dd <- d * (d - 1) / 2
    J <- ncol(X)

    stopifnot(length(beta) + length(u) == (J * d + dd))

    mu <- X %*% rvec(beta, d)
    if (what != "grad_u") {
      Sigma_chol <- wm_chol(Sigma_model, latent = u)
    }
    if (what == "loglike") {
      L <- -sum(u)^2 / 2 - prec_beta / 2 * sum(beta^2) +
        mp_loglike(
          Y = Y, mu = mu, Sigma_chol = Sigma_chol,
          lower_chol = Sigma_model$lower_chol,
          ...
        )
    } else {
      if (what %in% c("grad", "grad_beta")) {
        dL_dbeta <- -prec_beta * beta +
          mp_loglike_gradient_beta(
            Y = Y, X = X, mu = mu,
            Sigma_chol = Sigma_chol,
            lower_chol = Sigma_model$lower_chol,
            ...
          )
      }
      if (what %in% c("grad", "grad_u")) {
        dL_du <- -u + mp_loglike_gradient_u(
          Y = Y, u = u, mu = mu,
          Sigma_model = Sigma_model,
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

mp_logposterior_joint <- function(latent, Y, X, Sigma_model, prec_beta, ...,
                                  what = c("loglike", "grad")) {
  what <- match.arg(what)

  d <- ncol(Y)
  dd <- d * (d - 1) / 2
  J <- ncol(X)

  stopifnot(length(latent) == (J * d + dd))

  index_beta <- seq_len(J * d)
  index_u <- J * d + seq_len(dd)

  mp_logposterior(Y, X, Sigma_model, prec_beta,
    beta = latent[index_beta],
    u = latent[index_u],
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

mp_logposterior_fixed_beta <- function(latent, Y, X, Sigma_model, prec_beta,
                                       beta, ...,
                                       what = c("loglike", "grad")) {
  what <- match.arg(what)

  d <- ncol(Y)
  dd <- d * (d - 1) / 2

  stopifnot(length(latent) == (dd))

  mp_logposterior(Y, X, Sigma_model, prec_beta,
    beta = beta,
    u = latent,
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

mp_logposterior_fixed_u <- function(latent, Y, X, Sigma_model, prec_beta,
                                    u, ...,
                                    what = c("loglike", "grad")) {
  what <- match.arg(what)

  d <- ncol(Y)
  J <- ncol(X)

  stopifnot(length(latent) == (J * d))

  mp_logposterior(Y, X, Sigma_model, prec_beta,
    beta = latent,
    u = u,
    ...,
    what = ifelse(what == "loglike", "loglike", "grad_beta")
  )
}





# User interface ####


#' @title Setup multivariate probit model
#' @description Construct a multivariate probit model object that stores
#' information about model structure and parameters
#' @param model An optional existing model object to be updated.
#' @param response A matrix with n-by-d elements, where each row is a
#'   multivariate observation, see Details. A vector is interpreted as a
#'   single row matrix.
#' @param X An optimally precomputed n-by-J model matrix, where J is the number
#'   regression coeficcients for each of the d dimensions.
#' @param formula A formula interpretable by `model.matrix`.
#' @param data A `data.frame` containing the veriables needed by the
#'   formula.
#' @param df Degrees of freedom for the normalised Wishart prior for the
#'   correlation matrix. See Details.
#' @param prec_beta Prior precision for the regression coefficients
#' @return An object of class `mp_model`
#' @details The multivariate probit model has a multivariate binary
#' response variable, here denoted \eqn{Y}. The model is built from a
#' linear predictor
#' \deqn{M = X B}
#' where \eqn{X} is a n-by-J matrix of \eqn{J} predictors, and \eqn{B} is
#' a J-by-d matrix of regression coefficients.
#' Each row of \eqn{M} is the linear
#' predictor for one multivariate observation. The response variables \eqn{Y}
#' are linked to \eqn{M} by first defining latent Gaussian variables
#' \deqn{Z=M+E} where each row of \eqn{E} is a multivariate Normal vector,
#' \eqn{E \sim N(0,\Sigma)}. Then, \deqn{Y_{i,k}=I(Z_{i,k} > 0).}
#' Conditionally on \eqn{B}, each row of \eqn{Y} has a multinomial distribution
#' on the set of all \eqn{0/1} combinations, with each probability equal to a
#' hyperquadrant probability of a the multivariate Normal distribution
#' \eqn{N(\mu,\Sigma)}, where \eqn{\mu} is the corresponding row of \eqn{M}.
#'
#' Only the inequality \eqn{Y_{i,k} > 0} for the response variables is used,
#' so alternative data representations such as \eqn{-1/+1} will also work as
#' expected.
#'
#' The degrees of freedom for the normalised Wishart prior are linked to the
#' concentration pararameter \eqn{\eta} of the LKJ prior by the relation
#' `df = 2 * eta + d - 1`, which makes the two models equivalent.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @importFrom stats model.matrix
#' @rdname mp_model

mp_model <- function(model = NULL,
                     response = NULL,
                     X = NULL,
                     formula = NULL,
                     data = NULL,
                     df = NULL,
                     prec_beta = NULL) {
  if (!is.null(model)) {
    stopifnot(inherits(model, "mp_model"))
  } else {
    model <- list(
      Y = NULL,
      X = NULL,
      formula = NULL,
      data = NULL,
      Sigma_model = NULL,
      prec_beta = NULL
    )
  }
  if (!is.null(df)) {
    # Placeholder until model dimension is known
    model$Sigma_model <-
      wm_model(type = "nwishart", V_chol = 1, df = df)
  }
  if (!is.null(response)) {
    model$Y <- as.matrix(response) > 0
  }

  # New X provided, remove formula & data, formula|data provided is error
  # New formula provided, compute new X if data present, otherwise clear X
  # New data provided, compute new X if formula present, otherwise clear X
  # X not existing, and none provided, do nothing
  if (!is.null(X)) {
    model$X <- as.matrix(X)
    model$formula <- NULL
    model$data <- NULL
    if (!is.null(formula) || !is.null(data)) {
      stop("When X is provided, do not also provide formula & data")
    }
  } else if (!is.null(formula) || !is.null(data)) {
    if (!is.null(formula)) {
      model$formula <- formula
    }
    if (!is.null(data)) {
      model$data <- data
    }
    if (is.null(model$formula) || is.null(model$data)) {
      model$X <- NULL
    } else {
      model$X <- model.matrix(model$formula, data = model$data)
    }
  } else {
    # Do nothing.
  }
  # Extract covariate names
  if (is.null(model$X)) {
    model$x_names <- NULL
  } else {
    model$x_names <- colnames(model$X)
    if (is.null(model$x_names)) {
      model$x_names <- paste0(
        mp_options_get(
          "x_name_prefix",
          include_default = TRUE
        ),
        seq_len(ncol(model$X))
      )
      colnames(model$X) <- model$x_names
    }
  }

  if (!is.null(prec_beta)) {
    model$prec_beta <- prec_beta
  }
  if (is.null(model$prec_beta)) {
    model$prec_beta <- 0.0
  }

  if (!is.null(model$Y)) {
    stopifnot(nrow(model$X) == nrow(model$Y))

    if (is.null(model$Sigma_model)) {
      warning("Model degrees of freedom parameter 'df' is not set.")
    } else {
      model$Sigma_model <-
        wm_model(
          type = "nwishart",
          V_chol = sparse_identity(ncol(model$Y)),
          df = model$Sigma_model$df
        )
    }
  }

  # Extract covariate names
  if (is.null(model$Y)) {
    model$y_names <- NULL
  } else {
    model$y_names <- colnames(model$Y)
    if (is.null(model$y_names)) {
      model$y_names <- paste0(
        mp_options_get("y_name_prefix",
          include_default = TRUE
        ),
        seq_len(ncol(model$Y))
      )
      colnames(model$Y) <- model$y_names
    }
  }

  class(model) <- "mp_model"
  model
}




#' @aliases multiprobit mp_estimate
#' @title Estimate a multivariate probit model
#' @description Estimate a multivariate probit model from multivariate binary
#'   data in a Bayesian generalised linear model framework
#' @param model An [`mp_model`] model object.
#' @param ... Parameters passed on to [mp_model()]
#' @param options An [`mp_options`] object or a list that can be coerced into
#' an `mp_options` object. Options set here will override the global options.
#' @return An `mp_estimate` object.
#' @details For details on the multivariate probit model, see
#' [mp_model()]. The `multiprobit` function estimates
#' the model for observations stored in `response`
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   N <- 6
#'   d <- 2
#'   J <- 1
#'
#'   set.seed(1L)
#'   X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
#'   B <- matrix(0.5, J, d)
#'   Y <- matrix(rnorm(N * d, mean = as.vector(X %*% B)) > 0, N, d)
#'   df <- d + 1
#'   prec_beta <- 0.1
#'
#'   model <- mp_model(
#'     response = Y, X = X,
#'     df = df, prec_beta = prec_beta
#'   )
#'   opt <- multiprobit(
#'     model = model,
#'     options =
#'       mp_options(
#'         gaussint = list(max.threads = 1),
#'         strategy = "stepwise"
#'       )
#'   )
#' }
#' }
#' @export
#' @importFrom stats model.matrix
#' @importFrom stats optim
#' @importFrom stats optimHess
#' @rdname multiprobit

multiprobit <- function(model = NULL,
                        ...,
                        options = NULL) {
  model <- mp_model(model, ...)
  if (is.null(model$Y)) {
    stop(
      paste0(
        "No response/observations provided:\n",
        "Add to an existing model with mp_model(model, response = ...)."
      )
    )
  }
  options <- mp_options(mp_options_default(), mp_options_get(), options)
  mp_options_check(options)

  mp_log_message("Starting estimation with strategy '", options$strategy, "'",
                 verbosity = 1,
                 verbose = options[["verbose"]])

  if (options$strategy == "alternating") {
    opt <- optim_alternating(model, options = options)
  } else if (options$strategy == "joint") {
    opt <- optim_joint(model, options = options)
  } else if (options$strategy == "stepwise") {
    opt <- optim_stepwise(model, options = options)
  } else {
    stop("Unknown strategy; This cannot happen.")
  }

  mp_log_message("Finished estimation with strategy '", options$strategy, "'",
                 verbosity = 1,
                 verbose = options[["verbose"]])

  if (!is.null(options[["hessian"]])) {
    mp_log_message("Calculating hessian",
                   verbosity = 1,
                   verbose = options[["verbose"]])
    hessian <- calc_hessian(opt$estimate, model, options)
  } else {
    hessian <- NULL
  }

  est <- c(list(model = model, options = options, hessian = hessian), opt)
  class(est) <- "mp_estimate"
  est
}



# Optimisation methods ####


fn <- function(x, model, ..., opt_type = c("joint", "beta", "u"),
               options = NULL) {
  opt_type <- match.arg(opt_type)
  mp_log_message(
    paste0("f,latent(", opt_type, ") = (", paste0(x, collapse = ", "), ")"),
    verbosity = 3,
    verbose = options[["verbose"]]
  )
  what <- "loglike"
  f <- tryCatch(
    switch(opt_type,
      "joint" = do.call(
        mp_logposterior_joint,
        c(
          list(latent = x, what = what, ...),
          model
        )
      ),
      "beta" = do.call(
        mp_logposterior_fixed_u,
        c(
          list(latent = x, what = what, ...),
          model
        )
      ),
      "u" = do.call(
        mp_logposterior_fixed_beta,
        c(
          list(latent = x, what = what, ...),
          model
        )
      )
    ),
    error = function(e) {
      # Return -Inf for invalid parameter combinations
      -Inf
    }
  )
  mp_log_message(
    paste0("f(", opt_type, ") = ", paste0(f, collapse = ", "), ""),
    verbosity = 3,
    verbose = options[["verbose"]]
  )
  f
}
gr <- function(x, model, ..., opt_type = c("joint", "beta", "u"),
               options = NULL) {
  opt_type <- match.arg(opt_type)
  mp_log_message(
    paste0("g,latent(", opt_type, ") = (", paste0(x, collapse = ", "), ")"),
    verbosity = 3,
    verbose = options[["verbose"]]
  )
  what <- "grad"
  g <- switch(opt_type,
    "joint" = do.call(
      mp_logposterior_joint,
      c(
        list(latent = x, what = what, ...),
        model
      )
    ),
    "beta" = do.call(
      mp_logposterior_fixed_u,
      c(
        list(latent = x, what = what, ...),
        model
      )
    ),
    "u" = do.call(
      mp_logposterior_fixed_beta,
      c(
        list(latent = x, what = what, ...),
        model
      )
    )
  )
  mp_log_message(
    paste0("g(", opt_type, ") = (", paste0(g, collapse = ", "), ")"),
    verbosity = 3,
    verbose = options[["verbose"]]
  )
  g
}


calc_hessian <- function(latent, model, options) {
  if (is.null(options$hessian) ||
    (options$hessian == "none")) {
    hessian <- NULL
  } else {
    if (options$hessian == "full") {
      hessian <-
        do.call(
          optimHess,
          c(
            list(
              par = latent,
              fn = fn, gr = gr,
              opt_type = "joint",
              model = model,
              gaussint_options = options$gaussint,
              options = options
            ),
            options$optim
          )
        )
    } else if (options$hessian == "block") {
      d <- ncol(model$Y)
      N_beta <- ncol(model$X) * d
      N_u <- model$Sigma_model$N_latent
      hessian <-
        Matrix::bdiag(
          do.call(
            optimHess,
            c(
              list(
                par = latent[seq_len(N_beta)],
                fn = fn, gr = gr,
                opt_type = "beta",
                u = latent[N_beta + seq_len(N_u)],
                model = model,
                gaussint_options = options$gaussint,
                options = options
              ),
              options$optim
            )
          ),
          do.call(
            optimHess,
            c(
              list(
                par = latent[N_beta + seq_len(N_u)],
                fn = fn, gr = gr,
                opt_type = "u",
                beta = latent[seq_len(N_beta)],
                model = model,
                gaussint_options = options$gaussint,
                options = options
              ),
              options$optim
            )
          )
        )
    } else if (options$hessian == "diagonal") {
      warning("'diagonal' hessian not implemented.")
      hessian <- NULL
    } else {
      warning("Unknown hessian type '", options$hessian, "'")
    }
  }
  hessian
}

optim_alternating <- function(model, options = NULL) {
  options <- mp_options(
    mp_options_default(),
    # Local defaults:
    mp_options(gaussint = list(seed = 1L)),
    # Caller options:
    options
  )

  N <- nrow(model$Y)
  d <- ncol(model$Y)
  J <- ncol(model$X)
  N_beta <- J * d
  N_u <- model$Sigma_model$N_latent
  N_latent <- N_beta + N_u
  opt_beta <- list()
  opt_u <- list()
  result <-
    data.frame(
      index = rep(seq_len(N_latent),
        times = options$max_iter
      ),
      latent = 0,
      iteration = rep(seq_len(options$max_iter), each = N_latent),
      type = "alternating"
    )
  counts <-
    data.frame(
      iteration = seq_len(options$max_iter),
      fn = numeric(options$max_iter),
      gr_beta = numeric(options$max_iter),
      gr_u = numeric(options$max_iter)
    )
  beta <- rep(0, N_beta)
  u <- rep(0, N_u)
  for (loop in seq_len(options$max_iter)) {
    if (loop > 1) {
      counts[loop, c("fn", "gr_beta", "gr_u")] <-
        counts[loop - 1, c("fn", "gr_beta", "gr_u")]
    }
    opt_beta[[loop]] <-
      do.call(
        optim,
        c(
          list(
            par = beta,
            fn = fn, gr = gr,
            opt_type = "beta",
            model = model,
            u = u,
            gaussint_options = options$gaussint,
            options = options
          ),
          options$optim
        )
      )
    beta <- opt_beta[[loop]]$par
    counts$fn[loop] <-
      counts$fn[loop] + opt_beta[[loop]]$counts["function"]
    counts$gr_beta[loop] <-
      counts$gr_beta[loop] + opt_beta[[loop]]$counts["gradient"]
    opt_u[[loop]] <-
      do.call(
        optim,
        c(
          list(
            par = u,
            fn = fn, gr = gr,
            opt_type = "u",
            model = model,
            beta = beta,
            gaussint_options = options$gaussint,
            options = options
          ),
          options$optim
        )
      )
    u <- opt_u[[loop]]$par
    counts$fn[loop] <- counts$fn[loop] + opt_u[[loop]]$counts["function"]
    counts$gr_u[loop] <- counts$gr_u[loop] + opt_u[[loop]]$counts["gradient"]

    result$latent[(loop - 1) * N_latent + seq_len(N_latent)] <- c(beta, u)
  }
  loop <- options$max_iter
  list(
    estimate = result$latent[(options$max_iter - 1) * N_latent + seq_len(N_latent)],
    result = result, counts = counts, opt_beta = opt_beta, opt_u = opt_u
  )
}

optim_joint <- function(model, options = NULL) {
  options <- mp_options(
    mp_options_default(),
    # Local defaults:
    mp_options(gaussint = list(seed = 1L)),
    # Caller options:
    options
  )

  N <- nrow(model$Y)
  d <- ncol(model$Y)
  J <- ncol(model$X)
  N_beta <- J * d
  N_u <- model$Sigma_model$N_latent
  N_latent <- N_beta + N_u
  beta <- rep(0, N_beta)
  u <- rep(0, N_u)
  opt_joint <-
    do.call(
      optim,
      c(
        list(
          par = c(beta, u),
          fn = fn, gr = gr,
          opt_type = "joint",
          model = model,
          gaussint_options = options$gaussint,
          options = options
        ),
        options$optim
      )
    )
  result <- data.frame(
    index = seq_len(N_latent),
    latent = opt_joint$par,
    iteration = 1,
    type = "joint"
  )
  counts <- data.frame(
    iteration = 1,
    fn = opt_joint$counts["function"],
    gr_beta = opt_joint$counts["gradient"],
    gr_u = opt_joint$counts["gradient"]
  )
  list(
    estimate = result$latent,
    result = result, counts = counts, opt_joint = opt_joint
  )
}


optim_stepwise <- function(model, options = NULL) {
  options <- mp_options(
    mp_options_default(),
    # Local defaults:
    mp_options(gaussint = list(seed = 1L)),
    # Caller options:
    options
  )

  N <- nrow(model$Y)
  d <- ncol(model$Y)
  J <- ncol(model$X)
  N_beta <- J * d
  N_u <- model$Sigma_model$N_latent
  N_latent <- N_beta + N_u
  beta <- rep(0, N_beta)
  u <- rep(0, N_u)
  opt_beta <-
    do.call(
      optim,
      c(
        list(
          par = beta,
          fn = fn, gr = gr,
          opt_type = "beta",
          model = model,
          u = u,
          gaussint_options = options$gaussint,
          options = options
        ),
        options$optim
      )
    )
  beta <- opt_beta$par
  opt_joint <-
    do.call(
      optim,
      c(
        list(
          par = c(beta, u),
          fn = fn, gr = gr,
          opt_type = "joint",
          model = model,
          gaussint_options = options$gaussint,
          options = options
        ),
        options$optim
      )
    )
  result <- data.frame(
    index = seq_len(N_latent),
    latent = opt_joint$par,
    iteration = 1,
    type = "stepwise"
  )
  counts <- data.frame(
    iteration = 1,
    fn = opt_beta$counts["function"] +
      opt_joint$counts["function"],
    gr_beta = opt_beta$counts["gradient"] +
      opt_joint$counts["gradient"],
    gr_u = opt_joint$counts["gradient"]
  )
  list(
    estimate = result$latent,
    result = result, counts = counts,
    opt_beta = opt_beta, opt_joint = opt_joint
  )
}



# Summary output ####


#' @title Multiprobit Summary Information
#' @description Organise information about multiprobit parameter estimates
#' and models
#' @param object An [`mp_estimate`] or [`mp_model`]object
#' @param ... Additional parameters, currently unused.
#' @return An `mp_estimate_summary` or `mp_model_summary` list object with
#' elements
#' \describe{
#' \item{beta}{Summary information about the \eqn{B} coefficients
#' (see [mp_model()] for model definition).
#' The elements of the matrix \eqn{B} are listed row-wise, so that the
#' values for each covariate are grouped together, for all Y-dimensions.}
#' \item{u}{Summary information about the latent \eqn{u} parameters (see
#' [mp_model()] for model definition). This is the internal scale
#' representation of the probit correlation matrix parameter \eqn{\Sigma}.}
#' \item{Sigma}{Summary information for \eqn{\Sigma}. Includes a basic linear
#' error propagation estimate of the elementwise standard deviations, as
#' `[prior_]sd_linear`; the calculation for the prior in `mp_model`
#' indicates that this is an overestimate of the actual standard deviations
#' of the \eqn{\Sigma} elements.}
#' }
#' For model summaries, the prior mean and standard deviations are provided.
#' For estimate summaries, the available quantities depend on what was computed
#' and stored in the input object. Currently, the latent MAP estimates and
#' their Hessian-based uncertainties are provided, and the corresponding point
#' estimate of the probit correlation matrix parameter \eqn{\Sigma}.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#'
#' @aliases mp_estimate_summary
#' @method summary mp_estimate
#' @export
#' @rdname summary.mp_estimate

summary.mp_estimate <- function(object, ...) {
  out <- list()
  J <- ncol(object$model$X)
  d <- ncol(object$model$Y)
  N_beta <- J * d
  N_u <- object$model$Sigma_model$N_latent
  S <- solve(-object$hessian)
  out$beta <- data.frame(
    x_name = rep(object$model$x_names, each = d),
    y_name = rep(object$model$y_names, times = J),
    estimate = object$estimate[seq_len(N_beta)],
    sd = Matrix::diag(S)[seq_len(N_beta)]^0.5
  )
  out$u <- data.frame(
    name = paste0("u", seq_len(N_u)),
    estimate = object$estimate[N_beta + seq_len(N_u)],
    sd = Matrix::diag(S)[N_beta + seq_len(N_u)]^0.5
  )

  Sigma <-
    wm_matrix(object$model$Sigma_model, latent = out$u$estimate)
  Sigma_moments_linear <-
    wm_moments_linear(
      object$model$Sigma_model,
      mean_latent = out$u$estimate,
      cov_latent = S[N_beta + seq_len(N_u), N_beta + seq_len(N_u), drop = FALSE]
    )
  out$Sigma <- list(
    estimate = Sigma,
    mean_linear = Sigma_moments_linear$mean,
    sd_linear = Sigma_moments_linear$sd
  )

  class(out) <- "mp_estimate_summary"
  out
}



#' @aliases mp_model_summary
#' @method summary mp_model
#' @export
#' @rdname summary.mp_estimate

summary.mp_model <- function(object, ...) {
  out <- list()
  J <- ncol(object$X)
  d <- ncol(object$Y)
  N_beta <- J * d
  N_u <- object$Sigma_model$N_latent
  out$beta <- data.frame(
    x_name = rep(object$x_names, each = d),
    y_name = rep(object$y_names, times = J),
    prior_mean = rep(0, N_beta),
    prior_sd = rep(object$prec_beta^(-0.5), N_beta)
  )
  out$u <- data.frame(
    name = paste0("u", seq_len(N_u)),
    prior_mean = rep(0, N_u),
    prior_sd = rep(1, N_u)
  )

  Sigma_mean <- wm_matrix(object$Sigma_model, latent = out$u$prior_mean)
  Sigma_sd <- matrix(1 / object$Sigma_model$df^0.5, d, d)
  Matrix::diag(Sigma_sd) <- 0
  Sigma_moments_linear <-
    wm_moments_linear(
      object$Sigma_model,
      mean_latent = out$u$prior_mean,
      cov_latent = sparse_identity(object$Sigma_model$N_latent),
    )
  out$Sigma <- list(
    prior_mean = Sigma_mean,
    prior_sd = Sigma_sd,
    prior_mean_linear = Sigma_moments_linear$mean,
    prior_sd_linear = Sigma_moments_linear$sd
  )

  class(out) <- "mp_model_summary"
  out
}



#' @title Print multiprobit objects
#' @param x an [`mp_model`] or [`mp_estimate_summary`] object to be printed.
#' @param ... further arguments passed to or from other methods.
#' @method print mp_estimate_summary
#' @export
#' @rdname print_mp_objects

print.mp_estimate_summary <- function(x, ...) {
  cat("beta:\n")
  print(x$beta)
  cat("u:\n")
  print(x$u)
  cat("Sigma:\n")
  print(x$Sigma)
}


#' @method print mp_model_summary
#' @export
#' @rdname print_mp_objects

print.mp_model_summary <- function(x, ...) {
  cat("beta:\n")
  print(x$beta)
  cat("u:\n")
  print(x$u)
  cat("Sigma:\n")
  print(x$Sigma)
}
