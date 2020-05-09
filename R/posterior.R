
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param expectation a vector or data frame including variables `value` and
#' `std_err`
#' @param post_scaling An optional multiplicative scaling to be applied after
#' the square root operation, Default: 1
#' @return The input `expectation` object with the `value` and `std_error`
#' elements modified to reflect a square root operation, see Details.
#' @details Matches the mean (the value) and standard deviation (for the
#' standard error) of a positive random variable to those of a log-Normal
#' distribution, and computes the expectation and standard deviation of the
#' square root of that log-Normal distribution.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname sqrt_expectation
#' @keywords internal
#' @export

sqrt_expectation <- function(expectation, post_scaling = 1) {
  expectation_sqrt <- function(mu, sigma) {
    sqrt(mu) * (1 + (sigma / mu)^2)^(-1 / 8)
  }
  std_dev_sqrt <- function(mu, sigma) {
    sqrt(mu * (1 - (1 + (sigma / mu)^2)^(-1 / 4)))
  }
  value <- expectation_sqrt(expectation[["value"]], expectation[["std_err"]])
  std_err <- std_dev_sqrt(expectation[["value"]], expectation[["std_err"]])
  expectation[["value"]] <- value * post_scaling
  expectation[["std_err"]] <- std_err * post_scaling
  expectation
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param estimate PARAM_DESCRIPTION
#' @param N_mc PARAM_DESCRIPTION
#' @param options PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname calc_importance_sampling
#' @export

calc_importance_sampling <- function(estimate, N_mc, options, ...) {
  initial_random_seed <- get(".Random.seed", envir = globalenv())
  set.seed(1L)
  Q_R <- base::chol(-estimate$hessian)
  log_det_Q <- 2 * sum(log(diag(Q_R)))

  N_latent <- nrow(Q_R)

  z_norm <- numeric(N_mc)
  log_imp <- numeric(N_mc)
  log_post <- numeric(N_mc)
  z <- matrix(0, N_mc, N_latent)
  for (loop in seq_len(N_mc)) {
    #  (R^T R)^{-1} = R^{-1} R^{-T}
    #  Cov(R^{-1} z, R^{-1} z) = R^{-1} I R^{-T}
    z[loop, ] <- rnorm(N_latent)
    latent <- estimate$estimate + multiprobit:::tri_solve(Q_R, z[loop, ])

    z_norm[loop] <- sqrt(sum(z[loop, ]^2))
    log_imp[loop] <-
      log_det_Q / 2 - sum(z[loop, ]^2) / 2 - N_latent / 2 * log(2 * pi)
    log_imp[loop] <-
      log_det_Q / 2 - sum(z[loop, ]^2) / 2 - N_latent / 2 * log(2 * pi)
    log_post[loop] <-
      multiprobit:::fn(
        latent, estimate$model, gaussint_options = options$gaussint,
        ..., opt_type = "joint", options = options
      )
  }
  log_r <- log_post - log_imp

  list(
    initial_random_seed = initial_random_seed,
    N_mc = N_mc,
    Q_R = Q_R,
    data = data.frame(
      Index = seq_along(z_norm),
      z_norm = z_norm,
      log_imp = log_imp,
      log_post = log_post,
      log_r = log_r
    ),
    z = z
  )
}

#' @title Pareto Smoothed Importance Sampling Weights
#' @description Wrapper for `loo::psis` with `multiprobit` specific
#' requirements handled.
#' @param log_r logarithm of importance ratios for importance sampling
#' @return Pareto smoothed importance log-ratios. `-Inf` values are converted
#' to `-log(.Machine$double.eps`, and finite values are fed through `loo::psis`,
#' and the output is shifted to match the original offset. `loo::psis` on
#' its own shifts the log-ratios by the maximum input value, but we need the
#' offset to be the same as in the input, so that the resulting
#' \eqn{z^*=\sum_k r_k} doesn't change its interpretation before and after
#' applying the pareto smoothing.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[loo]{weights.importance_sampling}}
#' @rdname log_psis_raw
#' @export
#' @importFrom loo psis

log_psis_raw <- function(log_r) {
  ok <- is.finite(log_r)
  ps <- loo::psis(log_r[ok], r_eff = NA)
  log_r_psis <- numeric(length(log_r))
  log_r_psis[ok] <- ps$log_weights - min(ps$log_weights) + min(log_r[ok])
  log_r_psis[!ok] <- log(.Machine$double.eps)
  log_r_psis
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param log_r PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: "psis"
#' @return `positive`, `log_r`, `log_w`, `r`, `w`
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[matrixStats]{logSumExp}}
#' @rdname calc_imp_weights
#' @export
#' @importFrom matrixStats logSumExp

calc_imp_weights <- function(log_r, method = c("psis", "is", "average")) {
  method <- match.arg(method)
  if (method == "is") {
    ok <- is.finite(log_r)
    log_r <- ifelse(ok, log_r, log(.Machine$double.eps))
  } else if (method == "average") {
    ok <- is.finite(log_r)
    log_r <- rep(matrixStats::logSumExp(log_r[ok]) - log(length(log_r)),
                 length(log_r))
    ok <- rep(TRUE, length(log_r))
  } else {
    ok <- is.finite(log_r)
    log_r = log_psis_raw(log_r)
  }
  r <- exp(log_r)
  r[!ok] <- 0
  log_w <- log_r
  log_w[ok] <- log_r[ok] - matrixStats::logSumExp(log_r[ok])
  w <- exp(log_w)# / sum(exp(log_w))
  w[!ok] <- 0
  data.frame(positive = ok,
             log_r = log_r,
             log_w = log_w,
             r = r,
             w = w)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param imp_weights PARAM_DESCRIPTION
#' @param quantity PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname calc_expectation
#' @export

calc_expectation <- function(imp_weights, quantity) {
  N_q <- length(quantity)
  if (is.vector(imp_weights)) {
    w <- imp_weights
    if ((length(w) == 1) && (N_q > 1)) {
      w <- rep(w, N_q)
    }
  } else {
    w <- imp_weights[["w"]]
  }
  N <- length(w)
  integral <- sum(w * quantity)
  q2 <- (quantity - integral)^2
  std_err <- sum(w^2 * q2)^0.5
  data.frame(
    value = integral,
    std_err = std_err,
    R_eff =  sum(w * q2) / std_err^2 / N
  )
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param iw PARAM_DESCRIPTION
#' @param include_extras Whether to include transformed quantities used to
#' calculate the main quantities, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname calc_diagnostics
#' @export

calc_diagnostics <- function(iw, include_extras = FALSE) {
  N <- nrow(iw)
  N_positive <- sum(iw[["positive"]])
  result <-
    rbind(
      data.frame(
        deviance = "R_eff",
        # sum(w)^2 not needed for normalised weights
        value = sum(iw$w)^2 / sum(iw$w^2) / N,
        std_err = NA_real_,
        R_eff = NA_real_
      ),
      data.frame(
        deviance = "KLD_p_tp",
        calc_expectation(iw, iw$log_w + log(N))
      ),
      data.frame(
        deviance = "KLD_tp_p",
        calc_expectation(1 / N,
                         -iw$log_w[iw$positive] - log(N) -
                           min(iw$log_w[iw$positive]) * (N - N_positive))
      ),
      data.frame(
        deviance = "Hellinger2",
        calc_expectation(1 / N,
                         1 - sqrt(iw$w * N))
      ),
      data.frame(
        deviance = "Hellinger",
        sqrt_expectation(
          calc_expectation(1 / N,
                           1 - sqrt(iw$w * N))
        )
      ),
      data.frame(
        deviance = "JensenShannon",
        calc_expectation(1 / N,
                         N * iw$w / 2 * (iw$log_w + log(2)) +
                           1 / 2 * log(2 / N) -
                           (N * iw$w + 1) / 2 * log(iw$w + 1 / N))
      )
    )
  # Transformations
  result <-
    rbind(
      result,
      data.frame(
        deviance = "JensenShannon_over_log2",
        value = result$value[result$deviance == "JensenShannon"] / log(2),
        std_err = result$std_err[result$deviance == "JensenShannon"] / log(2),
        R_eff = NA_real_
      ),
      data.frame(
        deviance = "Hellinger_sqrt2",
        value = result$value[result$deviance == "Hellinger"] * sqrt(2),
        std_err = result$std_err[result$deviance == "Hellinger"] * sqrt(2),
        R_eff = NA_real_
      ),
      data.frame(
        deviance = "KLD_p_tp_over2_sqrt",
        sqrt_expectation(result[result$deviance == "KLD_p_tp", -1, drop = FALSE],
                         1 / sqrt(2))
      ),
      data.frame(
        deviance = "KLD_tp_p_over2_sqrt",
        sqrt_expectation(result[result$deviance == "KLD_tp_p", -1, drop = FALSE],
                         1 / sqrt(2))
      )
    )
  result <-
    rbind(
      result,
      data.frame(
        deviance = "TV_min",
        value = max(result$value[result$deviance == "Hellinger2"],
                    result$value[result$deviance == "JensenShannon_over_log2"]),
        std_err = max(result$std_err[result$deviance == "Hellinger2"],
                      result$std_err[result$deviance == "JensenShannon_over_log2"]),
        R_eff = NA_real_
      ),
      data.frame(
        deviance = "TV_max",
        value = min(result$value[result$deviance == "Hellinger_sqrt2"],
                    result$value[result$deviance == "KLD_p_tp_over2_sqrt"],
                    result$value[result$deviance == "KLD_tp_p_over2_sqrt"]),
        std_err = max(result$std_err[result$deviance == "Hellinger_sqrt2"],
                      result$std_err[result$deviance == "KLD_p_tp_over2_sqrt"],
                      result$std_err[result$deviance == "KLD_tp_p_over2_sqrt"]),
        R_eff = NA_real_
      )
    )
  main <- c("R_eff", "KLD_p_tp", "KLD_tp_p",
            "Hellinger", "JensenShannon", "TV_min", "TV_max")
  if (include_extras) {
    cbind(result,
          IsExtra = !(result$deviance %in% main))
  } else {
    result[result$deviance %in% main, , drop = FALSE]
  }
}
