% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{mp_logposterior}
\alias{mp_logposterior}
\alias{mp_logposterior_joint}
\alias{mp_logposterior_fixed_beta}
\alias{mp_logposterior_fixed_u}
\title{FUNCTION_TITLE}
\usage{
mp_logposterior(
  Y,
  X,
  Sigma_model,
  prec_beta,
  beta,
  u,
  ...,
  what = c("loglike", "grad", "grad_beta", "grad_u")
)

mp_logposterior_joint(
  latent,
  Y,
  X,
  Sigma_model,
  prec_beta,
  ...,
  what = c("loglike", "grad")
)

mp_logposterior_fixed_beta(
  latent,
  Y,
  X,
  Sigma_model,
  prec_beta,
  beta,
  ...,
  what = c("loglike", "grad")
)

mp_logposterior_fixed_u(
  latent,
  Y,
  X,
  Sigma_model,
  prec_beta,
  u,
  ...,
  what = c("loglike", "grad")
)
}
\arguments{
\item{Y}{PARAM_DESCRIPTION}

\item{X}{PARAM_DESCRIPTION}

\item{Sigma_model}{PARAM_DESCRIPTION}

\item{prec_beta}{PARAM_DESCRIPTION}

\item{beta}{PARAM_DESCRIPTION}

\item{u}{PARAM_DESCRIPTION}

\item{...}{PARAM_DESCRIPTION}

\item{what}{PARAM_DESCRIPTION, Options: "loglike", "grad", with additional
options "grad_beta", "grad_u" for \code{mp_logposterior} only}

\item{latent}{PARAM_DESCRIPTION}
}
\value{
OUTPUT_DESCRIPTION

OUTPUT_DESCRIPTION

OUTPUT_DESCRIPTION

OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
DETAILS

DETAILS

DETAILS

DETAILS
}
\examples{
\dontrun{
if (interactive()) {
  # EXAMPLE1
}
}
\dontrun{
if (interactive()) {
  # EXAMPLE1
}
}
\dontrun{
if (interactive()) {
  # EXAMPLE1
}
}
\dontrun{
if (interactive()) {
  # EXAMPLE1
}
}
}
