% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{summary.mp_estimate}
\alias{summary.mp_estimate}
\alias{mp_estimate_summary}
\alias{summary.mp_model}
\alias{mp_model_summary}
\title{Multiprobit Summary Information}
\usage{
\method{summary}{mp_estimate}(object, ...)

\method{summary}{mp_model}(object, ...)
}
\arguments{
\item{object}{An \code{\link{mp_estimate}} or \code{\link{mp_model}}object}

\item{...}{Additional parameters, currently unused.}
}
\value{
An \code{mp_estimate_summary} or \code{mp_model_summary} list object with
elements
\describe{
\item{beta}{Summary information about the \eqn{B} coefficients
(see \code{\link[=mp_model]{mp_model()}} for model definition).
The elements of the matrix \eqn{B} are listed row-wise, so that the
values for each covariate are grouped together, for all Y-dimensions.}
\item{u}{Summary information about the latent \eqn{u} parameters (see
\code{\link[=mp_model]{mp_model()}} for model definition). This is the internal scale
representation of the probit correlation matrix parameter \eqn{\Sigma}.}
\item{Sigma}{Summary information for \eqn{\Sigma}. Includes a basic linear
error propagation estimate of the elementwise standard deviations, as
\verb{[prior_]sd_linear}; the calculation for the prior in \code{mp_model}
indicates that this is an overestimate of the actual standard deviations
of the \eqn{\Sigma} elements.}
}
For model summaries, the prior mean and standard deviations are provided.
For estimate summaries, the available quantities depend on what was computed
and stored in the input object. Currently, the latent MAP estimates and
their Hessian-based uncertainties are provided, and the corresponding point
estimate of the probit correlation matrix parameter \eqn{\Sigma}.
}
\description{
Organise information about multiprobit parameter estimates
and models
}
\examples{
\dontrun{
if (interactive()) {
  # EXAMPLE1
}
}

}
