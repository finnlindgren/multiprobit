% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform.R
\name{latent_to_nwishart}
\alias{latent_to_nwishart}
\alias{latent_from_nwishart}
\title{Transform latent variables to Normalised Wishart}
\usage{
latent_to_nwishart(x, V_chol, df, lower_chol = FALSE)

latent_from_nwishart(W_chol, V_chol, df, lower_chol = FALSE)
}
\arguments{
\item{x}{A vector of latent variables identifying the Normalised Wishart
matrix, length \eqn{d(d-1)/2}}

\item{V_chol}{The Cholesky factor of the Wishart \eqn{V} parameter}

\item{df}{The Wishart degrees of freedom}

\item{lower_chol}{\code{TRUE} if lower triangular Cholesky factors are used
(default = \code{FALSE})}

\item{W_chol}{The Cholesky factor of the Normalised Wishart matrix}
}
\value{
A list with components
\item{W_chol}{The Cholesky factor of the Normalised Wishart matrix}
\item{s}{Vector of scaling values}
\item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
When \code{lower_chol} is \code{FALSE},
\verb{W_chol = (B_chol \\\%*\\\% V_chol) \\\%*\\\% diag(s)}.

When \code{lower_chol} is \code{TRUE},
\verb{W_chol = s * (V_chol \\\%*\\\% B_chol)}.

A list with components
\item{x}{A vector of latent variables identifying the Normalised Wishart
matrix, length \eqn{d(d-1)/2}}
\item{s}{Vector of scaling values}
\item{B_chol}{The inner Cholesky factor of the Barlett decomposition}
When \code{lower_chol} is \code{FALSE},
\verb{W_chol = (B_chol \\\%*\\\% V_chol) \\\%*\\\% diag(s)}.

When \code{lower_chol} is \code{TRUE},
\verb{W_chol = s * (V_chol \\\%*\\\% B_chol)}.
}
\description{
Transform latent iid \eqn{N(0,1)} to Normalised Wishart
matrices

Transform a Normalised Wishart matrix to latent iid
\eqn{N(0,1)} variables
}
\note{
This is an internal non-exported function

This is an internal non-exported function
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
}
\author{
Finn Lindgren

Finn Lindgren
}
\keyword{internal}
