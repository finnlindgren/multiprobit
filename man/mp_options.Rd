% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/environment.R
\name{mp_options}
\alias{mp_options}
\title{FUNCTION_TITLE}
\usage{
mp_options(...)
}
\arguments{
\item{...}{PARAM_DESCRIPTION}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
DETAILS
\describe{
\item{gaussint}{Optional list of options for
\code{excursions::gaussint}.
Specific relevant option:
\code{num_threads} The maximum number of allowed threads for parallel
computing by \code{excursions::gaussint}, Default: 0, meaning no limit.
}
\item{optim}{List of control options for \code{optim()}. Can also include
an element \code{method} for choosing the optimisation method.}
\item{strategy}{The estimation optimisation strategy. Options: "alternating", "joint",
"stepwise". Default: "stepwise"}
\item{max_iter}{The maximum number of steps for
\code{strategy == "alternating"}}
}
}
\examples{
\dontrun{
if (interactive()) {
  # EXAMPLE1
}
}
}