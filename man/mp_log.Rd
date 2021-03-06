% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/environment.R
\name{mp_log_reset}
\alias{mp_log_reset}
\alias{mp_log_get}
\alias{mp_log_message}
\title{multiprobit log message methods}
\usage{
mp_log_reset()

mp_log_get()

mp_log_message(
  ...,
  domain = NULL,
  appendLF = TRUE,
  verbosity = 1,
  allow_verbose = TRUE,
  verbose = NULL,
  verbose_store = NULL
)
}
\arguments{
\item{...}{Zero or more objects passed on to \code{\link[base:.makeMessage]{base::.makeMessage()}}}

\item{domain}{Domain for translations, passed on to \code{\link[base:.makeMessage]{base::.makeMessage()}}}

\item{appendLF}{logical; whether to add a newline to the message. Only
used for verbose output.}

\item{verbosity}{numeric value descibing the verbosity level of the message}

\item{allow_verbose}{Whether to allow verbose output. Must be set to FALSE
until the options object has been initialised.}

\item{verbose}{logical, numeric, or \code{NULL}; local override for verbose
output. If \code{NULL}, the global option or default value is used. If \code{FALSE},
no messages are printed. If \code{TRUE}, messages with \code{verbosity} \eqn{\leq 1}
are printed. If numeric, messages with \code{verbosity} \eqn{\leq} \code{verbose} are
printed.}

\item{verbose_store}{Same as \code{verbose}, but controlling what messages are
stored in the global log object.}
}
\value{
\code{mp_log_get} RETURN_VALUE

\code{mp_log_message} OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
\code{mp_log_reset()} clears the log contents.

\code{mp_log_message} DETAILS
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
