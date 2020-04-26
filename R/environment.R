#' @include aaaaa.R

#' @title Get access to the internal environment
#' @details The environment is defined in aaaaa.R which is loaded first.
#' @keywords internal
multiprobit_env_get <- function() {
  pkg_envir <- parent.env(environment())
  envir <- get0(".multiprobit_envir", envir = pkg_envir)
  if (!is.environment(envir)) {
    stop("Something went wrong: cannot find internal multiprobit environment.")
  }
  envir
}


# Documentation would clash with base .onLoad documentation
# @title Initialise log storage and global options
# @param libname a character string giving the library directory where the
#   package defining the namespace was found.
# @param pkgname a character string giving the name of the package.
# @aliases namespace_hooks
# @keywords internal
# @rdname namespace_hooks
.onLoad <- function(libname, pkgname) {
  mp_log_reset()
  mp_log_message("multiprobit attached", allow_verbose = FALSE)
  mp_log_message("Set default options", allow_verbose = FALSE)
  mp_options_reset()
}






#' @title multiprobit log message methods
#' @description FUNCTION_DESCRIPTION
#' @details `mp_log_reset()` clears the log contents.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_log

mp_log_reset <- function() {
  envir <- multiprobit_env_get()
  envir$log <- character(0)
  invisible()
}


#' @return `mp_log_get` RETURN_VALUE
#' @export
#' @rdname mp_log

mp_log_get <- function() {
  multiprobit_env_get()$log
}

#' @param ... Zero or more objects passed on to [`base::.makeMessage()`]
#' @param domain Domain for translations, passed on to [`base::.makeMessage()`]
#' @param appendLF logical; whether to add a newline to the message. Only
#'   used for verbose output.
#' @param allow_verbose Whether to allow verbose output. Must be set to FALSE
#' until the options object has been initialised.
#' @return `mp_log_message` OUTPUT_DESCRIPTION
#' @details `mp_log_message` DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_log

mp_log_message <- function(..., domain = NULL, appendLF = TRUE,
                           allow_verbose = TRUE) {
  if (allow_verbose) {
    if (mp_options_get("verbose")) {
      message(..., domain = domain, appendLF = appendLF)
    }
  }
  envir <- multiprobit_env_get()
  envir$log <- c(
    envir$log,
    .makeMessage(Sys.time(), ": ", ...,
      domain = domain,
      appendLF = FALSE
    )
  )
  invisible()
}





#' @title Create or update an options objects
#' @description Create a new options object, or merge information from several
#'   objects.
#'
#'   The `_get`, `_set`, and `_reset` functions operate on a global
#'   package options object. In many cases, setting options in specific calls
#'   to [multiprobit()] is recommended instead.
#' @param ... Named options, optionally started by one or more [`mp_options`]
#' objects. Options specified later override the previous options.
#' @return `mp_options()` returns an `mp_options` object.
#' @details For `mp_options` and `mp_options_set`, recognised options are:
#' \describe{
#' \item{verbose}{logical; if `TRUE`, log messages are printed with `message()`}
#' \item{gaussint}{List of options for
#'   `excursions::gaussint`.
#'   Specific relevant options:
#'   \describe{
#'     \item{num.threads}{The maximum number of allowed threads for parallel
#'       computing by `excursions::gaussint`, Default: 0, meaning no limit.}
#'     \item{seed}{The seed for the internal random number generator for
#'       `excursions::gaussint`.}
#'   }
#' }
#' \item{optim}{List of options for `optim()`. Relevant options are `method`
#'   for choosing the optimisation method, and `control` for setting
#'   convergence criteria and verbose output.}
#' \item{strategy}{The estimation optimisation strategy.
#'   Options: "alternating", "joint", "stepwise". Default: "stepwise"}
#' \item{max_iter}{The maximum number of steps for
#'   `strategy == "alternating"`}
#' \item{hessian}{The hessian computation style. Options: "none", "diagonal",
#'   "block", and "full". Default: "full"}
#'   }
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Compine global and user options:
#'   options1 <- mp_options(mp_options_get(), hessian = "block")
#'   # Create a proto-options object in two equivalent ways:
#'   options2 <- as.mp_options(hessian = "diagonal")
#'   options2 <- as.mp_options(list(hessian = "diagonal"))
#'   # Combine options objects:
#'   options3 <- mp_options(options1, options2)
#' }
#' }
#' @export
#' @rdname mp_options

mp_options <- function(...) {
  new_mp_options <- function() {
    options <- list()
    class(options) <- c("mp_options", "list")
    options
  }
  traverse <- function(options, override) {
    stopifnot(is.list(options) && is.list(override))
    for (x in names(override)) {
      if (inherits(options[[x]], "list") &&
        inherits(override[[x]], "list")) {
        options[[x]] <- traverse(
          options[[x]],
          override[[x]]
        )
      } else {
        options[[x]] <- override[[x]]
      }
    }
    options
  }

  input_options <- list(...)
  n_input <- length(input_options)
  if (n_input == 0) {
    # Return an empty options object
    return(new_mp_options())
  }
  if (inherits(input_options[[1]], "mp_options")) {
    options <- input_options[[1]]
    k <- 1
    while ((k < n_input) &&
      inherits(input_options[[k + 1]], "mp_options")) {
      k <- k + 1
      options <- traverse(options, input_options[[k]])
    }
    if (k < n_input) {
      if (any(vapply(
        input_options[-seq_len(k)],
        function(x) inherits(x, "mp_options"),
        TRUE
      ))
      ) {
        stop("All mp_options objects must come first in the parameter list.")
      }
      options <- traverse(options, input_options[-seq_len(k)])
    }
  } else {
    options <- traverse(new_mp_options(), input_options)
  }

  options
}

#' @param x An object to be converted to an `mp_options` object.
#' @return For `as.mp_options()`, `NULL` or no input returns an empty `mp_options`
#' object, a `list` is converted via `mp_options(...)`, and `mp_options` input
#' is passed through. Other types of input generates an error.
#'
#' @export
#' @rdname mp_options

as.mp_options <- function(x = NULL) {
  if (inherits(x, "mp_options")) {
    x
  } else if (is.null(x)) {
    mp_options()
  } else if (is.list(x)) {
    do.call(mp_options, x)
  } else {
    stop("Cannot coerce object to 'mp_options'")
  }
}

#' @return `mp_options_default()` returns an `mp_options` object containing
#'   default options.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_options

mp_options_default <- function() {
  mp_options(
    verbose = FALSE,
    optim = list(
      method = "BFGS",
      control = list(
        fnscale = -1
      )
    ),
    gaussint = list(),
    max_iter = 5,
    strategy = "stepwise",
    hessian = "full"
  )
}


#' @details `mp_options_check` checks for valid contents of an `mp_options`` object
#' @param options An `mp_options` object to be checked
#' @return `mp_options_check()` returns a `logical`; `TRUE` if the object
#'   contains valid options for use by other functions
#' @details `mp_options_check()` produces warnings for invalid options.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   mp_options_check(mp_options(invalid = "something"))
#' }
#' }
#' @export
#' @rdname mp_options

mp_options_check <- function(options) {
  options <- as.mp_options(options)
  ok <- TRUE
  # Check valid max_iter
  if (is.null(options[["max_iter"]])) {
    ok <- FALSE
    warning("'max_iter' should be a positive integer, not NULL.")
  } else if (!is.numeric(options[["max_iter"]]) ||
    !(options[["max_iter"]] > 0)) {
    ok <- FALSE
    warning("'max_iter' should be a positive integer.")
  }

  # Check valid strategy
  if (is.null(options[["strategy"]])) {
    ok <- FALSE
    warning("'strategy' should not be NULL.")
  } else if (!is.character(options[["strategy"]]) ||
    !(length(options[["strategy"]]) == 1) ||
    !(options[["strategy"]] %in%
      c(
        "alternating",
        "stepwise",
        "joint"
      ))) {
    ok <- FALSE
    warning("'strategy' should be a valid strategy string, see ?mp_options")
  }


  # Check valid hessian setting
  if (is.null(options[["hessian"]])) {
    ok <- FALSE
    warning("'hessian' should not be NULL.")
  } else if (!is.character(options[["hessian"]]) ||
    !(length(options[["hessian"]]) == 1) ||
    !(options[["hessian"]] %in%
      c(
        "none",
        "diagonal",
        "block",
        "full"
      ))) {
    ok <- FALSE
    warning("'hessian' should be a valid hessian string, see ?mp_options")
  }

  ok
}



#' @param name Either `NULL`, or single option name string, or a list of
#'   option namescharacter vector or list with option names,
#'   Default: NULL
#' @return `mp_options_get` returns either an [`mp_options`] object, for
#'   `name == NULL`, the contents of single option, if `name` is a options name
#'   string, or a named list of option contents, if `name` is a list of option
#'   name strings.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_options

mp_options_get <- function(name = NULL) {
  if (is.null(name)) {
    return(multiprobit_env_get()$options)
  }
  if (is.list(name)) {
    mget(unlist(name), as.environment(multiprobit_env_get()$options))
  } else {
    multiprobit_env_get()$options[[name]]
  }
}


#' @details `mp_options_set()` is used to set global package options.
#' @return `mp_options_set()` returns a copy of the global options, invisibly.
#' @seealso [mp_options()], [mp_options_default()], [mp_options_get()]
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   mp_options_set(
#'     gauusint = list(max.threads = 1),
#'     optim = list(control = list(trace = 5))
#'   )
#' }
#' }
#' @export
#' @rdname mp_options

mp_options_set <- function(...) {
  envir <- multiprobit_env_get()
  envir$options <- mp_options(envir$options, ...)
  invisible(mp_options_get())
}

#' @details After its first run, `mp_options_reset()` is equivalent to
#'   `mp_options_set(mp_options_default())`.
#' @export
#' @rdname mp_options

mp_options_reset <- function() {
  envir <- multiprobit_env_get()
  envir$options <- mp_options_default()
  invisible(mp_options_get())
}
