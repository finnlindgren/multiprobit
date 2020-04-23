#' @include aaaaa.R

#' @title Get access to the internal environment
#' @details The environment is defined in aaaaa.R which is loaded first.
#' @keywords internal
multiprobit_env_get <- function() {
  if (exists(".multiprobit_envir") && is.environment(.multiprobit_envir)) {
    return(.multiprobit_envir)
  } else {
    stop("Something went wrong: cannot find internal multiprobit environment.")
  }
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
  .multiprobit_envir$log <- character(0)
  mp_log_message("multiprobit attached", allow_verbose = FALSE)
  mp_log_message("Set default options", allow_verbose = FALSE)
  mp_options_set(default = TRUE)
}



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_log_get

mp_log_get <- function() {
  multiprobit_env_get()$log
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ... zero or more objects passed on to [`base::.makeMessage()`]
#' @param domain Domain for translations, passed on to [`base::.makeMessage()`]
#' @param appendLF logical; whether to add a newline to the message. Only
#'   used for verbose output.
#' @param allow_verbose Whether to allow verbose output. Must be set to FALSE
#' until the options object has been initialised.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_log_message

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
                 appendLF = FALSE)
  )
  invisible()
}




#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_options_get

mp_options_get <- function(x = NULL) {
  if (is.null(x)) {
    return(multiprobit_env_get()$options)
  }
  if (is.list(x)) {
    mget(unlist(x), as.environment(multiprobit_env_get()$options))
  } else {
    multiprobit_env_get()$options[[x]]
  }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' \describe{
#' \item{gaussint}{Optional list of options for
#'   `excursions::gaussint`.
#'   Specific relevant option:
#'   `num_threads` The maximum number of allowed threads for parallel
#' computing by `excursions::gaussint`, Default: 0, meaning no limit.
#' }
#' \item{optim}{List of control options for `optim()`. Can also include
#'   an element `method` for choosing the optimisation method.}
#' \item{strategy}{The estimation optimisation strategy. Options: "alternating", "joint",
#'   "stepwise". Default: "stepwise"}
#' \item{max_iter}{The maximum number of steps for
#'   `strategy == "alternating"`}
#'   }
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_options

mp_options <- function(...) {
  options <- list(...)
  class(options) <- c("mp_options", "list")
  options
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname as.mp_options

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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_options_default

mp_options_default <- function() {
  mp_options(
    verbose = FALSE,
    optim = list(method = "BFGS"),
    gaussint = list(),
    max_iter = 5,
    strategy = "stepwise",
    hessian = "full"
  )
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' This mid-level function does not access the global and/or default options
#' by itself.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname mp_options_merge

mp_options_merge <- function(...) {
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

  options <- list(...)
  if (length(options) == 0) {
    # Return an empty options object
    options <- mp_options()
  } else if (inherits(options[[1]], "mp_options")) {
    if (length(options) == 1) {
      options <- options[[1]]
    } else if ((length(options) == 2) &&
      inherits(options[[2]], "mp_options")) {
      options <- traverse(options[[1]], options[[2]])
    } else {
      if (any(vapply(
        options[-1],
        function(x) inherits(x, "mp_options"),
        TRUE
      ))) {
        stop("Either 0, 1, or 2 full mp_options inputs allowed")
      }
      options <- traverse(options[[1]], options[-1])
    }
  } else {
    options <- traverse(mp_options(), options)
  }

  options
}


#' @title Set global multiprobit options
#' @description Changes the global options. In many cases, setting options
#' in specific calls to [multiprobit()] is recommended instead.
#' @param ... Zero or more global multiprobit options, see [mp_options()] for details of
#' available options.
#' @param default logical; If `TRUE`, reinitialise the global options with
#' default options before setting new options. Default: `FALSE`
#' @return A copy of the global options is returned invisibly.`
#' @seealso [mp_options()], [mp_options_default()], [mp_options_merge()],
#' [mp_options_get()]
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   mp_options_set(gauusint = list(max.threads = 1),
#'                  optim = list(trace = 5))
#' }
#' }
#' @export
#' @rdname mp_options_set

mp_options_set <- function(..., default = FALSE) {
  envir <- multiprobit_env_get()
  if (default) {
    envir$options <- mp_options_default()
  }
  envir$options <- mp_options_merge(envir$options, ...)
  invisible(mp_options_get())
}
