# Accessor function
#' @export
getConfig <- function() {
    return(.pkg_env$config)
}

# Setter function
#' @export
setConfig <- function(...) {
    new_values <- list(...)
    .pkg_env$config <- utils::modifyList(.pkg_env$config, new_values)
}
