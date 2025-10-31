obj_is_na <- function(x) {
  return(is.atomic(x) && length(x) == 1 && is.na(x))
}

obj_is_na_or_NULL <- function(x) {
  return(obj_is_na(x) || is.null(x))
}

is.integer <- function(x) {
  return(is.numeric(x) && x %% 1 == 0)
}

is.boolean <- function(x) {
  return(is.atomic(x) &&
           !is.na(x) &&
           !is.null(x) &&
           (x == T || x == F))
}

absolute_package_filename <- function(filename) {
  return(system.file(filename, package = config$package_name))
}

absolute_package_directory <- function(directory) {
  return(paste0(system.file(directory, package = config$package_name),"/"))
}

package_readRDS <- function(filename) {
  absolute_filename <- system.file(filename, package = config$package_name)
  return(readRDS(absolute_filename))
}