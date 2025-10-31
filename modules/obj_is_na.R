obj_is_na <- function(x) {
  is.atomic(x) && length(x) == 1 && is.na(x)
}