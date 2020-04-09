#' @export
factorToNum <- function(x) {
  as.numeric(levels(x))[x]
}
