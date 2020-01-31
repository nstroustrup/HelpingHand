#' @name castToMatrix
#'
#' @title Cast a data.table to a matrix
#'
#' @description A wrap-around the \code{dcast} function that returns a \code{matrix} instead of a \code{data.table}
#'
#' @inheritParams data.table::dcast
#'
#' @return A matrix with columns and rows corresponding to the formula.
#'
#' @seealso \link[data.table]{dcast}
#'
#' @import data.table
#' @import magrittr
#' @export
#'
castToMatrix <- function(data,
                         formula,
                         value.var,
                         fun.aggregate,
                         fill = NA,
                         ...) {
  data <- data.table::as.data.table(data)
  if (missing(fun.aggregate)) {
    mat <- data.table::dcast(data = data,
                             formula = formula,
                             value.var = value.var,
                             fill = fill)
  } else {
    mat <- data.table::dcast(data = data,
                             formula = formula,
                             fun.aggregate = fun.aggregate,
                             value.var = value.var,
                             fill = fill)
  }
  var <- as.character(formula)[2]
  rn <- mat[[var]]
  mat <- mat[, (var) := NULL]
  mat <- as.matrix(mat)
  rownames(mat) <- rn
  mat
}
