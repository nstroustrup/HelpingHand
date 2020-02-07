#' @name correctedCountsSVA
#'
#' @title Matrix of batch corrected expression values
#'
#' @description Computes batch corrected counts using SVA.
#'
#' @param y matrix; Counts
#' @param mod The model matrix being used to fit the data
#' @param svaseq \link[sva]{svaseq} object
#'
#' @return A matrix of batch corrected counts.
#'
#' @seealso \link[sva]{svaseq}
#'
#' @export
#'
correctedCountsSVA = function(y, mod, svaobj) {
  X = cbind(mod, svaobj$sv)
  Hat = solve(t(X)%*%X)%*%t(X)
  beta = (Hat%*%t(y))
  P = ncol(mod)
  cleany = y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
