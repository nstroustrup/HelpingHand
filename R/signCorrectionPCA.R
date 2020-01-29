#' @name signCorrectionPCA
#'
#' @title Sign correction for principal component analysis
#'
#' @description Corrects the sign ambiguity that exists in PCA by correlating scores with the data.
#'
#' @details The sign of the jth component is defined as:
#' \code{sum(cor(scores[, j], data) * loadings[, j]**2))}
#'
#' @param data matrix; The data matrix with samples in the rows and variables in the columns.
#' @param loadings matrix; The loading matrix with variables in the rows and components in the columns.
#' @param scores matrix; The score matrix with samples in the rows and components in the columns.
#'
#' @return A list containing the sign corrected loadings and scores.
#'
#' @seealso \link[stats]{prcomp}
#'
#' @export
#'
signCorrectionPCA <- function(data, loadings, scores) {

  signs <- sapply(1:ncol(scores), function(j)
    sum(cor(scores[, j], data) * loadings[, j]**2))

  # plot(data.frame(data[, loadings[, 3] != 0], scores[, 3]))

  for (i in 1:ncol(scores)) {
    if (signs[i] < 0) {
      loadings[, i] <- - loadings[, i]
      scores[, i] <- - scores[, i]
    }
  }

  # plot(data.frame(data[, loadings[, 3] != 0], scores[, 3]))

  return(list(loadings = loadings, scores = scores, signs = signs))
}
