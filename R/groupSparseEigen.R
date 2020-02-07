#' @name groupSparseEigen
#'
#' @title Group-wise sparse spectral decomposition of a matrix
#'
#' @description Computes group-wise sparse eigenvalues and eigenvectors of a matrix.
#'
#' @details For each group defined in the \code{groups} parameter,
#' variables not pertaining to the group see their covariance set to zero.
#' Eigendecomposition is then performed on that matrix.
#'
#' @param cv covariance matrix;
#' @param groups boolean matrix; affiliation matrix. Does the variable in row i pertain to the group in column j?
#' @param trace logical; should progress be printed?
#'
#' @return A list containing the sparse eigenvector and eigenvalues.
#'
#' @seealso \link[base]{eigen}
#'
#' @export
#'
groupSparseEigen <- function(cv,
                             groups,
                             trace = FALSE) {

  # number of elements
  nvariables <- nrow(cv)
  ngroups <- ncol(groups)

  # init sparse eigen vectors and values to be returned
  vectors           <- matrix(0, nvariables, ngroups)
  rownames(vectors) <- rownames(cv)
  colnames(vectors) <- colnames(groups)
  values            <- numeric(ngroups)

  # perform set of nested eigen decompositions
  if (trace) pb <- utils::txtProgressBar(min = 0L, max = ngroups, initial = 0L, style = 3L)
  for (k in 1L:ngroups) {
    if (trace) utils::setTxtProgressBar(pb, k)

    # find subset
    in_group <- groups[, k]

    # set out of group correlations to zero
    cvk <- cv
    cvk <- cvk[in_group, ]
    cvk <- cvk[, in_group]

    # perform eigendecomposition and select first vector
    e <- eigen(cvk, symmetric = TRUE)
    vectors[in_group, k] <- e$vectors[, 1L]
    values[k] <- e$values[1L]
  }

  results <- list(vectors = vectors, values = values)
  results
}
