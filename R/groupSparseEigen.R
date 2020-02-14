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
#' @inheritParams gpca
#'
#' @return A list containing the sparse eigenvector and eigenvalues.
#'
#' @seealso \link[base]{eigen}
#'
#' @export
#'
groupSparseEigen <- function(x,
                             groups,
                             y = NULL,
                             lambda = NULL,
                             kernel = "linear",
                             bandwidth = NULL,
                             trace = FALSE) {

  # check whether the number of varaibles  match
  if (ncol(x)!=nrow(groups))
    stop("The numbers of variables in x and groups do not match")

  if (!is.null(y)) {

    if (is.null(dim(y)))
      y <- matrix(y, ncol = 1)

    # check whether the number of samples match
    if (nrow(x)!=nrow(y))
      stop("The numbers of samples in x and y do not match")

    # check whether lambda is non-negative
    if (lambda < 0)
      stop("lambda should be non-negative")

  }

  # number of elements
  nsamples <- nrow(x)
  nvariables <- ncol(x)
  ngroups <- ncol(groups)

  # init sparse eigen vectors and values to be returned
  vectors           <- matrix(0, nvariables, ngroups)
  rownames(vectors) <- colnames(x)
  colnames(vectors) <- colnames(groups)
  values            <- numeric(ngroups)

  # kernel for AC-PCA
  if (! is.null(y))
    kern <- .calKernel(y, kernel, bandwidth)

  # perform set of nested eigen decompositions
  if (trace) pb <- utils::txtProgressBar(min = 0L, max = ngroups, initial = 0L, style = 3L)
  for (k in 1L:ngroups) {
    if (trace) utils::setTxtProgressBar(pb, k)

    # find subset
    in_group <- groups[, k]
    xk <- x[, in_group]

    if (is.null(y)) {
      e <- RSpectra::eigs_sym(cov(xk),
                              k = 1L,
                              which = "LA")
    } else {
      nk <- sum(in_group)
      e <- RSpectra::eigs_sym(.calAv,
                              k = 1L,
                              which = "LA",
                              n = nk,
                              args = list(x = scale(xk, scale = FALSE),
                                          kern = kern,
                                          lambda = lambda))
      e$values <- e$values/(nsamples-1L)
    }

    # perform eigendecomposition and select first vector
    vectors[in_group, k] <- e$vectors[, 1L]
    values[k] <- e$values[1L]
  }

  results <- list(vectors = vectors, values = values)
  results
}

.calKernel <- function(y, kernel, bandwidth){

  if (kernel=="linear"){
    kern <- tcrossprod(y)
  } else if (kernel=="gaussian"){
    if (is.null(bandwidth)==T){
      stop("For gaussian kernel, please specify the bandwidth")
    } else{
      kern <- as.matrix(dist(y, method = "euclidean"))
      kern <- exp(-kern^2/2/bandwidth^2)
    }
  } else {
    stop("Please select a valid kernel, linear kernel or gaussian kernel")
  }
  kern
}

.calAv <- function(v, args){
  x <- args$x
  kern <- args$kern
  lambda <- args$lambda
  crossprod(x, (diag(dim(kern)[1])-lambda*kern)%*%(x%*%matrix(v, ncol=1)))
}
