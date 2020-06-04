#' @name gpca
#'
#' @title GPCA
#'
#' @description Group-wise sparse principal component analysis
#'
#' @param x matrix;
#' @param x numeric matrix; data with variables as columns and samples as rows.
#' @param groups boolean matrix; affiliation matrix. Does the variable in row i pertain to the group in column j?
#' @param trace logical; should progress be printed?
#' @param y numeric matrix; confounder matrix with rows as samples and columns as confounding factors.
#' @param lambda the tuning parameter, non-negative.
#' @param scale logical; should the data be scaled?
#' @param sign_correction logical; should \link{signCorrectionPCA} be called?
#' @param return_data logical; should data be returned?
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional.
#'
#' @return A list with scores and loadings.
#'
#' @export
gpca <- function(x,
                 groups,
                 y = NULL,
                 lambda = NULL,
                 scale = FALSE,
                 center = TRUE,
                 sign_correction = TRUE,
                 return_data = FALSE,
                 trace = FALSE,
                 kernel = "linear",
                 bandwidth = NULL) {

  # compute number of elements
  nsamples <- nrow(x)
  nvariables <- ncol(x)
  ngroups <- ncol(groups)
  nloadings   <- sum(groups)

  # check arguments
  if (nloadings == 0) stop("Invalid 'groups' matrix.")
  if (! is.logical(groups)) stop("Argument 'groups' must be a matrix of FALSE and TRUE")
  if (nrow(groups) != nvariables) stop("Argument 'groups' must have same rows as rows in 'x'.")

  # scale and center data
  z <- scale(x, center = center, scale = scale)

  # perform sparse eigendecomposition
  e <- groupSparseEigen(x = x,
                        groups = groups,
                        trace = trace,
                        y = y,
                        lambda = lambda,
                        kernel = kernel,
                        bandwidth = bandwidth)

  # loadings are eigenvectors
  #   %*% diag(sqrt(e$values), nrow = length(e$values), ncol = length(e$values))
  loadings <- e$vectors
  rownames(loadings) <- colnames(x)
  colnames(loadings) <- colnames(groups)

  # estimate scores from data and loadings
  scores <- z %*% MASS::ginv(t(loadings))
  rownames(scores) <- rownames(x)
  colnames(scores) <- colnames(groups)

  # mean squared deviation
  msd <- sqrt(matrixStats::colMeans2((scores %*% t(loadings) - z)**2))
  names(msd) <- colnames(x)

  # correct sign
  if (sign_correction) {
    s <- signCorrectionPCA(data = z, loadings = loadings, scores = scores)
    loadings <- s$loadings
    scores <- s$scores
  }

  # return results
  results <- list(nsamples = nsamples,
                  loadings = loadings,
                  scores = scores,
                  msd = msd,
                  scaling = attr(z, "scaled:scale"))
  if (return_data) results <- c(results, list(data = z))
  if (sign_correction) results <- c(results, list(signs = s$signs))

  return(results)
}

