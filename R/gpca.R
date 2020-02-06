#' @export
gpca <- function(x,
                 groups,
                 scale = FALSE,
                 sign_correction = TRUE,
                 return_data = FALSE,
                 trace = FALSE) {

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
  z <- scale(x, center = TRUE, scale = scale)

  # estimate covariance and perform sparse eigendecomposition
  cv <- cov(z)
  e <- groupSparseEigen(cv = cv, groups = groups, trace = trace)

  # loadings are eigenvectors
  loadings <- e$vectors  %*% diag(sqrt(e$values))
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
