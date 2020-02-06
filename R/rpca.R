#' @export
rpca <- function(x,
                 k,
                 scale = FALSE,
                 sign_correction = TRUE,
                 rotation = c("none", "varimax", "oblimin"),
                 normalize = TRUE,
                 maxit = 1000,
                 eps = 1e-5,
                 return_data = FALSE) {

  rotation <- match.arg(rotation)

  # compute number of elements
  nsamples <- nrow(x)
  nvariables <- ncol(x)

  # scale and center data
  z <- scale(x, center = TRUE, scale = scale)

  # estimate covariance and perform SVD
  cv <- cov(z)
  # e <- svd(cv, nu = 0L, nv = k)
  #loadings <- e$v %*% diag(sqrt(e$d[1:k]))
  e <- eigen(cv, symmetric = TRUE)
  loadings <-  e$vectors[, 1L:k] %*% diag(sqrt(e$values[1:k]))

  # rotate loadings
  if (rotation == "none") {
    rot <- list(loadings = loadings, rotmat = diag(1, k, k))
  } else if (rotation == "varimax") {
    rot <- varimax(loadings, normalize = normalize, eps = eps)
  } else {
    arglist <-list(loadings,
                   normalize = normalize,
                   maxit = maxit,
                   eps = eps)
    rot <- do.call(getFromNamespace(rotation, 'GPArotation'), arglist)
    rot$rotmat <- t(solve(rot$Th))
  }
  rotmat <- rot$rotmat
  loadings <- unclass(loadings(rot))
  rownames(loadings) <- colnames(x)
  colnames(loadings) <- paste0("Component", 1:k)

  # compute scores
  scores <- z %*% MASS::ginv(t(loadings))
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("Component", 1:k)

  # mean squared deviation
  msd <- sqrt(matrixStats::colMeans2((scores %*% t(loadings) - z)**2))
  names(msd) <- colnames(x)

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
