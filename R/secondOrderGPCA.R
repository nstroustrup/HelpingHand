#' @export
secondOrderGPCA <- function(x,
                            groups,
                            k,
                            scale = FALSE,
                            sign_correction = TRUE,
                            rotation = c("none", "varimax", "oblimin"),
                            normalize = TRUE,
                            maxit = 1000,
                            eps = 1e-5,
                            return_data = FALSE) {

  # group-wise sparse pca
  g <- gpca(x = x,
            groups = groups,
            scale = scale,
            sign_correction = sign_correction,
            return_data = return_data)

  # rotated pca
  rotation <- match.arg(rotation)
  r <- rpca(x = g$scores,
            k = k,
            rotation = rotation,
            scale = FALSE,
            sign_correction = sign_correction,
            normalize = normalize,
            maxit = maxit,
            eps = eps,
            return_data = FALSE)


  # return results
  results <- list(beta = g$loadings,
                  t = g$scores,
                  theta = g$msd,
                  lambda = r$loadings,
                  s = r$scores,
                  psi = r$msd,
                  rotmat = r$rotmat,
                  scaling = g$scaling)
  if (return_data) results <- c(results, list(data = g$data))

  return(results)
}
