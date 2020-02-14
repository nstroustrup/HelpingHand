#' @export
simulateData <- function(nsamples = 1000L,
                         nfactors = 2L,
                         groups_per_factor = 2L,
                         variables_per_group = 10L,
                         psi_rate = 1e2,
                         theta_rate = 1e2,
                         overlap_variables = 0) {

  # number of groups and varaibles
  ngroups    <- nfactors * groups_per_factor
  nvariables <- ngroups * variables_per_group

  # names
  samplesnames <- paste0("Sample", 1:nsamples)
  varnames <- paste0("Variable", 1:nvariables)
  groupnames <- paste0("Group", 1:ngroups)
  factornames <- paste0("Factor", 1:nfactors)

  # simulate beta, loadings from variables to groups
  np <- round(nvariables/ngroups)
  beta <- sapply(1:ngroups, function(i) {
    a <- round((i-1)*np)
    b <- numeric(nvariables)
    s <- sample(c(-1, +1), np, size = np, prob = c(0.2, 0.8))
    v <- rbeta(np, 2, 1)
    b[(a+1):min(a+np, nvariables)] <- s * v
    b
  })

  # certain variables are allowed to be in two groups
  if (overlap_variables > 0 & ngroups > 1) {
    for (i in 1:(ngroups-1)) {
      sel <- sample(which(beta[, i] == 0), size = overlap_variables)
      beta[sel, i] <- rbeta(overlap_variables, 2, 1)
    }
  }
  # beta <- apply(beta, 2, function(x) x/sqrt(sum(x**2))) # make vectors normal
  rownames(beta) <- varnames
  colnames(beta) <- groupnames

  # simulate lambda, loadings from groups to factors
  np <- round(ngroups/nfactors)
  lambda <- sapply(1:nfactors, function(i) {
    a <- round((i-1)*np)
    c(rbeta(a, 1, 10),
      rbeta(np, 10, 1),
      rbeta(ngroups - a - np, 1, 10))
  })
  # lambda <- apply(lambda, 2, function(x) x/sqrt(sum(x**2))) # make vectors normal
  rownames(lambda) <- groupnames
  colnames(lambda) <- factornames

  # simulate stochastic error using exponential
  psi     <- rexp(ngroups, rate = psi_rate)
  names(psi) <- groupnames
  theta   <- rexp(nvariables, rate = theta_rate)
  names(theta) <- varnames

  epsilon <- sapply(1:nvariables, function(i) rnorm(nsamples, 0, theta[i]))
  zeta    <- sapply(1:ngroups, function(i) rnorm(nsamples, 0, psi[i]))

  # simulate factor scores as normal distribution
  s <- replicate(nfactors, rnorm(nsamples, mean = 0, sd = 1))
  rownames(s) <- samplesnames
  colnames(s) <- factornames

  # compute group scores and data
  t <- s %*% t(lambda) + zeta
  rownames(t) <- samplesnames
  colnames(t) <- groupnames

  y <- t %*% t(beta)   + epsilon
  rownames(y) <- samplesnames
  colnames(y) <- varnames

    # compute covariance matrices
  omega <- lambda %*% t(lambda) + diag(psi)
  rownames(omega) <- groupnames
  colnames(omega) <- groupnames

  sigma <-  beta %*% omega %*% t(beta) + diag(theta)
  rownames(sigma) <- varnames
  colnames(sigma) <- varnames

  results <- list(y = y,
                  beta = beta, lambda = lambda,
                  s = s, t = t,
                  omega = omega, sigma = sigma,
                  theta = theta, psi = psi)
  return(results)
}
