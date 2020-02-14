library(corrplot)
devtools::load_all()

nsamples <- 1000
nfactors <- 4
groups_per_factor <- 5
variables_per_group <- 20


sim <- simulateData(nsamples = nsamples,
                    nfactors = nfactors,
                    groups_per_factor = groups_per_factor,
                    variables_per_group = variables_per_group, overlap_variables = 0,
                    psi_rate = 7,
                    theta_rate = 7)
groups <- sim$beta != 0
dim(sim$y)
colSums(sim$beta**2)
colSums(sim$lambda**2)
det(cov(sim$y))
dim(cor(sim$y))
t(groups) %*% groups
# dg <- findDisjointGroups(groups = groups)
# t(dg) %*% dg

# plot(sim$t)

# replicate(2, secondOrderGPCA(x = sim$y, groups = groups, nfactors = nfactors, rotation = "oblimin"))

# a <- groupSparseEigen(cov(sim$y), group = groups, deflation = TRUE, return_deflation = TRUE)
# b <- groupSparseEigen(cov(sim$y), group = groups, deflation = FALSE)
# sum(a$values - b$values)
#
# sum_eigen <- sum(eigen(cov(sim$y), symmetric = TRUE, only.values = TRUE)$value)
# sum_sparse_eigen <- sum(eigen(a$cv, symmetric = TRUE, only.values = TRUE)$value)
# sum_residual_eigen <- sum(a$values)
# c(sum_eigen = sum_eigen,
#   sum_sparse_eigen = sum_sparse_eigen,
#   sum_residual_eigen = sum_residual_eigen,
#   sum_sum = sum_sparse_eigen + sum_residual_eigen,
#   sum_no_deflation = sum(b$values))

# b <- groupSparseEigen2(cov(sim$y), group = groups)
# sum(a$values - b$values)

# corrplot(cor(sim$y))

# rpca_res <- rpca(x = sim$y, nfactors = 3, rotation = "varimax")
# plot(rpca_res$msd, sim$theta); abline(b=1,a=0)
# cor(as.numeric(abs(rpca_res$scores)), as.numeric(abs(sim$s)))
# plot(rpca_res$scores[, 1], sim$t[, 1])

# gpca_no_res <- gpca(x = sim$y, groups = groups, scale = FALSE, trace = TRUE)
# plot(gpca_no_res$msd, sim$theta); abline(b=1,a=0)
# plot(gpca_no_res$loadings, sim$beta); abline(b=1,a=0); abline(b=-1,a=0)
# plot(as.numeric(abs(gpca_no_res$scores)), as.numeric(abs(sim$t)))
# cor(as.numeric(abs(gpca_no_res$scores)), as.numeric(abs(sim$t)))

library(acPCA)
Y <- matrix(rnorm(nrow(sim$y)), ncol = 1)
acpca <- acPCA::acPCA(X = sim$y,
                      Y = Y,
                      centerX = TRUE, centerY = TRUE,
                      scaleX = FALSE, scaleY = FALSE,
                      lambda = 0.5,
                      kernel = "linear",
                      nPC = 1)
acpca2 <- gpca(x = sim$y,
               y = Y,
               lambda = 0.5,
               kernel = "linear", sign_correction = FALSE,
               groups = matrix(rep_len(TRUE, ncol(sim$y))),
               scale = FALSE)

plot(acpca$v, acpca2$loadings); abline(b=1,a=0); abline(b=-1,a=0)
plot(acpca$Xv, acpca2$scores); abline(b=1,a=0); abline(b=-1,a=0)

# gpca_no_res <- gpca(x = sim$y, y = rep(0, nrow(sim$y)), lambda = 0, groups = groups, scale = FALSE, trace = TRUE)
# plot(gpca_no_res$msd, sim$theta); abline(b=1,a=0)
# plot(gpca_no_res$loadings, sim$beta); abline(b=1,a=0); abline(b=-1,a=0)
# plot(as.numeric(abs(gpca_no_res$scores)), as.numeric(abs(sim$t)))
# cor(as.numeric(abs(gpca_no_res$scores)), as.numeric(abs(sim$t)))

# library(lavaan)
#
# pathway_model <- paste(sapply(colnames(groups), function(g) {
#   p <- groups[, g]
#   paste(g, "=~", paste(names(p)[p], collapse = " + "), "\n")
# }), collapse = "")
# pathway_model
# cat(pathway_model)
#
# lavaan_model <- lavaanify(pathway_model,
#                           meanstructure = FALSE,
#                           auto.fix.first = FALSE,
#                           int.ov.free = TRUE,
#                           int.lv.free = TRUE,
#                           auto.var = TRUE)
# lavaan_model[1:sum(groups), "ustart"] <- gpca_no_res$loadings[groups]
# lavaan_model[(nrow(lavaan_model)-ncol(groups)):nrow(lavaan_model), "ustart"] <- 1
# lavaan_model[(nrow(lavaan_model)-ncol(groups)):nrow(lavaan_model), "free"] <- 0
#
# # fit the model
# dat <- scale(sim$y, center = TRUE, scale = FALSE)
# fit <- lavaan(lavaan_model, data=dat, check.gradient = FALSE)
# fit
# # summary(fit)
# lavaan_beta <- matrix(0, nrow(sim$beta), ncol(sim$beta))
# lavaan_beta[groups] <- fit@ParTable$est[1:sum(groups)]
# plot(gpca_no_res$loadings, sim$beta); abline(b=1,a=0); abline(b=-1,a=0)
#
# plot(lavaan_beta, sim$beta); abline(b=1,a=0)
# plot(lavaan_beta, gpca_no_res$loadings); abline(b=1,a=0)

# gpca_res <- gpca(x = sim$y, groups = groups, scale = FALSE, trace = TRUE, deflation = TRUE)
# plot(gpca_res$msd, sim$theta); abline(b=1,a=0)
# plot(gpca_res$loadings, sim$beta); abline(b=1,a=0); abline(b=-1,a=0)
# cor(as.numeric(abs(gpca_res$scores)), as.numeric(abs(sim$t)))
# plot(gpca_res$scores[, 1], sim$t[, 1])
#
# df <- data.table(NoDeflation = as.numeric(gpca_no_res$loadings),
#            Deflation = as.numeric(gpca_res$loadings),
#            RealValue = as.numeric(sim$beta))
# ggplot(df, aes(Deflation, NoDeflation, color = RealValue)) +
#   geom_point()

# plot(sqrt(matrixStats::colMeans2((gpca_res$scores %*% t(gpca_res$loadings) - sim$y)**2)), sim$theta); abline(b=1,a=0)


# sogpca_res <- secondOrderGPCA(x = sim$y, groups = groups, nfactors = nfactors, rotation = "varimax")
#
# b <- ncol(sogpca_res$lambda)
# cv <- sapply(1:b, function(i) cor(sim$lambda, sogpca_res$lambda[, i]**2))
# o <- apply(cv, 1, which.max)
# sogpca_res$lambda <- sogpca_res$lambda[, o]
# sogpca_res$s <- sogpca_res$s[, o]
#
# plot(sogpca_res$theta, sim$theta); abline(b=1,a=0)
# plot(sogpca_res$beta, sim$beta); abline(b=1,a=0); abline(b=-1,a=0)
# cor(as.numeric(abs(sogpca_res$t)), as.numeric(abs(sim$t)))
# plot(sogpca_res$psi, sim$psi); abline(b=1,a=0)
# plot(sogpca_res$lambda, sim$lambda); abline(b=1,a=0); abline(b=-1,a=0)
# cor(as.numeric(abs(sogpca_res$s)), as.numeric(abs(sim$s)))
#
# signCorrection(sim$y, loadings = sogpca_res$beta, scores = sogpca_res$t)
#
# zeta <- with(sim, lambda %*% t(lambda) + diag(psi))
# sigma <-  with(sim, beta %*% zeta %*% t(beta) + diag(theta))
# zeta_hat <- with(sogpca_res, lambda %*% t(lambda) + diag(psi))
# s_hat <- with(sogpca_res, beta %*% zeta_hat %*% t(beta) + diag(theta))
#
# corrplot(cor(sim$y))
# corrplot(cov2cor(sigma))
# corrplot(cov2cor(s_hat))
#
# corrplot(cov2cor(zeta), is.corr = FALSE)
# corrplot(cov2cor(zeta_hat), is.corr = FALSE)
