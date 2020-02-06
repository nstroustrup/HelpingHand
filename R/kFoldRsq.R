#' @name kFoldRsq
#'
#' @title Cross validated R2 from linear models
#'
#' @description Estimates out of sample R2 and adjusted R2 using using cross-validation
#'
#' @param lmfit object of class \code{lm}
#' @inheritParams bootstrap::crossval
#'
#' @return Raw and cross-validated R2 and adjusted R2
#'
#' @export
#'
kFoldRsq <- function(lmfit, ngroup=10) {

  # adapted from http://www.statmethods.net/stats/regression.html

  mydata <- lmfit$model
  outcome <- names(lmfit$model)[1]
  predictors <- names(lmfit$model)[-1]

  theta.fit <- function(x,y){lsfit(x,y)}
  theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}
  X <- as.matrix(mydata[predictors])
  y <- as.matrix(mydata[outcome])
  n <- nrow(X)
  p <- ncol(X)

  results <- bootstrap::crossval(X, y, theta.fit, theta.predict, ngroup=ngroup)

  raw_rsq <- stats::cor(y, lmfit$fitted.values)**2 # raw R2
  cv_rsq <- stats::cor(y,results$cv.fit)**2 # cross-validated R2
  raw_arsq <- 1 - ((1-raw_rsq) * (n-1) / (n - p - 1))
  cv_arsq <- 1 - ((1-cv_rsq) * (n-1) / (n - p - 1))

  c(raw_rsq=raw_rsq, cv_rsq=cv_rsq, raw_arsq=raw_arsq, cv_arsq=cv_arsq)
}
