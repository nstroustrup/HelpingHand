#' Overlapping estimation
#'
#' Estimates overlapping area of two distributions.
#'
#' @param x numerical vector
#' @param y numerical vector
#' @param nbins number of equally spaced points at which the overlapping density is evaluated see density for details
#' @param ... optional arguments to be passed to function density
#'
#' @return Estimates of overlapped areas relative to each pair of distributions
#'
#' @seealso \link[overlapping]{overlap}
#'
#' @describeIn overlap
#'
#' @export
#'
overlap <- function(x, y, nbins = 1024, ...) {

  dx <- density(x, n = nbins)
  fx <- approxfun(dx$x, y = dx$y)

  dy <- density(y, n = nbins)
  fy <- approxfun(dy$x, y = dy$y)

  dat <- data.table::data.table(X = c(dx$x, dy$x))
  dat[, Y1 := fx(X)]
  dat[, Y2 := fy(X)]
  dat[is.na(dat)] <- 0
  dat[, MIN := pmin(Y1, Y2)]
  dat[, MAX := pmax(Y1, Y2)]

  over <- with(dat, sum(MIN, na.rm = TRUE)/sum(MAX, na.rm = TRUE))
  over
}


#' @param mean mean of normal distribution
#' @param sigma standard deviation of normal distribution
#'
#' @describeIn overlap
#'
#' @export
#'
overlapNormal <- function(x, mean, sd, nbins = 1024, ...) {

  dx <- density(x, n = nbins)
  fx <- approxfun(dx$x, y = dx$y)

  fy <- function(z) dnorm(z, mean, sd)
  dyx <- seq(mean - 5*sd, mean + 5*sd, nbins)

  dat <- data.table::data.table(X = c(dx$x, dyx))
  dat[, Y1 := fx(X)]
  dat[, Y2 := fy(X)]
  dat[is.na(dat)] <- 0
  dat[, MIN := pmin(Y1, Y2)]
  dat[, MAX := pmax(Y1, Y2)]

  over <- with(dat, sum(MIN, na.rm = TRUE)/sum(MAX, na.rm = TRUE))
  over
}

#' @param meanlog mean of log-normal distribution in log scale
#' @param sdlog standard deviation of log-normal distribution in log scale
#'
#' @describeIn overlap
#'
#' @export
#'
overlapLogNormal <- function(x, meanlog, sdlog, nbins = 1024, ...) {

  dx <- density(log(x), n = nbins)
  fx <- approxfun(dx$x, y = dx$y)

  fy <- function(z) dlnorm(exp(z), meanlog, sdlog)
  dyx <- seq(meanlog - 5*sdlog, meanlog + 5*sdlog, nbins)

  dat <- data.table::data.table(X = c(dx$x, dyx))
  dat[, Y1 := fx(X)]
  dat[, Y2 := fy(X)]
  dat[is.na(dat)] <- 0
  dat[, MIN := pmin(Y1, Y2)]
  dat[, MAX := pmax(Y1, Y2)]

  over <- with(dat, sum(MIN, na.rm = TRUE)/sum(MAX, na.rm = TRUE))
  over
}
