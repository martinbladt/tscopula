#' KPACF of ARFIMA Process B
#'
#' @param k vector of lags
#' @param theta list with components ar, ma and H specifying the ARFIMA parameters
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_arfima1 <- function(k, theta){
  if (is.list(theta))
    theta <- unlist(theta)
  phi0 <- numeric()
  theta0 <- numeric()
  H <- numeric()
  nm <- substring(names(theta), 1, 5)
  if ("phi" %in% nm)
    phi0 <- theta[nm == "phi"]
  if ("theta" %in% nm)
    theta0 <- theta[nm == "theta"]
  if ("H" %in% nm)
    H <- plogis(theta[nm == "H"])
  acvf <- arfima::tacvfARFIMA(phi = phi0, theta = theta0, H = H, maxlag = k)
  if (is.null(acvf))
    return(rep(NA, k))
  acf <- acvf[-1]/acvf[1]
  pacf <- acf2pacf(acf)
  (2/pi)*asin(pacf)
}

#' KPACF of Exponential Type
#'
#' @param k vector of lags
#' @param theta parameters of exponential decay function
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_exp <- function(k, theta){
  if (is.list(theta))
    theta <- unlist(theta)
  arg <- pmin(theta[1] + theta[2] * (1:k), 0)
  exp(arg)
}

#' KPACF of Power Type
#'
#' @param k vector of lags
#' @param theta parameters of power decay function
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_pow <- function(k, theta){
  if (is.list(theta))
    theta <- unlist(theta)
  arg <- pmin(theta[1] + theta[2] * log(1:k), 0)
  exp(arg)
}

#' KPACF of Transformed Exponential Type
#'
#' @param k vector of lags
#' @param theta parameters of exponential decay function
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_exp2 <- function(k, theta){
  if (is.list(theta))
    theta <- unlist(theta)
  arg <- pmin(theta[1] + theta[2] * (1:k), 0)
  acf <- exp(arg)
  acf2pacf(acf)
}

#' KPACF of Transformed Power Type
#'
#' @param k vector of lags
#' @param theta parameters of power decay function
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_pow2 <- function(k, theta){
  if (is.list(theta))
    theta <- unlist(theta)
  arg <- pmin(theta[1] + theta[2] * log(1:k), 0)
  acf <- exp(arg)
  acf2pacf(acf)
}

