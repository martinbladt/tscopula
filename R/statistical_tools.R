#' New Generic for carrying out a likelihood ratio test between two Time Series Models
#'
#' Methods are available for objects of class \linkS4class{tscopulafit},
#' \linkS4class{tscmfit}
#'
#' @param x,y objects of the model class.
#' @param ...
#'
#' @return a likelihood ratio test result.
#' @export
#'
#'
setGeneric("LRT", function(x, y, ...) {
  standardGeneric("LRT")
})

#' LRT Method for tscmfit Class
#'
#' @param x,y objects of class \linkS4class{tscfit}.
#'
#' @return LRT between the models.
#' @export
#'
setMethod("LRT", c(x = "tscmfit", y = "tscmfit"), function(x, y) {
  xnames <- names(x@fit$par)
  ynames <- names(y@fit$par)
  if (!(all(xnames %in% ynames) || all(ynames %in% xnames))) {
    stop("Non-nested models")
  }
  LR <- 2 * abs(logLik(y) - logLik(x))
  degrees <- abs(attributes(logLik(y))$df - attributes(logLik(x))$df)
  return(c(LR = LR, p.val = pchisq(LR, df = degrees, lower.tail = FALSE)))
})

#' LRT Method for tscopulafit Class
#'
#' @param x,y objects of class \linkS4class{tscopulafit}.
#'
#' @return LRT between the models.
#' @export
#'
setMethod("LRT", c(x = "tscopulafit", y = "tscopulafit"), function(x, y) {
  LRT(as(x, "tscmfit"), as(y, "tscmfit"))
})

#' Rank Autocorrelation Function
#'
#' @param data a vector of data to which the autocorrelation function is to be computed
#' @param method one of three ways to compute the acf: "kendall", "spearman" or "pearson"
#' @param nlags the number of lags considered
#' @param plot logical specifing whether to plot the acf or not
#' @param main title of the plot
#' @param xlab x-axis label of the plot
#' @param ylab y-axis label of the plot
#' @param lwd width of lines
#'
#' @return a vector containing the autocorrelation values
#' @export
#'
#' @examples
#' racf(runif(1000))
#' y <- sim(dvinecopula("gauss", list(0.5)))
#' racf(y)
racf <- function(data, method = "kendall",
                 nlags = floor(5 * log10(length(data))),
                 plot = TRUE,
                 main = paste("Series", deparse(substitute(data)), sep = " "),
                 xlab = "Lag",
                 ylab = "ACF",
                 lwd = 2) {
  n <- length(data)
  result <- numeric(nlags + 1)
  for (i in 1:(nlags + 1)) result[i] <- cor(data[1:(n - i + 1)], data[i:n], method = method)
  if (plot) {
    plot(0:nlags, result,
      type = "h",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ylim = c(min(-0.02, min(result) - 0.02), 1),
      lwd = lwd
    )
    abline(h = 0, lwd = lwd)
    if (method == "kendall") {
      abline(h = 1.96 * sqrt(2 * (2 * n + 5) / (9 * n * (n - 1))), col = "blue", lty = 2)
      abline(h = -1.96 * sqrt(2 * (2 * n + 5) / (9 * n * (n - 1))), col = "blue", lty = 2)
    }
    if (method %in% c("spearman", "pearson")) {
      abline(h = tanh(1.96 / sqrt(n - 3)), col = "blue", lty = 2)
      abline(h = tanh(-1.96 / sqrt(n - 3)), col = "blue", lty = 2)
    }
  }
  invisible(result)
}

#' Partial Rank Autocorrelation Function
#'
#' @param data a vector of data to which the autocorrelation function is to be computed
#' @param method one of three ways to compute the acf: "kendall", "spearman" or "pearson"
#' @param nlags the number of lags considered
#' @param plot logical specifing whether to plot the acf or not
#' @param main title of the plot
#' @param xlab x-axis label of the plot
#' @param ylab y-axis label of the plot
#' @param lwd width of lines
#' @param exact logical indicating whether the rpacf should be exact or fast
#'
#' @return a list containing the partial autocorrelation values and pvalues
#' @export
#'
#' @examples
#' pracf(runif(1000))
#' y <- sim(dvinecopula("gauss", list(0.5)))
#' pracf(y)
pracf <- function(data, method = "kendall",
                  nlags = floor(5 * log10(length(data))),
                  plot = TRUE,
                  main = paste("Series", deparse(substitute(data)), sep = " "),
                  xlab = "Lag",
                  ylab = "PACF",
                  lwd = 2) {
  n <- length(data)
  result <- numeric(nlags + 1)
  for (i in 0:nlags) {
    d <- matrix(NA, n - i, i + 1)
    for (j in 1:(n - i)) d[j, ] <- data[j:(j + i)]
    Pcor <- ppcor::pcor(d, method = method)
    result[i + 1] <- Pcor$estimate[1, i + 1]
  }
  if (plot) {
    plot(0:nlags, result,
      type = "h",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ylim = c(min(-0.02, min(result) - 0.02), 1),
      lwd = lwd
    )
    abline(h = 0, lwd = lwd)
    if (method == "kendall") {
      abline(h = 1.96 * sqrt(2 * (2 * n + 5) / (9 * n * (n - 1))), col = "blue", lty = 2)
      abline(h = -1.96 * sqrt(2 * (2 * n + 5) / (9 * n * (n - 1))), col = "blue", lty = 2)
    }
    if (method %in% c("spearman", "pearson")) {
      abline(h = tanh(1.96 / sqrt(n - 3)), col = "blue", lty = 2)
      abline(h = tanh(-1.96 / sqrt(n - 3)), col = "blue", lty = 2)
    }
  }
  invisible(list(result = result))
}

#' logLik Method for uGARCHfit Class
#'
#' @param object an object of class \linkS4class{uGARCHfit}.
#'
#' @return an object of class logLik
#' @import rugarch
#'
setMethod("logLik", "uGARCHfit", function(object) {
  ll <- likelihood(object)
  attr(ll, "nobs") <- object@model$modeldata$T
  attr(ll, "df") <- length(coef(object))
  class(ll) <- "logLik"
  ll
})

#' Extract BiCop objects from Fitted tscopula
#'
#' @param fit an object of class \linkS4class{tscopulafit} or \linkS4class{tscmfit}.
#'
#' @return a list of BiCop objects.
#' @export
#'
getBiCop <- function(fit) {
  copula <- fit@tscopula
  if (is(copula, "vtscopula")) {
    copula <- copula@Vcopula
  }
  specs <- copula@modelspec
  pars <- copula@pars
  output <- vector("list", length(specs))
  for (i in seq_along(output)) {
    output[[i]] <- BiCop_wrapper(specs[[i]]$family, pars[[i]])
  }
  output
}
