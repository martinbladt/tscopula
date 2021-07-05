#' Generalized hyperbolic distribution
#'
#' @param x vector of values
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param lambda overarching shape parameter
#' @param shape shape parameter
#' @param gamma skewness parameter
#' @param mu location parameter
#' @param sigma scale parameter
#' @name Ghyp
#'
NULL
#> NULL
#'
#' @rdname Ghyp
#' @export
pGhyp <- function(q, lambda = 0, shape = 1, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::ghyp(lambda = lambda, alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::pghyp(q, object = obj)
}
#' @rdname Ghyp
#' @export
qGhyp <- function(p, lambda, shape = 1, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::ghyp(lambda = lambda, alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::qghyp(p, object = obj)
}
#' @rdname Ghyp
#' @export
dGhyp <- function(x, lambda, shape = 1, gamma = 0, mu = 0, sigma = 1, log = FALSE) {
  if ((shape < 0) | (sigma < 0)) {
    return(NA)
  }
  obj <- ghyp::ghyp(lambda = lambda, alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  out <- ghyp::dghyp(x, object = obj, logvalue = TRUE)
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}
#' @rdname Ghyp
#' @export
rGhyp <- function(n, lambda, shape = 1, gamma = 0, mu = 0, sigma = 0) {
  obj <- ghyp::ghyp(lambda = lambda, alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::rghyp(n, object = obj)
}

#' Hyperbolic distribution
#'
#' @param x vector of values
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param shape shape parameter
#' @param gamma skewness parameter
#' @param mu location parameter
#' @param sigma scale parameter
#' @name hyp
#'
NULL
#> NULL
#'
#' @rdname hyp
#' @export
phyp <- function(q, shape = 1, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::hyp(alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::pghyp(q, object = obj)
}
#' @rdname hyp
#' @export
qhyp <- function(p, shape = 1, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::hyp(alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::qghyp(p, object = obj)
}
#' @rdname hyp
#' @export
dhyp <- function(x, shape = 1, gamma = 0, mu = 0, sigma = 1, log = FALSE) {
  if ((shape < 0) | (sigma < 0)) {
    return(NA)
  }
  obj <- ghyp::hyp(alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  out <- ghyp::dghyp(x, object = obj, logvalue = TRUE)
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}
#' @rdname hyp
#' @export
rhyp <- function(n, shape = 1, gamma = 0, mu = 0, sigma = 0) {
  obj <- ghyp::hyp(alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::rghyp(n, object = obj)
}

#' Hyperbolic t distribution
#'
#' @param x vector of values
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param nu degree-of-freedom parameter
#' @param gamma skewness parameter
#' @param mu location parameter
#' @param sigma scale parameter
#' @name hypt
#'
NULL
#> NULL
#'
#' @rdname hypt
#' @export
phypt <- function(q, nu = 4, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::student.t(nu = nu, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::pghyp(q, object = obj)
}
#' @rdname hypt
#' @export
qhypt <- function(p, nu = 4, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::student.t(nu = nu, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::qghyp(p, object = obj)
}
#' @rdname hypt
#' @export
dhypt <- function(x, nu = 4, gamma = 0, mu = 0, sigma = 1, log = FALSE) {
  if ((sigma < 0) | (nu < 0)) {
    return(NA)
  }
  obj <- ghyp::student.t(nu = nu, mu = mu, sigma = sigma, gamma = gamma)
  out <- ghyp::dghyp(x, object = obj, logvalue = TRUE)
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}
#' @rdname hypt
#' @export
rhypt <- function(n, nu = 4, gamma = 0, mu = 0, sigma = 0) {
  obj <- ghyp::student.t(nu = nu, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::rghyp(n, object = obj)
}

#' Normal inverse Gaussian distribution
#'
#' @param x vector of values
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param shape shape parameter
#' @param gamma skewness parameter
#' @param mu location parameter
#' @param sigma scale parameter
#' @name NIG
#'
NULL
#> NULL
#'
#' @rdname NIG
#' @export
pNIG <- function(q, shape = 1, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::NIG(alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::pghyp(q, object = obj)
}
#' @rdname NIG
#' @export
qNIG <- function(p, shape = 1, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::NIG(alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::qghyp(p, object = obj)
}
#' @rdname NIG
#' @export
dNIG <- function(x, shape = 1, gamma = 0, mu = 0, sigma = 1, log = FALSE) {
  if ((shape < 0) | (sigma < 0)) {
    return(NA)
  }
  obj <- ghyp::NIG(alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  out <- ghyp::dghyp(x, object = obj, logvalue = TRUE)
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}
#' @rdname NIG
#' @export
rNIG <- function(n, shape = 1, gamma = 0, mu = 0, sigma = 0) {
  obj <- ghyp::NIG(alpha.bar = shape, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::rghyp(n, object = obj)
}

#' Distribution Function of weibull2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pweibull2 <- function(q, sh1 = 1, sc1 = 1, sh2 = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * pweibull(q[q > 0], shape = sh1, scale = sc1)
  result[q <= 0] <- (1 - p0) * pweibull(-q[q <= 0], shape = sh2, scale = sc2, lower.tail = F)
  return(result)
}

#' Quantile Function of weibull2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qweibull2 <- function(p, sh1 = 1, sc1 = 1, sh2 = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- qweibull((p[p > (1 - p0)] - (1 - p0)) / p0, shape = sh1, scale = sc1)
  result[p <= (1 - p0)] <- -qweibull(1 - p[p <= (1 - p0)] / (1 - p0), shape = sh2, scale = sc2)
  return(result)
}

#' Density of weibull2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dweibull2 <- function(x, sh1 = 1, sc1 = 1, sh2 = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * dweibull(x[x > 0], shape = sh1, scale = sc1)
  result[x <= 0] <- (1 - p0) * dweibull(-x[x <= 0], shape = sh2, scale = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for weibull2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rweibull2 <- function(n, sh1 = 1, sc1 = 1, sh2 = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- rweibull(sum(mix), shape = sh1, scale = sc1)
  result[mix == 0] <- -rweibull(sum(1 - mix), shape = sh2, scale = sc2)
  return(result)
}

#' Distribution Function of burr2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pburr2 <- function(q, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * actuar::pburr(q[q > 0], shape1 = sh1, shape2 = sh1b, scale = sc1)
  result[q <= 0] <- (1 - p0) * actuar::pburr(-q[q <= 0], shape1 = sh2, shape2 = sh2b, scale = sc2, lower.tail = F)
  return(result)
}

#' Quantile Function of burr2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qburr2 <- function(p, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- actuar::qburr((p[p > (1 - p0)] - (1 - p0)) / p0, shape1 = sh1, shape2 = sh1b, scale = sc1)
  result[p <= (1 - p0)] <- -actuar::qburr(1 - p[p <= (1 - p0)] / (1 - p0), shape1 = sh2, shape2 = sh2b, scale = sc2)
  return(result)
}

#' Density of burr2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dburr2 <- function(x, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * actuar::dburr(x[x > 0], shape1 = sh1, shape2 = sh1b, scale = sc1)
  result[x <= 0] <- (1 - p0) * actuar::dburr(-x[x <= 0], shape1 = sh2, shape2 = sh2b, scale = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for burr2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rburr2 <- function(n, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- actuar::rburr(sum(mix), shape1 = sh1, shape2 = sh1b, scale = sc1)
  result[mix == 0] <- -actuar::rburr(sum(1 - mix), shape1 = sh2, shape2 = sh2b, scale = sc2)
  return(result)
}

#' Distribution Function of lgamma2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shapelog parameter of right branch
#' @param r1 ratelog parameter of right branch
#' @param sh2 shapelog parameter of left branch
#' @param r2 ratelog parameter of left branch
#'
#' @return
#' @keywords internal
#'
plgamma2 <- function(q, sh1 = 1, r1 = 2, sh2 = 1, r2 = 2, p0 = 0) {
  sh1 <- abs(sh1)
  r1 <- abs(r1)
  sh2 <- abs(sh2)
  r2 <- abs(r2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * actuar::plgamma(q[q > 0] + 1, shapelog = sh1, ratelog = r1)
  result[q <= 0] <- (1 - p0) * actuar::plgamma(-q[q <= 0] + 1, shapelog = sh2, ratelog = r2, lower.tail = F)
  return(result)
}

#' Quantile Function of lgamma2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shapelog parameter of right branch
#' @param r1 ratelog parameter of right branch
#' @param sh2 shapelog parameter of left branch
#' @param r2 ratelog parameter of left branch
#'
#' @return
#' @keywords internal
#'
qlgamma2 <- function(p, sh1 = 1, r1 = 2, sh2 = 1, r2 = 2, p0 = 0) {
  sh1 <- abs(sh1)
  r1 <- abs(r1)
  sh2 <- abs(sh2)
  r2 <- abs(r2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- actuar::qlgamma((p[p > (1 - p0)] - (1 - p0)) / p0, shapelog = sh1, ratelog = r1) - 1
  result[p <= (1 - p0)] <- -actuar::qlgamma(1 - p[p <= (1 - p0)] / (1 - p0), shapelog = sh2, ratelog = r2) + 1
  return(result)
}

#' Density of lgamma2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shapelog parameter of right branch
#' @param r1 ratelog parameter of right branch
#' @param sh2 shapelog parameter of left branch
#' @param r2 ratelog parameter of left branch
#'
#' @return
#' @keywords internal
#'
dlgamma2 <- function(x, sh1 = 1, r1 = 2, sh2 = 1, r2 = 2, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  r1 <- abs(r1)
  sh2 <- abs(sh2)
  r2 <- abs(r2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * actuar::dlgamma(x[x > 0] + 1, shapelog = sh1, ratelog = r1)
  result[x <= 0] <- (1 - p0) * actuar::dlgamma(-x[x <= 0] + 1, shapelog = sh2, ratelog = r2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for lgamma2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shapelog parameter of right branch
#' @param r1 ratelog parameter of right branch
#' @param sh2 shapelog parameter of left branch
#' @param r2 ratelog parameter of left branch
#'
#' @return
#' @keywords internal
#'
rlgamma2 <- function(n, sh1 = 1, r1 = 2, sh2 = 1, r2 = 2, p0 = 0) {
  sh1 <- abs(sh1)
  r1 <- abs(r1)
  sh2 <- abs(sh2)
  r2 <- abs(r2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- actuar::rlgamma(sum(mix), shapelog = sh1, ratelog = r1)
  result[mix == 0] <- -actuar::rlgamma(sum(1 - mix), shapelog = sh2, ratelog = r2)
  return(result)
}

#' Distribution Function of gweibull2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pgweibull2 <- function(q, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * rmutil::pgweibull(q[q > 0], s = sh1, m = sh1b, f = sc1)
  result[q <= 0] <- (1 - p0) * (1 - rmutil::pgweibull(-q[q <= 0], s = sh2, m = sh2b, f = sc2))
  return(result)
}

#' Quantile Function of gweibull2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qgweibull2 <- function(p, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- rmutil::qgweibull((p[p > (1 - p0)] - (1 - p0)) / p0, s = sh1, m = sh1b, f = sc1)
  result[p <= (1 - p0)] <- -rmutil::qgweibull(1 - p[p <= (1 - p0)] / (1 - p0), s = sh2, m = sh2b, f = sc2)
  return(result)
}

#' Density of gweibull2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dgweibull2 <- function(x, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * rmutil::dgweibull(x[x > 0], s = sh1, m = sh1b, f = sc1)
  result[x <= 0] <- (1 - p0) * rmutil::dgweibull(-x[x <= 0], s = sh2, m = sh2b, f = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for gweibull2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rgweibull2 <- function(n, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- rmutil::rgweibull(sum(mix), s = sh1, m = sh1b, f = sc1)
  result[mix == 0] <- -rmutil::rgweibull(sum(1 - mix), s = sh2, m = sh2b, f = sc2)
  return(result)
}

#' Distribution Function of ggamma2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pggamma2 <- function(q, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * rmutil::pggamma(q[q > 0], s = sh1, m = sh1b, f = sc1)
  result[q <= 0] <- (1 - p0) * (1 - rmutil::pggamma(-q[q <= 0], s = sh2, m = sh2b, f = sc2))
  return(result)
}

#' Quantile Function of ggamma2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qggamma2 <- function(p, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- rmutil::qggamma((p[p > (1 - p0)] - (1 - p0)) / p0, s = sh1, m = sh1b, f = sc1)
  result[p <= (1 - p0)] <- -rmutil::qggamma(1 - p[p <= (1 - p0)] / (1 - p0), s = sh2, m = sh2b, f = sc2)
  return(result)
}

#' Density of ggamma2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dggamma2 <- function(x, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * rmutil::dggamma(x[x > 0], s = sh1, m = sh1b, f = sc1)
  result[x <= 0] <- (1 - p0) * rmutil::dggamma(-x[x <= 0], s = sh2, m = sh2b, f = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for ggamma2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rggamma2 <- function(n, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- rmutil::rggamma(sum(mix), s = sh1, m = sh1b, f = sc1)
  result[mix == 0] <- -rmutil::rggamma(sum(1 - mix), s = sh2, m = sh2b, f = sc2)
  return(result)
}

#' Distribution Function of burrgamma
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pburrgamma <- function(q, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * rmutil::pggamma(q[q > 0], s = sh1, m = sh1b, f = sc1)
  result[q <= 0] <- (1 - p0) * actuar::pburr(-q[q <= 0], shape1 = sh2, shape2 = sh2b, scale = sc2, lower.tail = F)
  return(result)
}

#' Quantile Function of burrgamma
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qburrgamma <- function(p, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- rmutil::qggamma((p[p > (1 - p0)] - (1 - p0)) / p0, s = sh1, m = sh1b, f = sc1)
  result[p <= (1 - p0)] <- -actuar::qburr(1 - p[p <= (1 - p0)] / (1 - p0), shape1 = sh2, shape2 = sh2b, scale = sc2)
  return(result)
}

#' Density of burrgamma
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dburrgamma <- function(x, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * rmutil::dggamma(x[x > 0], s = sh1, m = sh1b, f = sc1)
  result[x <= 0] <- (1 - p0) * actuar::dburr(-x[x <= 0], shape1 = sh2, shape2 = sh2b, scale = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for burrgamma
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rburrgamma <- function(n, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- rmutil::rggamma(sum(mix), s = sh1, m = sh1b, f = sc1)
  result[mix == 0] <- -actuar::rburr(sum(1 - mix), shape1 = sh2, shape2 = sh2b, scale = sc2)
  return(result)
}


#' Empirical Distribution Function
#'
#' Compute a version of the empirical distribution function that is
#' standardized by (n+1) to lie strictly in (0,1).
#'
#' @param q vector of values at which to evaluate the empirical distribution
#' function.
#' @param data vector of data values specifying the empirical distribution
#' function
#'
#' @return
#' @export
#'
#' @examples
#' data <- rnorm(100)
#' pedf(c(-2, -1, 0, 1, 2), data)
pedf <- function(q, data) {
  n <- length(data)
  (ecdf(data)(q) * n + 0.5) / (n + 1)
}

#' Quantiles of Empirical Distribution Function
#'
#' @param p vector of probabilities at which to evaluate the quantiles
#' of the empirical distribution function.
#' @param data vector of data values specifying the empirical distribution
#' function
#'
#' @return
#' @export
#'
#' @examples
#' data <- rnorm(100)
#' qedf(c(0.25, 0.5, 0.75), data)
qedf <- function(p, data) {
  quantile(data, p, type = 6)
}
