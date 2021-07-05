#' Automatic IFM fit function for vt-D-vines by AIC
#'
#' @param y vector or time series of data to be used for fitting
#' @param margins a vector of marginal distributions to consider (selected using largest likelihood)
#' @param families a vector of families to consider for each pair copula of the d-vine
#' @param k maximum number of generalized lags to consider
#' @param fulcrum a pre-specified fulcrum value. If FALSE, automatic selection using profile likelihood is performed
#' @param omethod optimisation method for the marginal and final fit
#' @param maxit maximum number of iterations for omethod
#'
#' @return An object of class \linkS4class{tscopulafit}.
#' @export
#'
#'
automatic_fit <- function(y,
                          margins = c("norm", "sst", "NIG", "weibull2", "burr2"),
                          families = c(1, 13, 4, 5, 6),
                          k = 4,
                          fulcrum = FALSE,
                          omethod = "Nelder-Mead",
                          maxit = 5000) {
  best_margin <- 0
  marg_lik <- -Inf
  cat(" finding margin...")
  for (m in margins) {
    marg_fit <- fit(margin(m), y, list(method = omethod), control = list(maxit = maxit))
    current_lik <- logLik(marg_fit)
    if (current_lik > marg_lik) {
      marg_lik <- current_lik
      best_marg_fit <- marg_fit
    }
  }
  cat("\r", "finding margin: ", best_marg_fit@margin@name, "\n")
  aic_marg <- AIC(best_marg_fit)
  U <- pmarg(best_marg_fit, y)
  n <- length(U)
  if (fulcrum == F) {
    cat(" finding fulcrum...")
    V0 <- abs(2 * U - 1)
    copsel <- VineCopula::BiCopSelect(u1 = V0[1:(n - 1)], u2 = V0[2:n], familyset = families)
    vtdvine <- vtscopula(dvinecopula(
      family = translate(copsel$family), rotation = translate(copsel$family, what = "rotation"),
      pars = list(c(copsel$par))
    ),
    Vtransform = Vsymmetric()
    )
    pf <- profilefulcrum(U, tscopula = vtdvine, locations = seq(from = 0, to = 1, length = 41), plot = FALSE)
    fulcrum <- as.numeric(pf[which.max(pf[, 2]), 1])
  }
  cat("\r", "finding fulcrum: ", fulcrum, "\n")
  cat(" finding copula families:\n")
  data <- Vlinear(delta = fulcrum)@Vtrans(U, delta = fulcrum)
  taus <- pracf(data, nlags = k, plot = FALSE)$result[2:(k + 1)]
  fam <- c()
  aic <- Inf
  initial_pars <- c()
  try_errors <- 0
  for (lag in 1:k) {
    fam <- c(fam, NA)
    local_best <- -1
    local_aic <- Inf
    for (f in families) {
      fam[lag] <- f
      initial_par_local <- c(initial_pars, VineCopula::BiCopTau2Par(family = f, taus[lag]))
      
      dvmod <- dvinecopula(
        family = translate(fam), rotation = translate(fam, what = "rotation"),
        pars = as.list(initial_par_local)
      )
      copmod <- vtscopula(dvmod, Vtransform = Vlinear(delta = fulcrum))
      cop_fit <- try(fit(copmod, U,
                         list(fulcrum = fulcrum, hessian = FALSE, method = "BFGS"), # internal should always be BFGS for speed
                         control = list(trace = F, maxit = 500)
      ), silent = TRUE)
      
      if (!("try-error" %in% is(cop_fit))) {
        aic_curr <- AIC(cop_fit)
        if (aic_curr < local_aic) {
          local_best <- f
          local_aic <- aic_curr
          local_par_best <- initial_par_local
        }
        if (aic_curr < aic) {
          aic <- aic_curr
          best_fam <- fam
          best_cop_fit <- cop_fit
        }
        cat("\r", "trying lag:", lag,
            ", best family:", translate(best_fam),
            ", prelim AIC:", aic + aic_marg,
            ", try errors:", try_errors,
            sep = " "
        )
      } else {
        try_errors <- try_errors + 1
      }
    }
    if (local_best == -1) {
      local_best <- 6
      local_par_best <- c(initial_pars, 1.00001)
    }
    fam[lag] <- local_best
    initial_pars <- local_par_best
  }
  
  cat("\r", "best lag:", length(best_fam), ", best family:", translate(best_fam),
      ", prelim AIC:", aic + aic_marg, ", try errors:", try_errors,
      sep = " "
  )
  cat("\n")
  cat(" fitting final model...")
  margmod <- margin(best_marg_fit@margin@name)
  x <- tscm(best_cop_fit, margmod)
  fit_f <- fit(x, y,
               list(fulcrum = fulcrum, hessian = TRUE, method = omethod),
               control = list(trace = FALSE, maxit = maxit),
               method = "IFM"
  )
  cat("\r", "fitting final model: done")
  cat("\n")
  invisible(fit_f)
}

#' Automatic IFM fit function for vt-D-vines by MOM
#'
#' @param y vector or time series of data to be used for fitting
#' @param margins a vector of marginal distributions to consider (selected using largest likelihood)
#' @param families a vector of families to consider for each pair copula of the D-vine
#' @param k maximum number of generalized lags to consider
#' @param fulcrum a pre-specified fulcrum value. If FALSE, automatic selection using profile likelihood is performed
#'
#' @return An object of class \linkS4class{tscopula}.
#' @export
#'
#'
automatic_mom <- function(y,
                          margins = c("norm", "sst", "NIG", "weibull2", "burr2"),
                          families = c(1, 13, 4, 5, 6),
                          k = 4,
                          fulcrum = FALSE) {
  best_margin <- 0
  marg_lik <- -Inf
  cat(" finding margin...")
  for (m in margins) {
    marg_fit <- fit(margin(m), y, list(method = "BFGS"), control = list(maxit = 5000))
    current_lik <- logLik(marg_fit)
    if (current_lik > marg_lik) {
      marg_lik <- current_lik
      best_marg_fit <- marg_fit
    }
  }
  cat("\r", "finding margin: ", best_marg_fit@margin@name, "\n")
  aic_marg <- -as.numeric(logLik(best_marg_fit))
  U <- pmarg(best_marg_fit, y)
  n <- length(U)
  if (fulcrum == F) {
    cat(" finding fulcrum...")
    V0 <- abs(2 * U - 1)
    copsel <- VineCopula::BiCopSelect(u1 = V0[1:(n - 1)], u2 = V0[2:n], familyset = families)
    vtdvine <- vtscopula(dvinecopula(family = translate(copsel$family), rotation = translate(copsel$family, what = "rotation"), pars = list(copsel$par)), Vtransform = Vsymmetric())
    pf <- profilefulcrum(U, tscopula = vtdvine, locations = seq(from = 0, to = 1, length = 41), plot = FALSE)
    fulcrum <- as.numeric(pf[which.max(pf[, 2]), 1])
  }
  cat("\r", "finding fulcrum: ", fulcrum, "\n")
  cat(" finding copula families:\n")
  data <- Vlinear(delta = fulcrum)@Vtrans(U, delta = fulcrum)
  taus <- pracf(data, nlags = k, plot = FALSE)$result[2:(k + 1)]
  fam <- c()
  aic <- Inf
  initial_pars <- c()
  
  for (lag in 1:k) {
    fam <- c(fam, NA)
    local_best <- -1
    local_aic <- Inf
    effective_families <- if(taus[lag] > 0) families else c(1, 5)
    for (f in effective_families) {
      fam[lag] <- f
      initial_par_local <- c(initial_pars, VineCopula::BiCopTau2Par(family = f, taus[lag]))
      
      dvmod <- dvinecopula(
        family = translate(fam), rotation = translate(fam, what = "rotation"),
        pars = as.list(initial_par_local)
      )
      copmod <- vtscopula(dvmod, Vtransform = Vlinear(delta = fulcrum))
      cop_fit <- copmod
      aic_curr <- -joint(cop_fit, U)
      if (aic_curr < local_aic) {
        local_best <- f
        local_aic <- aic_curr
        local_par_best <- initial_par_local
      }
      if (aic_curr < aic) {
        aic <- aic_curr
        best_fam <- fam
        best_cop_fit <- cop_fit
      }
      cat("\r", "trying lag:", lag,
          ", best family:", translate(best_fam),
          ", prelim likelihood:", -aic - aic_marg,
          sep = " "
      )
    }
    fam[lag] <- local_best
    initial_pars <- local_par_best
  }
  
  cat("\r", "best lag:", length(best_fam), ", best family:", translate(best_fam),
      ", best likelihood:", -aic - aic_marg,
      sep = " "
  )
  cat("\n")
  x <- tscm(best_cop_fit, best_marg_fit)
  invisible(x)
}

#' Automatic IFM fit function for vt-D-vines by stepwise data transformation
#'
#' @param y vector or time series of data to be used for fitting
#' @param margins a vector of marginal distributions to consider (selected using largest likelihood)
#' @param families a vector of families to consider for each pair copula of the D-vine
#' @param k maximum number of generalized lags to consider
#' @param fulcrum a pre-specified fulcrum value. If FALSE, automatic selection using profile likelihood is performed
#'
#' @return An object of class \linkS4class{tscopula}.
#' @export
#'
#'
automatic_step <- function(y,
                           margins = c("norm", "sst", "NIG", "weibull2", "burr2"),
                           families = c(1, 13, 4, 5, 6),
                           k = 4,
                           fulcrum = FALSE) {
  best_margin <- 0
  marg_lik <- -Inf
  cat(" finding margin...")
  for (m in margins) {
    marg_fit <- fit(margin(m), y, list(method = "BFGS"), control = list(maxit = 5000))
    current_lik <- logLik(marg_fit)
    if (current_lik > marg_lik) {
      marg_lik <- current_lik
      best_marg_fit <- marg_fit
    }
  }
  cat("\r", "finding margin: ", best_marg_fit@margin@name, "\n")
  aic_marg <- -as.numeric(logLik(best_marg_fit))
  U <- pmarg(best_marg_fit, y)
  n <- length(U)
  if (fulcrum == F) {
    cat(" finding fulcrum...")
    V0 <- abs(2 * U - 1)
    copsel <- VineCopula::BiCopSelect(u1 = V0[1:(n - 1)], u2 = V0[2:n], familyset = families)
    vtdvine <- vtscopula(dvinecopula(family = translate(copsel$family), rotation = translate(copsel$family, what = "rotation"), pars = list(copsel$par)), Vtransform = Vsymmetric())
    pf <- profilefulcrum(U, tscopula = vtdvine, locations = seq(from = 0, to = 1, length = 41), plot = FALSE)
    fulcrum <- as.numeric(pf[which.max(pf[, 2]), 1])
  }
  cat("\r", "finding fulcrum: ", fulcrum, "\n")
  cat(" finding copula families:\n")
  data <- Vlinear(delta = fulcrum)@Vtrans(U, delta = fulcrum)
  fam <- c()
  aic <- Inf
  initial_pars <- c()
  taus <- rep(NA, k)
  data <- cbind(data[1:(n - 1)], data[2:n])
  taus[1] <- cor(data, method = "kendall")[1, 2]
  for (lag in 1:k) {
    fam <- c(fam, NA)
    local_best <- -1
    local_aic <- Inf
    effective_families <- if(taus[lag] > 0) families else c(1, 5)
    for (f in effective_families) {
      fam[lag] <- f
      initial_par_local <- c(initial_pars, VineCopula::BiCopTau2Par(family = f, taus[lag]))
      dvmod <- dvinecopula(
        family = translate(fam), rotation = translate(fam, what = "rotation"),
        pars = as.list(initial_par_local)
      )
      copmod <- vtscopula(dvmod, Vtransform = Vlinear(delta = fulcrum))
      cop_fit <- copmod
      aic_curr <- -joint(cop_fit, U)
      if (aic_curr < local_aic) {
        local_best <- f
        local_aic <- aic_curr
        local_par_best <- initial_par_local
        local_cop_best <- cop_fit
      }
      if (aic_curr < aic) {
        aic <- aic_curr
        best_fam <- fam
        best_cop_fit <- cop_fit
      }
      cat("\r", "trying lag:", lag,
          ", best family:", translate(best_fam),
          ", prelim likelihood:", -aic - aic_marg,
          sep = " "
      )
    }
    vinemodel <- local_cop_best@Vcopula
    modelnames <-
      sapply(vinemodel@modelspec, function(v) {
        tolower(v$family)
      })
    pars <- vinemodel@pars
    n <- dim(data)[1]
    model <- rvinecopulib::bicop_dist(modelnames[lag], parameters = pars[[lag]])
    data <-
      cbind(rvinecopulib::hbicop(data[(1:(n - 1)), ], model, cond_var = 2),
            rvinecopulib::hbicop(data[(2:n), ], model, cond_var = 1))
    taus[lag + 1] <- cor(data, method = "kendall")[1, 2]
    
    fam[lag] <- local_best
    initial_pars <- local_par_best
  }
  cat("\r", "best lag:", length(best_fam), ", best family:", translate(best_fam),
      ", best likelihood:", -aic - aic_marg,
      sep = " "
  )
  cat("\n")
  x <- tscm(best_cop_fit, best_marg_fit)
  invisible(x)
}

translate <- function(x, what = "family") {
  result <- x
  VCfam <- c(1, 13, 4, 5, 6)
  RCLfam <- c("gaussian", "clayton", "gumbel", "frank", "joe")
  RCLrot <- c(0, 180, 0, 0, 0)
  if (what == "family") {
    for (i in seq_along(x)) {
      result[i] <- RCLfam[VCfam == x[i]]
    }
  }
  if (what == "rotation") {
    for (i in seq_along(x)) {
      result[i] <- RCLrot[VCfam == x[i]]
    }
  }
  return(result)
}

translate_back <- function(x, rot) {
  result <- x
  VCfam <- c(1, 13, 4, 5, 6)
  RCLfam <- c("gaussian", "clayton", "gumbel", "frank", "joe")
  RCLrot <- c(0, 180, 0, 0, 0)
  for (i in seq_along(x)) {
    result[i] <- VCfam[(RCLfam == x[i]) & (RCLrot == rot[i])]
  }
  return(as.numeric(result))
}
