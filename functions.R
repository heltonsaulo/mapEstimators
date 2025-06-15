
################################################################################
##  Packages
################################################################################

library(ggplot2)
library(tidyverse)
library(dplyr)
library(MASS)
library(nleqslv)


################################################################################
##  Density of the Transformed Exponential Family
################################################################################


# y         : data
# mu, sigma : parameters
# p         : power-transformation parameter
# T_fun     : transformation T(x)
# T_prime   : derivative T'(x)
pdf_exp_family_transformed <- function(y, mu, sigma, p, T_fun, T_prime) {
  y_p       <- y^p
  T_val     <- T_fun(y_p)
  Tp_val    <- T_prime(y_p)
  const     <- p * (mu * sigma)^mu / gamma(mu)
  density   <- const * y^(p - 1) * abs(Tp_val) *
    T_val^(mu - 1) * exp(-mu * sigma * T_val)
  return(density)
}


################################################################################
##  General Analytic Estimator of mu and sigma
################################################################################

estimate_general_sigma_mu <- function(X, p,
                                      T_fun,
                                      T_prime,
                                      T_double_prime,
                                      alpha2, beta2) {
  n            <- length(X)
  log_X        <- log(X)
  T_vals       <- T_fun(X)
  Tp_vals      <- T_prime(X)
  Tpp_vals     <- T_double_prime(X)
  
  # Helpers h1 … h5 (eqs. 7 & 9)
  h1 <- T_vals
  h2 <- log_X + (Tpp_vals / Tp_vals - Tp_vals / T_vals) * X * log_X
  h3 <- (Tp_vals / T_vals) * X * log_X
  h4 <- Tp_vals * X * log_X
  h5 <- beta2 * h3 + (alpha2 - 1) * h4
  
  # Sample means X1 … X5
  X1 <- mean(h1); X2 <- mean(h2)
  X3 <- mean(h3); X4 <- mean(h4)
  X5 <- mean(h5)
  
  # Intermediate terms (eq. 8)
  t1    <- (1 / n) * X5 - (1 + X2) * X1
  t2sq  <- t1^2
  t3    <- 4 * beta2^2 / n * X4 * ((alpha2 - 1) / n * X3 - (1 + X2))
  denom <- 2 * beta2^2 / n * X4
  
  # Sigma_hat (eq. 8)
  sigma_hat <- (t1 + sqrt(t2sq - t3)) / denom
  
  # Mu_hat (eq. 10)
  mu_hat    <- (1 + X2) / (sigma_hat * X4 - X3)
  
  return(list(mu = mu_hat, sigma = sigma_hat, p = p))
}


################################################################################
##  Analytic Estimation for Models in Table 1 (Gamma, Weibull,…) 
################################################################################


estimate_from_table1_model_analytic <- function(X,
                                                model = c("gamma","inv_gamma","weibull","inv_weibull"),
                                                delta = NULL,
                                                p    = 1,
                                                alpha2 = 2, beta2 = 1) {
  model <- match.arg(model)
  if (model %in% c("weibull","inv_weibull") && is.null(delta))
    stop("‘delta’ is required for the chosen model.")
  
  # Define T, T′, T″
  funcs <- switch(model,
                  gamma = list(
                    T =       function(x) x,
                    Tp =      function(x) rep(1, length(x)),
                    Tpp =     function(x) rep(0, length(x))
                  ),
                  inv_gamma = list(
                    T =       function(x) 1 / x,
                    Tp =      function(x) -1 / x^2,
                    Tpp =     function(x) 2 / x^3
                  ),
                  weibull = list(
                    T =       function(x) x^delta,
                    Tp =      function(x) delta * x^(delta - 1),
                    Tpp =     function(x) delta * (delta - 1) * x^(delta - 2)
                  ),
                  inv_weibull = list(
                    T =       function(x) x^(-delta),
                    Tp =      function(x) -delta * x^(-delta - 1),
                    Tpp =     function(x) delta * (delta + 1) * x^(-delta - 2)
                  )
  )
  
  res <- estimate_general_sigma_mu(
    X, p,
    T_fun         = funcs$T,
    T_prime       = funcs$Tp,
    T_double_prime= funcs$Tpp,
    alpha2        = alpha2,
    beta2         = beta2
  )
  res$model <- model
  if (!is.null(delta)) res$delta <- delta
  return(res)
}



################################################################################
##  Random Generation for Each Special Case (Table 1) 
################################################################################

rgeneric_expfam <- function(n, model = c("gamma","inv_gamma","weibull","inv_weibull"),
                            mu, sigma, delta = NULL, p = 1) {
  model <- match.arg(model)
  X <- switch(model,
              gamma =        rgamma(n, shape = mu, rate = mu * sigma),
              inv_gamma =    1 / rgamma(n, shape = mu, rate = mu * sigma),
              weibull = {
                if (is.null(delta)) stop("‘delta’ needed for Weibull.")
                rgamma(n, shape = mu, rate = mu * sigma)^(1 / delta)
              },
              inv_weibull = {
                if (is.null(delta)) stop("‘delta’ needed for Inverse Weibull.")
                1 / rgamma(n, shape = mu, rate = mu * sigma)^(1 / delta)
              }
  )
  # Back-transform to Y = X^(1/p)
  return(X^(1 / p))
}

################################################################################
##  Numerical MAP Estimator via Score Equations + nleqslv 
################################################################################

estimate_map_numerical <- function(
    X,
    T_fun, T_prime, T_double_prime,
    alpha1 = 1, beta1 = 1,
    alpha2 = 1, beta2 = 1,
    init   = c(mu = 1, sigma = 1)
) {
  n <- length(X)
  score  <- function(par) {
    mu    <- par[1]; sigma <- par[2]
    if (mu <= 0 || sigma <= 0) return(rep(1e10,2))
    Yp    <- X; logYp <- log(Yp)
    Tv    <- T_fun(Yp); Tp <- T_prime(Yp); Tpp <- T_double_prime(Yp)
    
    # ∂/∂μ
    s_mu <- n*log(mu)+n*log(sigma)+n - n*digamma(mu) -
      sigma * sum(Tv) + sum(log(Tv)) + (alpha1-1)/mu - beta1
    # ∂/∂σ
    s_sigma <- n*mu/sigma - mu*sum(Tv) + (alpha2-1)/sigma - beta2
    c(s_mu, s_sigma)
  }
  
  sol <- nleqslv(init, score, method="Broyden", control=list(xtol=1e-8))
  if (sol$termcd != 1) warning("MAP did not converge")
  list(mu = sol$x[1], sigma = sol$x[2], convergence = sol$termcd)
}


################################################################################
## MLE for Gamma Model (μ, σ) via L-BFGS-B 
################################################################################

mle_gamma_mu_sigma <- function(x, init = c(mu=1, sigma=1)) {
  if (any(x <= 0)) stop("All x must be positive.")
  nll <- function(par) {
    mu    <- par[1]; sigma <- par[2]
    if (mu <= 0 || sigma <= 0) return(Inf)
    ll <- sum(mu * log(mu * sigma) - lgamma(mu) +
                (mu - 1) * log(x) - mu * sigma * x)
    -ll
  }
  fit <- optim(init, nll, method="L-BFGS-B", lower = c(1e-6,1e-6))
  list(mu_hat      = fit$par[1],
       sigma_hat   = fit$par[2],
       convergence = fit$convergence,
       logLik      = -fit$value)
}


