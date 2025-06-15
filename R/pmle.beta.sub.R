#' pmle.beta.sub
#'
#' @description A sub-function for \code{pmle.beta}, performs the core computation for 
#'              penalized maximum likelihood estimation (PMLE) of the beta mixture model.
#'
#'
#' @param x A numeric vector of observed values assumed to follow a beta mixture distribution.
#' @param m0 The number of components in the finite mixture model.
#' @param para0 A numeric vector containing initial values for mixing proportions and 
#'        the \eqn{\alpha} and \eqn{\beta} parameters for all components.
#' @param an A penalty control parameter, usually chosen as \eqn{n^{-1/2}}, where \eqn{n} is the sample size.
#' @param epsilon A smoothing parameter for the mixing proportions in the E-step. Default is \code{1}.
#'
#' @return A numeric vector containing updated estimates of the mixing proportions, 
#'         \eqn{\alpha} and \eqn{\beta} values for each component, the log-likelihood, 
#'         and the penalized log-likelihood.
#' @export

pmle.beta.sub <- function(x, m0, para0, an, epsilon) {
  mix_prop <- para0[1:m0]
  alpha    <- para0[(m0 + 1):(2 * m0)]
  beta     <- para0[(2 * m0 + 1):(3 * m0)]
  theta    <- para0[(m0 + 1):(3 * m0)]  
  
  n <- length(x)
  
  # E-step
  pdf.sub     <- t(mapply(function(pp, aa, bb) pp * dbeta(x, aa, bb),
                          mix_prop, alpha, beta))
  pdf.mixture <- colSums(pdf.sub) + 1e-100
  ww          <- sweep(pdf.sub, 2, pdf.mixture, FUN = "/")
  
  mix_prop <- (rowSums(ww) + epsilon) / (n + m0 * epsilon)
  
  # M-step
  theta <- Pen.M.Step(x, t(ww), mix_prop, theta, an)
  alpha <- theta[1:m0]
  beta  <- theta[(m0 + 1):(2 * m0)]
  
  dens <- dmix.beta(x, mix_prop, alpha, beta)
  loglike <- sum(log(dens + 1e-100))
  ploglik <- sum(log(dens + 1e-100)) + sum((log(alpha) - alpha) + (log(beta) - beta)) / an
  
  ind <- order(alpha)
  
  c(mix_prop[ind], alpha[ind], beta[ind], loglike, ploglik)
}

#' Pen.M.Step
#'
#' @description A sub-function for \code{pmle.beta.sub}, performs the M-step in the 
#'              EM algorithm for beta mixture.
#'
#'
#' @param x A numeric vector of observed values assumed to follow a beta mixture distribution.
#' @param ww A matrix of posterior probabilities with dimensions \code{n x m0}, representing 
#'           the probability that each observation belongs to each mixture component.
#' @param mix_prop A numeric vector of mixing proportions for each component.
#' @param theta A numeric vector of current estimates for \eqn{\alpha} and \eqn{\beta} 
#'              parameters for all components.
#' @param an A penalty control parameter, typically \eqn{n^{-1/2}}, that influences the strength 
#'           of the regularization.
#'
#' @return A numeric vector of updated parameter estimates for \eqn{\alpha} and \eqn{\beta}.
#' @export

Pen.M.Step <- function(x, ww, mix_prop, theta, an) {
  k <- length(mix_prop)
  
  ploglikelihood <- function(theta, x) {
    if (any(theta <= 0)) return(NA)
    alpha <- theta[1:k]
    beta  <- theta[(k + 1):(2 * k)]
    logpdf <- vapply(seq_len(k), function(j) {
      log(dbeta(x, alpha[j], beta[j]))
    }, numeric(length(x)))      
    
    ll <- sum(ww * logpdf) +
      sum((log(alpha) - alpha) + (log(beta) - beta)) / an
    ll
  }
  
  pgradlik <- function(theta, x) {
    if (any(theta <= 0)) return(NA)
    alpha <- theta[1:k]
    beta  <- theta[(k + 1):(2 * k)]
    
    g1j <- sapply(1:k, function(j) {
      sum(ww[, j] * (log(x) + digamma(alpha[j] + beta[j]) - digamma(alpha[j]))) -
        ((alpha[j] - 1) / alpha[j]) / an
    })
    
    g2j <- sapply(1:k, function(j) {
      sum(ww[, j] * (log(1 - x) + digamma(alpha[j] + beta[j]) - digamma(beta[j]))) -
        ((beta[j] - 1) / beta[j]) / an
    })
    c(g1j, g2j)
  }
  
  phesslik <- function(theta, x) {
    if (any(theta <= 0)) return(NA)
    alpha <- theta[1:k]
    beta  <- theta[(k + 1):(2 * k)]
    Hess  <- matrix(0, nrow = 2 * k, ncol = 2 * k)
    for (j in 1:k) {
      Hess[j, j]         <- sum(ww[, j] * (trigamma(alpha[j] + beta[j]) - trigamma(alpha[j]))) -
        (1 / alpha[j]^2) / an
      Hess[j + k, j + k] <- sum(ww[, j] * (trigamma(alpha[j] + beta[j]) - trigamma(beta[j]))) -
        (1 / beta[j]^2) / an
      Hess[j, j + k]     <- sum(ww[, j] * trigamma(alpha[j] + beta[j]))
      Hess[j + k, j]     <- Hess[j, j + k]
    }
    Hess
  }
  
  new_est_param <- maxLik::maxNR(fn = ploglikelihood,
                                 grad = pgradlik,
                                 hess = phesslik,
                                 start = theta,
                                 x = x)
  new_est_param$estimate
}
