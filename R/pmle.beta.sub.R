#' pmle.beta.sub
#'
#' @description A sub function for pmle.beta, does the actual work of PMLE of the beta mixture.
#' It is used in the pmle.beta function.
#' @param x A vector of the observed values.
#' @param m0 The order of the finite mixture model.
#' @param para0 A vector of estimated alpha and beta for every component. 
#' @param an A size control parameter that controls the severity of the penalty. The recommended value is n^{-3/2}.
#' @param epsilon The size of the penalized function of the mixing distribution, default value: epsilon = 1.
#' 
#' @export
pmle.beta.sub <- function(x, m0, para0, an, epsilon) {
  # Extract initial parameters
  mix_porp <- para0[1:m0]
  alpha    <- para0[(m0 + 1):(2 * m0)]
  beta     <- para0[(2 * m0 + 1):(3 * m0)]
  theta    <- para0[(m0 + 1):(3 * m0)]  # concatenation of alpha and beta
  
  n <- length(x)
  
  # E-step: compute component densities and responsibilities
  pdf.sub     <- t(mapply(function(pp, aa, bb) pp * dbeta(x, aa, bb),
                          mix_porp, alpha, beta))
  pdf.mixture <- colSums(pdf.sub) + 1e-100
  ww          <- sweep(pdf.sub, 2, pdf.mixture, FUN = "/")
  
  # Update mixing proportions with penalty
  mix_porp <- (rowSums(ww) + epsilon) / (n + m0 * epsilon)
  
  # M-step: update alpha and beta using penalized maximization
  theta <- Pen.M.Step(x, t(ww), mix_porp, theta, an)
  alpha <- theta[1:m0]
  beta  <- theta[(m0 + 1):(2 * m0)]
  
  # Compute log-likelihood and penalized log-likelihood
  dens <- dmix.beta(x, mix_porp, alpha, beta)
  loglike <- sum(log(dens + 1e-100))
  ploglik <- sum(log(dens) + 1e-100) + sum((log(alpha) - alpha) + (log(beta) - beta)) / an
  
  # Order components by increasing alpha
  ind <- order(alpha)
  
  c(mix_porp[ind], alpha[ind], beta[ind], loglike, ploglik)
}

#' Pen_M_Step
#'
#' @description A sub function for pmle.beta.sub, does the actual work of the M step in the EM algorithm.
#' It is used in the pmle.beta.sub function.
#' @param x A vector of the observed values.
#' @param ww The probability that data point i belongs to component j (matrix with dimensions n x m0).
#' @param mix_porp The mixing proportions for each component.
#' @param theta A vector of estimated alpha and beta for every component.
#' @param an A size control parameter that controls the severity of the penalty.
#' 
#' @export
Pen.M.Step <- function(x, ww, mix_porp, theta, an) {
  k <- length(mix_porp)
  
  ploglikelihood <- function(theta, x) {
    if (any(theta <= 0)) return(NA)
    alpha <- theta[1:k]
    beta  <- theta[(k + 1):(2 * k)]
    # Penalized log-likelihood: sum over all data points and components
    ll <- sum(ww * log(dmix.beta(x, mix_porp, alpha, beta))) +
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
