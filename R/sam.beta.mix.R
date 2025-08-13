#' sam.beta.mix
#'
#' @description A general closed-form estimator for finite beta mixture models, based on the 
#'              Score-Adjusted Moment (SAM) method extended to mixtures. This function is also 
#'              used as an initializer for the EM algorithm in the pmle.beta function.
#'
#' @param x A numeric vector of observations in (0, 1).
#' @param m0 Integer; the number of components in the mixture.
#' @param tol The tolerance threshold for convergence based on change in penalized log-likelihood. 
#'            Default is \code{1e-6}.
#' @param max.iter The maximum number of iterations allowed in the final EM optimization phase.
#' @param epsilon A regularization parameter controlling the smoothing in the E-step. Default is \code{1}.
#' @param an A size control parameter that determines the severity of the penalty. 
#'           The recommended value is \eqn{n^{-1/2}}. If \code{NULL}, it is set to \code{length(x)^(1/2)}.
#'
#' @return A list containing:
#' \item{mix_prop}{Estimated mixing proportions for each component.}
#' \item{alpha}{Estimated alpha (shape1) parameters.}
#' \item{beta}{Estimated beta (shape2) parameters.}
#'
#'
#' @export
sam.beta.mix <- function(x, m0, tol = 1e-6, max.iter = 5000, epsilon=1, an=NULL) {
  if (is.null(an)) an <- length(x)^(1/2)
  k0 <- kmeans(x, m0)
  mix_prop <- k0$size / sum(k0$size)
  th <- sam.calculation(x, k0$cluster, m0)
  alpha <- th[1:m0]; beta <- th[(m0 + 1):(2 * m0)]
  dens <- dmix.beta(x, mix_prop, alpha, beta)
  ploglik <- sum(log(dens + 1e-100)) + 
    sum((log(alpha) - alpha) + (log(beta) - beta)) / an
  diff <- Inf; tt <- 0
  while (diff > tol && tt < max.iter) {
    pdf.sub <- t(mapply(function(p, a, b) p * dbeta(x, a, b),
                        mix_prop, alpha, beta))
    pdf.mix <- colSums(pdf.sub) + 1e-100
    w <- sweep(pdf.sub, 2, pdf.mix, "/")
    para <- sam.beta.mix.sub(x, w, epsilon)
    mix_prop.new <- para[1:m0]
    alpha.new <- para[(m0 + 1):(2 * m0)]
    beta.new <- para[(2 * m0 + 1):(3 * m0)]
    dens.new <- dmix.beta(x, mix_prop.new, alpha.new, beta.new)
    ploglik.new <- sum(log(dens.new + 1e-100)) + 
      sum((log(alpha.new) - alpha.new) + (log(beta.new) - beta.new)) / an
    diff <- ploglik.new - ploglik
    if (diff > tol) {
      mix_prop <- mix_prop.new
      alpha <- alpha.new
      beta <- beta.new
      ploglik <- ploglik.new
      tt <- tt + 1
    } else break
  }
  o <- order(alpha)
  list(mix_prop = rousignif(unname(mix_prop[o])),
       alpha = rousignif(unname(alpha[o])),
       beta = rousignif(unname(beta[o])),
       ploglik = rousignif(ploglik),
       iter.n = tt)
}

#' sam.beta.mix.sub
#'
#' @description A sub functino for the main function sam.beta.mix.
#'
#' @param x A numeric vector of observations in (0, 1).
#' @param w A matrix of component-wise weights.
#' @param epsilon A regularization parameter controlling the smoothing in the E-step. Default is \code{1}.
#'
#' @return A numeric vector of length 3m: concatenated estimates of
#' \item{pi}{Mixing proportions.}
#' \item{alpha}{Component-wise alpha parameters.}
#' \item{beta}{Component-wise beta parameters.}
#'
#'
#' @export
sam.beta.mix.sub <- function(x, w, epsilon) {
  W <- t(w)
  m <- ncol(W)
  y <- 1 - x
  Wj <- colSums(W)
  bar_xj <- colSums(W * x) / Wj
  bar_yj <- 1 - bar_xj
  log_x <- log(x)
  log_y <- log(y)
  mean_log_xj <- colSums(W * log_x) / Wj
  mean_log_yj <- colSums(W * log_y) / Wj
  mean_xlogxj <- colSums(W * (x * log_x)) / Wj
  mean_ylogyj <- colSums(W * (y * log_y)) / Wj
  denom <- mean_xlogxj - bar_xj * mean_log_xj + mean_ylogyj - bar_yj * mean_log_yj 
  alpha <- bar_xj / denom
  beta <- bar_yj / denom
  pi <- (Wj+epsilon) / (sum(Wj)+m*epsilon)
  c(pi = pi, alpha = alpha, beta = beta)
}

#' sam.calculation
#'
#' @description Computes unweighted SAM estimates of beta parameters for each mixture component, 
#'              based on initial hard assignments (e.g., from K-means clustering).
#'
#' @param x A numeric vector of data in (0, 1).
#' @param cluster_assignments Integer vector assigning each observation to one of the m0 components.
#' @param m0 Number of mixture components.
#'
#' @return A numeric vector: \code{c(alpha_1, ..., alpha_m0, beta_1, ..., beta_m0)}.
#'
#'
#' @export

sam.calculation <- function(x, cluster_assignments, m0) {
  alpha <- numeric(m0)
  beta <- numeric(m0)
  
  for (i in 1:m0) {
    cluster_data <- x[cluster_assignments == i]
    temp <- sam.beta(cluster_data)
    alpha[i] <- temp$alpha
    beta[i] <- temp$beta
  }
  
  c(alpha, beta)
}

#' sam.beta
#'
#' @description Estimates the parameters of a single beta distribution using the Score-Adjusted 
#'              Moment (SAM) estimator proposed by Chen and Xiao (2025).
#'
#' @param x A numeric vector sampled from a single beta distribution.
#'
#' @return A list with:
#' \item{alpha}{Estimated alpha (shape1) parameter.}
#' \item{beta}{Estimated beta (shape2) parameter.}
#'
#' @details This estimator replaces the classical second moment equation with a mixed moment 
#'          condition involving log terms to obtain robust, closed-form estimates. It performs 
#'          well for both small and large samples and avoids iterative optimization.
#'
#' @references Chen, P., & Xiao, X. (2025). Novel closed-form point estimators for the beta distribution.
#'             *Statistical Theory and Related Fields*, 9(1), 12–33.
#'
#' @export

sam.beta <- function(x) {
  y <- 1 - x
  bar_x <- mean(x)
  bar_y <- 1 - bar_x
  mean_log_x <- mean(log(x))
  mean_log_y <- mean(log(y))
  mean_xlogx <- mean(x * log(x))
  mean_ylogy <- mean(y * log(y))
  denom <- mean_xlogx - bar_x * mean_log_x + mean_ylogy - bar_y * mean_log_y
  list(alpha = bar_x / denom, beta = bar_y / denom)
}