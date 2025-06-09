#' sam.beta.mix
#'
#' @description A general closed-form estimator for finite beta mixture models, based on the 
#'              Score-Adjusted Moment (SAM) method extended to mixtures. This function is also 
#'              used as an initializer for the EM algorithm in the pmle.beta function.
#'
#' @param x A numeric vector of observations in (0, 1).
#' @param m0 Integer; the number of components in the mixture.
#'
#' @return A list containing:
#' \item{mix_prop}{Estimated mixing proportions for each component.}
#' \item{alpha}{Estimated alpha (shape1) parameters.}
#' \item{beta}{Estimated beta (shape2) parameters.}
#'
#'
#' @export

sam.beta.mix <- function(x,m0) {
  kmeans_init <- kmeans(x, m0)
  cluster_assignments <- kmeans_init$cluster
  mix_prop <- kmeans_init$size / sum(kmeans_init$size)
  theta <- sam.calculation(x, cluster_assignments, m0)
  alpha <- theta[1:m0]; beta <- theta[(m0+1):(2*m0)]
  dens <- dmix.beta(x, mix_prop, alpha, beta)
  loglike <- sum(log(dens + 1e-100))
  pdf.sub <- t(mapply(function(pp, aa, bb) pp * dbeta(x, aa, bb),
                      mix_prop, alpha, beta))
  pdf.mixture <- colSums(pdf.sub) + 1e-100
  ww          <- sweep(pdf.sub, 2, pdf.mixture, FUN = "/")
  para0 <- sam.beta.mix.sub(x,ww)
  mix_prop <- para0[1:m0]
  alpha <- para0[(m0 + 1):(2 * m0)]
  beta <- para0[(2 * m0 + 1):(3 * m0)]
  ind <- order(alpha)
  list(
    mix_prop = rousignif(unname(mix_prop[ind])),
    alpha = rousignif(unname(alpha[ind])),
    beta = rousignif(unname(beta[ind]))
  )
}

#' sam.beta.mix.sub
#'
#' @description A sub functino for the main function sam.beta.mix.
#'
#' @param x A numeric vector of observations in (0, 1).
#' @param w A matrix of component-wise weights (responsibilities), with dimensions components × observations.
#' @param epsilon Optional; unused in this implementation but reserved for extensions.
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
  beta  <- bar_yj / denom
  pi <- Wj / sum(Wj)
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
    beta[i]  <- temp$beta
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