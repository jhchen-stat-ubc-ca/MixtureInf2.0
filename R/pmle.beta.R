#' pmle.beta
#'
#' @description This function computes the penalized maximum likelihood estimate (PMLE) of the 
#' parameters under a finite mixture of beta distributions.
#'
#'
#' @param x A numeric vector of observed values.
#' @param m0 The number of components (order) in the beta mixture model.
#' @param n.init The number of additional random initializations to generate. Default is 5.
#' @param n.iter The number of EM iterations to perform for each initial value. The initialization 
#'               yielding the highest penalized log-likelihood will be further optimized.
#' @param max.iter The maximum number of iterations allowed in the final EM optimization phase.
#' @param tol The tolerance threshold for convergence based on change in penalized log-likelihood. 
#'            Default is \code{1e-6}.
#' @param epsilon A regularization parameter controlling the smoothing in the E-step. Default is \code{1}.
#' @param an A size control parameter that determines the severity of the penalty. 
#'           The recommended value is \eqn{n^{-3/2}}. If \code{NULL}, it is set to \code{length(x)^(3/2)}.
#' @param maxit The maximum number of iterations for the nonlinear solvers used in the 
#'              method of moments estimation during initialization.
#'
#' @return A list with the following components:
#' \item{mix_prop}{Estimated mixing proportions.}
#' \item{alpha}{Estimated \eqn{\alpha} parameters for each component.}
#' \item{beta}{Estimated \eqn{\beta} parameters for each component.}
#' \item{loglik}{Log-likelihood at the PMLE.}
#' \item{ploglik}{Penalized log-likelihood at the PMLE.}
#' \item{iter.n}{Number of EM iterations performed after initialization.}
#' \item{classification}{Component assignment for each observation based on posterior probabilities.}
#'
#' @examples
#' data <- c(rbeta(50, 10, 2), rbeta(50, 3, 18))
#' pmle.beta(data, 2)
#' @export
pmle.beta <- function(x, m0, n.init = 5, n.iter = 10, max.iter = 5000, tol = 1e-6,
                      epsilon = 1, an = NULL, maxit = 5000) {
  if (is.null(an)) an <- length(x)^(1/2)
  sar_init <- sar.beta.mix(x,m0)
  sam_init <- sam.beta.mix(x,m0)
  all_init_params <- rbind(c(sar_init$mix_prop,sar_init$alpha,sar_init$beta),
                           c(sam_init$mix_prop,sam_init$alpha,sam_init$beta))
  for (i in seq_len(n.init)) {
    group_assign <- sample(1:m0, size = length(x), replace = TRUE)
    mix_prop <- as.numeric(table(factor(group_assign, levels = 1:m0))) / length(x)
    shapes_sam <- sam.calculation(x, group_assign, m0)
    all_init_params <- rbind(all_init_params, c(mix_prop, shapes_sam))
  }
  for (i in seq_len(n.init)) {
    group_assign <- sample(1:m0, size = length(x), replace = TRUE)
    mix_prop <- as.numeric(table(factor(group_assign, levels = 1:m0))) / length(x)
    shapes_sar <- sar.calculation(x, group_assign, m0)
    all_init_params <- rbind(all_init_params, c(mix_prop, shapes_sar))
  }
  results <- vector("list", nrow(all_init_params))
  for (i in seq_len(nrow(all_init_params))) {
    para0 <- all_init_params[i,]
    out <- NULL
    for (j in seq_len(n.iter)) {
      output <- pmle.beta.sub(x, m0, para0, an, epsilon)
      para0 <- output[1:(3 * m0)]
      out <- output
    }
    results[[i]] <- out
  }
  candidates <- do.call(rbind, results)
  index <- which.max(candidates[, 3 * m0 + 2])
  para0 <- candidates[index, 1:(3 * m0)]
  ploglik <- candidates[index, 3 * m0 + 2]
  increment <- Inf
  tt <- 0
  while (increment > tol && tt < max.iter) {
    step <- pmle.beta.sub(x, m0, para0, an, epsilon)
    para0 <- step[1:(3 * m0)]
    increment <- step[3 * m0 + 2] - ploglik
    ploglik <- step[3 * m0 + 2]
    tt <- tt + 1
  }
  mix_prop <- para0[1:m0]
  alpha <- para0[(m0 + 1):(2 * m0)]
  beta <- para0[(2 * m0 + 1):(3 * m0)]
  pdf.sub <- matrix(NA_real_, nrow = m0, ncol = length(x))
  for (k in seq_len(m0)) pdf.sub[k, ] <- mix_prop[k] * dbeta(x, alpha[k], beta[k])
  pdf.mixture <- colSums(pdf.sub) + 1e-100
  ww <- sweep(pdf.sub, 2, pdf.mixture, "/")
  classification <- apply(t(ww), 1, which.max)
  loglik <- sum(log(dmix.beta(x, mix_prop, alpha, beta) + 1e-100))
  list(
    mix_prop = rousignif(mix_prop),
    alpha = rousignif(alpha),
    beta = rousignif(beta),
    loglik = rousignif(loglik),
    ploglik = rousignif(ploglik),
    iter.n = tt,
    classification = classification
  )
}