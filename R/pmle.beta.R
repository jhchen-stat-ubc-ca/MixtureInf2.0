#' pmle.beta
#'
#' @description This function computes the penalized maximum likelihood estimate (PMLE) of the 
#' parameters under a finite mixture of beta distributions.
#'
#'
#' @param x A numeric vector of observed values.
#' @param m0 The number of components (order) in the beta mixture model.
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
pmle.beta <- function(x, m0, n.iter = 10, max.iter = 5000, tol = 1e-6,
                      epsilon = 1, an = NULL, maxit = 5000) {
  if (is.null(an)) an <- length(x)^(3/2)
  uniq_init_params <- unique(t(mom.bmm(x, m0, maxit)))
  safe.pmle.beta.sub <- function(para0) suppressWarnings(
    tryCatch(pmle.beta.sub(x, m0, para0, an, epsilon), error = function(e) NULL)
  )
  results <- vector("list", nrow(uniq_init_params))
  for (i in seq_len(nrow(uniq_init_params))) {
    para0 <- c(
      uniq_init_params[i, 1:(m0 - 1)],
      1 - sum(uniq_init_params[i, 1:(m0 - 1)]),
      uniq_init_params[i, m0:(3 * m0 - 1)]
    )
    out <- NULL
    for (j in seq_len(n.iter)) {
      step <- safe.pmle.beta.sub(para0)
      if (is.null(step)) break
      para0 <- step[1:(3 * m0)]
      out   <- step
    }
    results[[i]] <- out
  }
  results <- Filter(Negate(is.null), results)
  if (!length(results)) stop("All initializations failed")
  output <- do.call(rbind, results)
  index <- which.max(output[, 3 * m0 + 2])
  para0 <- output[index, 1:(3 * m0)]
  ploglik <- output[index, 3 * m0 + 2]
  increment <- Inf
  tt <- 0
  while (increment > tol && tt < max.iter) {
    step <- safe.pmle.beta.sub(para0)
    if (is.null(step)) break
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