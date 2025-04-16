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
#' @param seed An optional integer seed for reproducibility. Default is \code{NULL}.
#' @param maxit The maximum number of iterations for the nonlinear solvers used in the 
#'              method of moments estimation during initialization.
#'
#' @return A list with the following components:
#' \item{mix_porp}{Estimated mixing proportions.}
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
pmle.beta <- function(x, m0, n.iter = 10, max.iter = 5000, tol = 1e-6, epsilon = 1, 
                      an = NULL, seed = NULL, maxit = 5000) {
  
  if (is.null(an)) an <- length(x)^(3/2)
  
  init_params <- mom.bmm(x, m0, seed, maxit)
  uniq_init_params <- unique(do.call(rbind, lapply(1:ncol(init_params), 
                                                   function(i) init_params[, i]$x)))
  
  output_list <- vector("list", nrow(uniq_init_params))
  for (i in 1:nrow(uniq_init_params)) {
    para0 <- c(uniq_init_params[i, 1:(m0 - 1)], 
               1 - sum(uniq_init_params[i, 1:(m0 - 1)]),
               uniq_init_params[i, m0:(3 * m0 - 1)])
    for (j in 1:n.iter) {
      outpara <- pmle.beta.sub(x, m0, para0, an, epsilon)
      para0 <- outpara[1:(3 * m0)]  
    }
    output_list[[i]] <- outpara
  }
  
  output_mat <- do.call(rbind, output_list)
  best_index <- which.max(output_mat[, (3 * m0 + 2)])
  para0 <- output_mat[best_index, 1:(3 * m0)]
  
  squarem_res <- SQUAREM::squarem(
    par = para0,
    fixptfn = function(p) {
      ret <- pmle.beta.sub(x, m0, p, an, epsilon)
      ret[1:(3 * m0)]
    },
    control = list(tol = tol, maxiter = max.iter)
  )
  
  para0 <- squarem_res$par
  
  outpara <- pmle.beta.sub(x, m0, para0, an, epsilon)
  
  mix_porp <- outpara[1:m0]
  alpha <- outpara[(m0 + 1):(2 * m0)]
  beta  <- outpara[(2 * m0 + 1):(3 * m0)]
  
  pdf.sub <- t(mapply(function(pp, aa, bb) pp * dbeta(x, aa, bb),
                      mix_porp, alpha, beta))
  pdf.mixture <- colSums(pdf.sub) + 1e-100
  ww <- sweep(pdf.sub, 2, pdf.mixture, "/")
  
  list(mix_porp = rousignif(mix_porp),
       alpha = rousignif(alpha),
       beta = rousignif(beta),
       loglik = rousignif(outpara[3 * m0 + 1]),
       ploglik = rousignif(outpara[3 * m0 + 2]),
       iter.n = squarem_res$iter,
       classification = apply(t(ww), 1, which.max))
}