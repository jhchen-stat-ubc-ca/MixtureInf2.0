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
pmle.beta <- function(x, m0, n.iter=10, max.iter = 5000, tol = 1e-6,
                      epsilon = 1, an = NULL, maxit=5000) {
  map <- map.beta.mixture(x, m0)
  if (is.null(an)) {
    n <- length(x)
    q <- 2 * m0
    P <- sum(log(map$alpha) - map$alpha + log(map$beta) - map$beta)
    xi <- -(q / 2) / P
    delta <- -log(xi) / log(n)
    an <- n^delta
  }
  uniq_init_params <- unique(t(mom.bmm(x, m0, maxit)))
  map_row <- c(map$mix_prop[1:(m0 - 1)], map$alpha, map$beta)
  starts_params <- rbind(uniq_init_params, MAP = map_row)
  safe.pmle.beta.sub <- function(para0) suppressWarnings(
    tryCatch(pmle.beta.sub(x, m0, para0, an, epsilon), error = function(e) NULL)
  )
  n.starts <- nrow(starts_params)
  results  <- vector("list", n.starts)
  for (i in seq_len(n.starts)) {
    
    para0 <- c(starts_params[i, 1:(m0 - 1)],
               1 - sum(starts_params[i, 1:(m0 - 1)]),
               starts_params[i, m0:(3 * m0 - 1)])
    
    out <- NULL
    for (j in seq_len(n.iter)) {
      step <- safe.pmle.beta.sub(para0)
      if (is.null(step)) { out <- NULL; break }
      para0 <- step[1:(3 * m0)]
      out   <- step
    }
    results[[i]] <- out
  }
  
  results <- Filter(Negate(is.null), results)
  
  output <- do.call(rbind, results)
  index <- which.max(output[, 3 * m0 + 2])
  para0 <- output[index, 1:(3 * m0)]
  ploglik0 <- output[index, 3 * m0 + 2]
  increment <- Inf
  tt <- 0
  while (increment > tol && tt < max.iter) {
    outpara    <- pmle.beta.sub(x, m0, para0, an, epsilon)
    ploglik1  <- outpara[3 * m0 + 2]
    increment  <- ploglik1 - ploglik0
    para0      <- outpara[1:(3 * m0)]
    ploglik0  <- ploglik1
    tt         <- tt + 1
  }
  
  mix_prop <- para0[1:m0]
  alpha    <- para0[(m0 + 1):(2 * m0)]
  beta     <- para0[(2 * m0 + 1):(3 * m0)]
  
  pdf.sub <- matrix(NA_real_, m0, length(x))
  for (k in seq_len(m0)) pdf.sub[k, ] <- mix_prop[k] * dbeta(x, alpha[k], beta[k])
  pdf.mix <- colSums(pdf.sub) + 1e-100
  ww  <- sweep(pdf.sub, 2, pdf.mix, "/")
  classification <- apply(t(ww), 1, which.max)
  
  loglik <- sum(log(dmix.beta(x, mix_prop, alpha, beta) + 1e-100))
  
  list(
    mix_prop      = unname(rousignif(mix_prop)),
    alpha         = unname(rousignif(alpha)),
    beta          = unname(rousignif(beta)),
    loglik        = unname(rousignif(loglik)),
    ploglik       = unname(rousignif(ploglik0)),
    iter.n        = tt,
    classification = classification
  )
}
