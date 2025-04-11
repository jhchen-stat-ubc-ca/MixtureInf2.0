#' mom.bmm
#'
#' @description A sub-function for \code{pmle.beta}, performs the method of moments (MoM) estimation 
#'              for the beta mixture model.
#'
#'
#' @param x The input data, assumed to follow a beta mixture distribution.
#' @param m0 The number of components in the mixture model.
#' @param seed An optional integer seed for reproducibility.
#' @param maxit Maximum number of iterations allowed in the nonlinear solver.
#'
#' @return A matrix of estimated initial parameters with valid solutions.
#' @export
mom.bmm <- function(x, m0, seed, maxit) {
  repeat {
    mix_porp <- sort(kmeans(x, m0)$size) / length(x)
    theta <- mom.calculation(x, mix_porp, seed)
    
    MM_BMM <- function(params) {
      m0 <- (length(params) + 1) / 3
      mix_porp <- params[1:(m0 - 1)]
      alpha <- params[m0:(2 * m0 - 1)]
      beta <- params[(2 * m0):(3 * m0 - 1)]
      mm <- numeric(3 * m0 - 1)
      
      for (i in 1:(3 * m0 - 1)) {
        for (j in 1:(m0 - 1)) {
          mm[i] <- mm[i] + mix_porp[j] * prod((alpha[j] + 0:(i - 1)) / (alpha[j] + beta[j] + 0:(i - 1)))
        }
        mm[i] <- mm[i] + (1 - sum(mix_porp)) * prod((alpha[m0] + 0:(i - 1)) / (alpha[m0] + beta[m0] + 0:(i - 1)))
      }
      
      observed_sample <- sapply(1:(3 * m0 - 1), function(k) mean(x^k))
      (observed_sample - mm)^2
    }
    
    test_result <- nleqslv::testnslv(c(mix_porp[1:(m0 - 1)], theta), MM_BMM, control = list(maxit = maxit))
    valid <- test_result$out$termcd < 4
    method <- test_result$out$Method[valid]
    global <- test_result$out$Global[valid]
    
    if (length(method) > 0 && length(global) > 0) break
  }
  
  possible_soln <- lapply(seq_along(method), function(i) {
    nleqslv::nleqslv(c(mix_porp[1:(m0 - 1)], theta), MM_BMM,
                     method = method[i], global = global[i],
                     control = list(maxit = maxit))
  })
  
  valid_soln <- sapply(possible_soln, function(sol) {
    (sum(sol$x[1:(m0 - 1)]) < 1) && all(sol$x > 0)
  })
  
  if (!any(valid_soln)) {
    repeat {
      mix_porp <- sort(kmeans(x, m0)$size) / length(x)
      theta <- MoM_Calculation(x, mix_porp, seed)
      possible_soln <- lapply(seq_along(method), function(i) {
        nleqslv::nleqslv(c(mix_porp[1:(m0 - 1)], theta), MM_BMM,
                         method = method[i], global = global[i],
                         control = list(maxit = maxit))
      })
      valid_soln <- sapply(possible_soln, function(sol) {
        (sum(sol$x[1:(m0 - 1)]) < 1) && all(sol$x > 0)
      })
      if (any(valid_soln)) break
    }
  }
  
  init_params <- do.call(cbind, possible_soln[valid_soln])
  return(init_params)
}

#' mom.calculation
#'
#' @description A sub-function for \code{mom.bmm}, generates starting values for the method 
#'              of moments estimation of a beta mixture. 
#'
#' @param x The input data, assumed to be from a beta mixture.
#' @param mix_porp A numeric vector of mixing proportions for each component.
#' @param seed An optional integer seed for reproducibility.
#'
#' @return A numeric vector containing concatenated alpha and beta estimates for each component.
#' @export

mom.calculation <- function(x, mix_porp, seed) {
  if (!is.null(seed)) set.seed(seed)
  
  n_comp <- length(mix_porp)
  alpha <- numeric(n_comp)
  beta <- numeric(n_comp)
  
  for (i in seq_along(mix_porp)) {
    temp_data <- sample(x, size = length(x) * mix_porp[i])
    temp <- mom.beta(temp_data)
    alpha[i] <- temp$alpha
    beta[i]  <- temp$beta
  }
  
  c(alpha, beta)
}

#' mom.beta
#'
#' @description A sub-function for \code{mom.calculation}, estimates the parameters of a single beta 
#'              distribution using the method of moments.
#'
#'
#' @param data A numeric vector assumed to be sampled from a beta distribution.
#'
#' @return A list with components:
#' \item{alpha}{The estimated alpha parameter.}
#' \item{beta}{The estimated beta parameter.}
#' @export

mom.beta <- function(data) {
  mu <- mean(data)
  var_data <- var(data)
  alpha <- mu * (mu * (1 - mu) / var_data - 1)
  beta  <- (1 - mu) * (mu * (1 - mu) / var_data - 1)
  return(list(alpha = alpha, beta = beta))
}
