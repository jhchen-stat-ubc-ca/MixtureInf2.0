#' mom.bmm
#'
#' @description A sub function for pmle.beta, does the actual work of MoM of the beta mixture.
#' It is used in the pmle.beta function.
#' @param x The input data.
#' @param m0 The order/component number of the mixture model.
#' @param seed For reproducible results.
#' @param maxit Maximum amount of iterations.
#' 
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
#' @description A sub function for mom.bmm, generates starting points for the mom.bmm function.
#' It is used in the mom.bmm function.
#' @param x The input data.
#' @param mix_porp The mixing proportion for each component/subpopulation.
#' @param seed For reproducible results.
#' 
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
#' @description A sub function for mom.calculation, calculates the MoM for the beta distribution.
#' It is used in the mom.calculation function.
#' @param data The input data.
#' 
#' @export
mom.beta <- function(data) {
  mu <- mean(data)
  var_data <- var(data)
  alpha <- mu * (mu * (1 - mu) / var_data - 1)
  beta  <- (1 - mu) * (mu * (1 - mu) / var_data - 1)
  return(list(alpha = alpha, beta = beta))
}
