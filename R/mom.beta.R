#' MoM_BMM
#'
#' @description A sub function for pmle.beta, does the actual work of MoM of the beta mixture.
#' It is used in the pmle.beta function.
#' @param x The input data.
#' @param m0 The order/component number of the mixture model.
#' @param seed For reproducible results.
#' @param maxit Maximum amount of iterations.
#' 
#' @export
MoM_BMM <- function(x,m0,seed,maxit) {
  rerun <- TRUE
  while (rerun) {
    mix_porp = sort(kmeans(x,m0)$size)/length(x)
    theta = MoM_Calculation(x, mix_porp, seed)
    
    MM_BMM <- function(params) {
      m0 = (length(params) + 1) / 3
      mix_porp = params[1:(m0 - 1)]
      alpha = params[m0:(2 * m0 - 1)]
      beta = params[(2 * m0):(3 * m0 - 1)]
      mm = numeric(3 * m0 - 1)
      
      for (i in 1:(3 * m0 - 1)) {
        for (j in 1:(m0 - 1)) {
          mm[i] <- mm[i] + mix_porp[j] * prod(sapply(0:(i - 1), function(r) (alpha[j] + r) / (alpha[j] + beta[j] + r)))
        }
        mm[i] <- mm[i] + (1 - sum(mix_porp)) * prod(sapply(0:(i - 1), function(r) (alpha[m0] + r) / (alpha[m0] + beta[m0] + r)))
      }
      
      observed_sample <- sapply(1:(3 * m0 - 1), function(k) mean(x^k))
      return((observed_sample - mm)^2)
    }
    
    test_result = nleqslv::testnslv(c(mix_porp[1:m0-1], theta), MM_BMM, control = list(maxit = maxit))
    
    method = test_result$out[as.numeric(rownames(test_result$out)[which(test_result$out$termcd < 4)]),]$Method
    global = test_result$out[as.numeric(rownames(test_result$out)[which(test_result$out$termcd < 4)]),]$Global
    
    rerun <- identical(method, character(0)) || identical(global, character(0))
  }
  possible_soln=sapply(1:length(method), function(i) {
    nleqslv::nleqslv(c(mix_porp[1:m0-1],theta),MM_BMM, method = method[i],global=global[i],control=list(maxit=maxit)) })
  init_params=possible_soln[,which(sapply(1:ncol(possible_soln), function(i) (sum(possible_soln[,i]$x[1:(m0-1)])<1)&
                                            (sum(possible_soln[, i]$x > 0)==(3*m0-1))))]
  if (is.null(ncol(init_params))==TRUE || ncol(init_params) == 0) {
    rerun_cond = TRUE
    while (rerun_cond) {
      mix_porp = sort(kmeans(x,m0)$size)/length(x)
      theta = MoM_Calculation(x, mix_porp, seed)
      possible_soln=sapply(1:length(method), function(i) {
        nleqslv(c(mix_porp[1:m0-1],theta),MM_BMM, method = method[i],global=global[i],control=list(maxit=maxit)) })
      init_params=possible_soln[,which(sapply(1:ncol(possible_soln), function(i) (sum(possible_soln[,i]$x[1:(m0-1)])<1)&
                                                (sum(possible_soln[, i]$x > 0)==(3*m0-1))))]
      rerun_cond = is.null(ncol(init_params)) || ncol(init_params) == 0
    }
  }
  return(init_params)
}

#' MoM_Calculation
#'
#' @description A sub function for MoM_BMM, generate a starting points for the MoM_BMM function.
#' It is used in the MoM_BMM function.
#' @param x The input data.
#' @param mix_porp The mixing proportion for each component/subpopulation.
#' @param seed For reproducible results.
#' 
#' @export
MoM_Calculation <- function(x, mix_porp, seed) {
  if (is.null(seed)==TRUE) {
    alpha = c()
    beta = c()
    for(i in 1:length(mix_porp)) {
      temp_data = sample(x,length(x)*mix_porp[i])
      temp = MoM_Beta(temp_data)
      alpha = c(alpha, temp$alpha)
      beta = c(beta, temp$beta)
    }
    MoM_Estimates = c(alpha, beta)
    return(MoM_Estimates)
  }
  else {
    set.seed(seed)
    alpha = c()
    beta = c()
    for(i in 1:length(mix_porp)) {
      temp_data = sample(x,length(x)*mix_porp[i])
      temp = MoM_Beta(temp_data)
      alpha = c(alpha, temp$alpha)
      beta = c(beta, temp$beta)
    }
    MoM_Estimates = c(alpha, beta)
    return(MoM_Estimates)
  }
}

#' MoM_Beta
#'
#' @description A sub function for MoM_Calculation, calculates the MoM for the beta distribution.
#' It is used in the MoM_Calculation function.
#' @param data The input data.
#' 
#' @export
MoM_Beta <- function(data) {
  alpha = mean(data)*(mean(data)*(1-mean(data))/var(data) - 1)
  beta = (1-mean(data))*(mean(data)*(1-mean(data))/var(data) - 1)
  return(list(alpha = alpha, beta = beta))
}