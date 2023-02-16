#' qmix.exp
#'
#' @description This function computes the quantile of the finite exponential mixture distribution
#' @param alpha A vector of the mixing probabilities.
#' @param mu A vector of the subpopulation means.
#' @param qq The level of the quantile to be computed.
#'
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples n = 4000
#' theta = c(3, 9, 18)
#' alpha = c(.5, .3, .2)
#' x = rmix.exp(n, alpha, mu)
#' qmix.exp(.99, alpha, theta)
qmix.exp <- function(qq, alpha, mu) {
  #alpha:  vector of mixture probabilities.
  #mu:  vector of subpopulation means.
  #  qq :  the level of the quantile to be computed.
  
  qmix.exp.sub <- function(x, alpha, mu, qq) 
  { sum(alpha*pexp(x, rate=1/mu))-qq }
  
  if((qq<0)|(qq>1)) stop("quantile level out of range")
  
  uniroot(qmix.exp.sub, c(0, qexp(qq, rate=1/max(mu))), 
          alpha = alpha, mu = mu, qq=qq, tol=1e-10)$root 
}