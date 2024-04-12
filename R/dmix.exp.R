#' dmix.exp
#'
#' @description This function computes the pdf of the finite exponential mixture distribution.
#' @param x The value at which we computed the cdf.
#' @param alpha A vector of the mixing probabilities.
#' @param mu A vector of the subpopulation means.
#' @param logObs logObs=T if the data is in logarithm
#' 
#' @examples alpha = c(.3, .7)
#' mu = c(1.8, 0.9)
#' n = 2000
#' x = rmix.exp(n, alpha = alpha, theta = mu)
#' dmix.exp(x,alpha,mu)
#' 
#' @export
#' 
#' @note In order to plot the mixture density together with its subpopulation density, one can use the plotmix.exp function.
dmix.exp <- function(x, alpha, mu, logObs = F) {
  dmix = 0
  for(i in 1:length(alpha)) {
    dmix = dmix + alpha[i]*dexp(x, rate = 1/mu[i])
  }
  if(logObs) {
    dmix = 0
    for(i in 1:length(alpha)) {
      dmix = dmix + alpha[i]*exp(x-exp(x)/mu[i])/mu[i]
    }
  }
  dmix
}