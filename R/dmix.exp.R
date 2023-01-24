#' dmix.exp
#'
#' @description This function computes the pdf of the finite exponential mixture distribution.
#' @param x The value at which we computed the cdf.
#' @param alpha A vector of the mixing probabilities.
#' @param mu A vector of the subpopulation means.
#' @param logObs logObs=T if the data is in logarithm
#'
#' @return
#' @export
#'
#' @examples
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