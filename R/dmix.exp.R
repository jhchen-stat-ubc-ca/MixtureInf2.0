#' dmix.exp
#'
#' @description This function computes the pdf of the finite exponential mixture distribution.
#' It is used in various functions within the Exponential folder.
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