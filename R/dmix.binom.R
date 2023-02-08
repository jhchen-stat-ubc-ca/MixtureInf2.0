#' dmix.binom
#'
#' @description This function computes the probability mass function points in vector x.
#' @param x A vector at which the probability mass function is evaluated.
#' @param size The number of trials.
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the probabilities of success of each component.
#'
#' @return
#' @export
#'
#' @examples > Compute a density function from a 3 component binomial mixture
#' alpha = c(.5, .1, .4)
#' theta = c(.1, .35, .8)
#' size = 25
#' x = rmix.binom(1000, size, alpha, theta)
#' dmix.binom(x,size,alpha,theta)
dmix.binom <- function(x, size, alpha, theta) {
  pmf.x = x*0
  m0 = length(alpha)
  for(i in 1:m0) { pmf.x = pmf.x + dbinom(x, size, theta[i])*alpha[i]}
  pmf.x
}
