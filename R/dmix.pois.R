#' dmix.pois
#'
#' @description This function computes the probability mass function of the Poisson mixture
#' @param xx A vector of x at which the probability is computed.
#' @param alpha A vector of the mixing proportions.
#' @param theta A vector of the subpopulation means.
#'
#' @return
#' @export
#'
#' @examples #Compute a density function from a 2 component poisson mixture
#'           theta = c(9, 10)
#'           alpha = c(.4, .6)
#'           dmix.pois(1:20, alpha, theta)
dmix.pois <- function(xx, alpha, theta) {
  pmf = xx*0
  for(i in 1:length(alpha)) {
    pmf = pmf + alpha[i]*dpois(xx, theta[i]) }
  pmf
}