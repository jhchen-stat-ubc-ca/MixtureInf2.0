#' dmix.exp
#'
#' @description This function computes the pdf of the finite exponential mixture distribution.
#' @param x The value at which we computed the cdf.
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the subpopulation means.
#'
#' @return
#' @export
#'
#' @examples
dmix.exp <- function(x, alpha, theta) {
  dmix = 0
  for(i in 1:length(alpha)) {
    dmix = dmix + alpha[i]*dexp(x, rate = 1/theta[i])
  }
  dmix
}