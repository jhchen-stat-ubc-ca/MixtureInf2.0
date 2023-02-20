#' dmix.pois
#'
#' @description This function computes the probability mass function of the Poisson mixture.
#' It is used in various functions within the Poisson folder.
dmix.pois <- function(xx, alpha, theta) {
  pmf = xx*0
  for(i in 1:length(alpha)) {
    pmf = pmf + alpha[i]*dpois(xx, theta[i]) }
  pmf
}