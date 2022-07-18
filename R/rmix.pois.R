# This function generates samples from finite Poisson mixture
#' rmix.pois
#'
#' @description This function generates samples from the finite Poisson mixture
#' @param n The sample size.
#' @param alpha A vector of the mixing proportions.
#' @param theta A vector of the subpopulation means.
#'
#' @return
#' @export
#'
#' @examples
rmix.pois <- function(n, alpha, theta) {
  #n:      sample size.
  #alpha:  vector of mixing proportions.
  #theta:  vector of subpopulation mean.
  xx =c()
  nindex =  c(rmultinom(1, n, alpha/sum(alpha)))
  for(i in 1:length(alpha)) xx =c(xx, rpois(nindex[i], theta[i]))
  sample(xx) ### avoid separation of subpopulations.
}