#' rmix.pois
#'
#' @description This function generates a random sample from a finite Poisson mixture
#' @param n The sample size.
#' @param alpha A vector of the mixing proportions.
#' @param theta A vector of the subpopulation means.
#'
#' @return
#' @export
#'
#' @examples #Generate a random sample from a 2 component Poisson mixture,
#'and computes the sample mean and variance. 
#'x <- rmix.pois(200,c(0.3,0.7),c(2,5))
#'mean(x)
#'var(x)

rmix.pois <- function(n, alpha, theta) {
  #n:      sample size.
  #alpha:  vector of mixing proportions.
  #theta:  vector of subpopulation mean.
  xx =c()
  nindex =  c(rmultinom(1, n, alpha/sum(alpha)))
  for(i in 1:length(alpha)) xx =c(xx, rpois(nindex[i], theta[i]))
  sample(xx) ### avoid separation of subpopulations.
}