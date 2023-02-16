#' rmix.pois
#'
#' @description This function generates a random sample from a finite Poisson mixture
#' @param n The sample size.
#' @param alpha A vector of the mixing proportions.
#' @param theta A vector of the subpopulation means.
#'
#' @return It returns a sample of size n from an m-component Poisson mixture.
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples n = 3000
#' mu = c(3, 9, 15)
#' alpha = c(.2, .3, .5)
#' xx = rmix.pois(n, alpha, mu)
rmix.pois <- function(n, alpha, theta) {
  #n:      sample size.
  #alpha:  vector of mixing proportions.
  #theta:  vector of subpopulation mean.
  xx =c()
  nindex =  c(rmultinom(1, n, alpha/sum(alpha)))
  for(i in 1:length(alpha)) xx =c(xx, rpois(nindex[i], theta[i]))
  sample(xx) ### avoid separation of subpopulations.
}