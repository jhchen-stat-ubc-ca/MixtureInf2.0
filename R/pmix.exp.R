#' pmix.exp
#'
#' @description This function computes the cdf of the finite exponential mixture distribution
#' @param x The value at which we computed the cdf.
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the subpopulation means.
#'
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples n = 4000
#' theta = c(3, 9, 18)
#' alpha = c(.5, .3, .2)
#' x = rmix.exp(n, alpha, mu)
#' pmix.exp(x, alpha, theta)
pmix.exp <- function(x, alpha, theta) {
  sum(alpha*pexp(x,rate=1/theta))
}