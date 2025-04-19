#' dmix.beta
#'
#' @description A function that computes the density of the univariate beta mixture.
#' @param x The uni-variate mixture density at x.
#' @param mix_porp A vector of the mixing proportions.
#' @param alpha A vector of the component alphas.
#' @param beta A vector of the component betas.
#' 
#' @examples 
#' dmix.beta(1:10,c(.5,.5),c(2,5),c(10,10))
#' 
#' @export
#' 
#' @note In order to plot the mixture density together with its subpopulation density, one can use the plotmix.beta function.

dmix.beta <- function(x, mix_porp, alpha, beta) {
  mix_porp <- as.numeric(mix_porp)
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  x <- as.vector(x)
  
  sapply(x, function(xi) {
    sum(mix_porp * dbeta(xi, alpha, beta))
  })
}

