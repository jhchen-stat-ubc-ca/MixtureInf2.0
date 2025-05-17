#' dmix.beta
#'
#' @description A function that computes the density of the univariate beta mixture.
#' @param x The uni-variate mixture density at x.
#' @param mix_prop A vector of the mixing proportions.
#' @param alpha A vector of the component alphas.
#' @param beta A vector of the component betas.
#' 
#' @examples 
#' dmix.beta(1:10,c(.5,.5),c(2,5),c(10,10))
#' 
#' @export
#' 
#' @note In order to plot the mixture density together with its subpopulation density, one can use the plotmix.beta function.

dmix.beta <- function(x, mix_prop, alpha, beta) {
  dens <- vapply(
    seq_len(length(mix_prop)),
    function(j) mix_prop[j] * dbeta(x, alpha[j], beta[j]),
    numeric(length(x))
  )
  rowSums(dens)
}


