#' pmix.beta
#'
#' @param x The cdf value of the beta mixture at x.
#' @param mix_prop A vector of the mixing proportions.
#' @param alpha A vector of the component alphas.
#' @param beta A vector of the component betas.
#'
#'
#' @examples 
#' x = c(rbeta(50,1,10),rbeta(50,9,0.5))
#' pmix.beta(x,c(.5,.5),c(1,9),c(10,0.5))
#' @export
pmix.beta <- function(x, mix_prop, alpha, beta) {
  vapply(x, function(xx) sum(mix_prop * pbeta(xx, alpha, beta)), numeric(1))
}
