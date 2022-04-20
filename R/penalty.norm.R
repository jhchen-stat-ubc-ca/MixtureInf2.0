#' penalty.norm.f
#'
#' @param alpha A vector of the mixture proportions of length m, m represents the order of the mixture
#' @param sigma A vector of the component standard deviations.
#' @param sigma0 
#' @param an 
#' @param lambda The level of penalty, default value is 0.
#'
#' @return
#' @export
#'
#' @examples
penalty.norm.f <- function(alpha,sigma,sigma0,an,lambda)
{
  -an*sum(sigma0/sigma + log(sigma/sigma0))+lambda*sum(log(alpha))
}