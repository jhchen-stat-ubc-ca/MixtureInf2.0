#' dmix.norm
#'
#' @description A function that computes the density of the univariate normal mixture
#' @param xx The uni-variate mixture density at x.
#' @param alpha A vector of the mixing proportions.
#' @param mu A vector of the component means.
#' @param sigma A vector of the component standard deviations.
#' @param tol A fixed value at 1e-8
#'
#' @return 
#' @export
#'
#' @examples mu = c(8, 9)
#' alpha = c(.2, .8)
#' sigma = rep(1,length(alpha))
#' n = 1000
#' x = rmix.norm(n,alpha,mu)
#' dmix.norm(1:10,alpha,mu,sigma)
dmix.norm <- function(xx, alpha, mu, sigma) {
  if(any(alpha<0)) stop("error: negative mixing proportion")
  if(any(sigma<0)) stop("error: negative standard deviation")
  m1 = length(alpha); m2 = length(mu); m3 = length(sigma)
  if((m1-m2)^2+(m1-m3)^2 > 0) 
    stop("error: differ lengths of alpha, mu and sigma")
  dd = xx*0
  for(i in 1:m1) dd = dd+alpha[i]*dnorm(x, mu[i], sigma[i])
  return(dd)
}