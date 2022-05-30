#' dmix.norm
#'
#' @param x uni-variate mixture density at x.
#' @param alpha A vector of the mixing proportions.
#' @param mu A vector of the component means.
#' @param sigma A vector of the component standard deviations.
#' @param tol A fixed value at 1e-8
#'
#' @return
#' @export
#'
#' @examples
dmix.norm <-
  function(x, alpha, mu, sigma, tol=1e-8)
  {
    if(any(alpha<0))
      stop("error: negative mixing proportion")
    if(any(sigma<0))
      stop("error: negative standard deviation")
    m1 = length(alpha); m2 = length(mu); m3 = length(sigma)
    if((m1-m2)^2+(m1-m3)^2 > 0)
      stop("error: differ lengths of alpha, mu and sigma")
    alpha = alpha/sum(alpha)
    dd = x*0
    for(i in 1:m1) dd = dd+alpha[i]*dnorm(x, mu[i], sigma[i])
    return(dd)
  }