#' dmix.norm
#'
#' @description A function that computes the density of the univariate normal mixture.
#' It is used in various functions within the Gaussian folder.
dmix.norm <- function(xx, alpha, mu, sigma) {
  if(any(alpha<0)) stop("error: negative mixing proportion")
  if(any(sigma<0)) stop("error: negative standard deviation")
  m1 = length(alpha); m2 = length(mu); m3 = length(sigma)
  if((m1-m2)^2+(m1-m3)^2 > 0) 
    stop("error: differ lengths of alpha, mu and sigma")
  dd = xx*0
  for(i in 1:m1) dd = dd+alpha[i]*dnorm(xx, mu[i], sigma[i])
  return(dd)
}