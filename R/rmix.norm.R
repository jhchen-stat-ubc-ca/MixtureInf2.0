#' rmix.norm
#'
#' @description This function generates a random sample from a mixture of normal.
#' @param n The sample size.
#' @param alpha A vector of the mixture proportions of length m, also the order of the mixture.
#' @param mu A vector of the component distributions' mean.
#' @param sigma A vector of the component distributions' standard deviation, default value: sigma =
#' rep(1,length(alpha)).
#'
#' @return It returns a samples of size n from an m-component normal mixture.
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples #generate a random sample from a 2 component normal mixture,
#' plot the histogram of the data.
#' x <- rmix.norm(200,c(0.3,0.7),c(-1,2),c(1,2))
#' hist(x)
#' @export
rmix.norm <- function (n, alpha, mu, sigma=rep(1,length(alpha))) 
  #n:      sample size.
  #alpha:  vector of mixing proportions.
  #mu:     vector of component means.
  #sigma:  vector of component standard deviations.
{
  if(any(alpha<0))
    stop("error: negative mixing proportion")
  if(any(sigma<0))
    stop("error: negative standard deviation")
  m1 = length(alpha)
  m2 = length(mu)
  m3 = length(sigma)
  if((m1-m2)^2+(m1-m3)^2 > 0)
    stop("error: differ lengths of alpha, mu and sigma")
  alpha = alpha/sum(alpha)
  data = c()
  nindex = rmultinom(1, n, alpha)
  for(i in 1:m1)
    data=c(data, rnorm(nindex[i], mu[i], sigma[i]))
  sample(data)   ### avoid separation of subpopulations.
}