pmix.norm <- function(x, alpha, mu, sigma)
{
  #cdf value of Normal mixture at x.
  #alpha:  vector of mixing proportions.
  #mu:     vector of component means.
  #sigma:  vector of component standard deviations.
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
  mixture.cdf = x*0
  for(i in 1:m1) {
    mixture.cdf = mixture.cdf+alpha[i]*pnorm(x,mu[i],sigma[i])
  }
	return(mixture.cdf)
}
### Curated Feb 7, 2022