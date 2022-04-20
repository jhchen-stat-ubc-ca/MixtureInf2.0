#' pmle.norm.sub.a
#'
#' @param x The data that can either be a vector or a matrix with the 1st column being the observed values
#' and the 2nd column being the corresponding frequencies.
#' @param m The order of the mixture
#' @param para0 
#' @param lambda The size of the penalized function of the mixing distribution, default value: lambda = 0.
#' @param an 
#'
#' @return
#' @export
#'
#' @examples
pmle.norm.sub.a <- function(x, m, para0, lambda, an)  {
  sn=var(x)
  n = length(x)
  pdf.sub = matrix(0, m, n);
  ww = x*0
  
  alpha = para0[1:m]
  mu = para0[(m+1):(2*m)]
  sigma = para0[(2*m+1):(3*m)]
  
  ## E-step
  for(j in 1:m) pdf.sub[j,] = dnorm(x,mu[j],sqrt(sigma[j]))*alpha[j]+1e-50
  pdf.mixture = apply(pdf.sub,2,sum)
  ## M-step
  for(j in 1:m) {
    ww = pdf.sub[j,]/pdf.mixture;  
    ww.total = sum(ww)
    alpha[j] = (ww.total+lambda)/(n+m*lambda)
    ## modification of Chen, chen and Kalbfleisch
    mu[j] = sum(ww*x)/ww.total
    sigma[j] = (sum(ww*(x-mu[j])^2)+2*sn*an)/(ww.total+2*an)
  }
  ### penalty of Chen, Tan, Zhang
  
  pdf= dmix.norm(x, alpha, mu,sqrt(sigma)) + 1e-100
  ### avoid under flow
  loglike = sum(log(pdf))
  ploglike = loglike + penalty.norm.f(alpha,sigma, sn, an, lambda)
  ###output, increasing subpop mean.
  ind=sort(mu,index.return = TRUE)$ix
  return(c(alpha[ind],mu[ind],sigma[ind],loglike,ploglike))
}