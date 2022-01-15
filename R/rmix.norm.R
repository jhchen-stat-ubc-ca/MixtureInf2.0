#' rmix.norm
#'
#' @description Generate a random sample from a mixture of normals.
#' @param n The sample size.
#' @param alpha Vector of mixture proportions of length m, the order of the mixture.
#' @param mu Vector of means of component distributions.
#' @param sigma vector of standard deviations of component distributions, default value: sigma =
#' rep(1,length(alpha)).
#'
#' @return It returns a samples of size n from an m-component normal mixture.
#' @export
#'
#' @examples #generate a random sample from a 2 component normal mixture,
#' plot the histogram of the data.
#' x <- rmix.norm(200,c(0.3,0.7),c(-1,2),c(1,2))
#' hist(x)
rmix.norm <-
  function (n,alpha,mu,sigma=rep(1,length(alpha))) 
    #n:      sample size.
    #alpha:  vector of mixing proportion.
    #mu:     vector of means of each component.
    #sigma:  vector of standard deviation of each component.
  {
    if(any(alpha)<0)
      stop("error:mixing proportion must be positive")
    m=length(alpha)
    alpha=alpha/sum(alpha)
    data=c()
    
    nindex=rmultinom(1,n,alpha)
    
    foreach ( i = 1:m) %dopar%
      data=c(data,rnorm(nindex[i],mu[i],sigma[i]))
    data
  }
# Note: add checking option
# check length(mu)=length(sd)=length(alpha)