#' EMtest.norm.penalty
#'
#' @description This function produces the recommended size of penalty of EMtest under the Gaussian mixture.
#' It is used in the emtest.norm function.  
EMtest.norm.penalty <- function(m0, para, nn) {
  if(m0 ==1) an = 0.25
  if(m0 > 3) an = 0.2
  
  if(m0==2) {
    alpha = para[[1]]; mu = para[[2]]; sig = para[[3]]
    omega12 = EMtest.norm.omega(alpha[1],mu[1],sig[1], 
                                alpha[2],mu[2],sig[2])
    omega21 = EMtest.norm.omega(alpha[2],mu[2],sig[2],
                                alpha[1],mu[1],sig[1])
    omega = (omega12 + omega21)/2
    an=-1.859-0.577*log(omega/(1-omega))-60.453/nn
    an=0.35*exp(an)/(1+exp(an))
  }
  
  if(m0==3) {
    alpha = para[[1]]; mu = para[[2]]; sig = para[[3]]
    omega12 = EMtest.norm.omega(alpha[1],mu[1],sig[1], 
                                alpha[2],mu[2],sig[2])
    omega21 = EMtest.norm.omega(alpha[2],mu[2],sig[2],
                                alpha[1],mu[1],sig[1])
    omega1 = (omega12 + omega21)/2
    omega32 = EMtest.norm.omega(alpha[3],mu[3],sig[3], 
                                alpha[2],mu[2],sig[2])
    omega23 = EMtest.norm.omega(alpha[2],mu[2],sig[2],
                                alpha[3],mu[3],sig[3])
    omega2 = (omega23 + omega32)/2
    odds.temp = log(omega1*omega2/(1-omega1)/(1-omega2))
    an = -1.602 - 0.240*odds.temp - 130.394/nn
    an = 0.35*exp(an)/(1+exp(an))
  }
  
  return(an)
}

### Compute overlapping probability of two subpopulations
###         of Gaussian mixture.
###   Chen, Li and Fu (2012 Jasa). Section 3.1
###   replace over1() in the previous package.

EMtest.norm.omega <-
  function(alpi, mui, sigi, alpj, muj, sigj, tol=1e-5)
  {
    sigi=sqrt(sigi); sigj=sqrt(sigj)
    
    if(abs(sigi/sigj-1) < tol) {
      delta = abs(mui-muj)/sigi
      if(delta < tol) delta = tol
      omega = pnorm(-delta/2 + log( alpj / alpi)/delta)
    } else {
      ncp=(mui-muj)*sigi/(sigi^2-sigj^2)
      temp1 = sigj^2*(mui-muj)^2/(sigj^2-sigi^2)^2
      temp2 = 2*sigj^2/(sigi^2-sigj^2)*log(alpi*sigj/alpj/sigi)
      sqrt.value = sqrt(max(temp1 - temp2,0))
      omega = pnorm(sqrt.value-ncp)-pnorm(-sqrt.value-ncp)
      if(sigi<sigj) omega = 1 - omega
    }
    return(omega)
  }