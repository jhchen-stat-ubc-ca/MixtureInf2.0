#' EMtest.norm.penalty
#'
#' @param m0 The order
#' @param para The list of mixing proportion and its mean and var
#' @param nn The sample size
#'
#' @return
#' @export
#'
#' @examples
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