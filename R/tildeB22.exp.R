#' tildeB22.exp
#'
#' @description This function computes tilde{B}_{22} under the Exponential mixture.
#' It is used in the emtest.exp function.
tildeB22.exp <- function(alpha, mu) {
  m0 = length(alpha)
  alpha = alpha/sum(alpha)
  mu = mu/sum(alpha*mu)
  N = 10000
  # N = 10000
  quan = matrix((0:(N-1)+0.5)/N, ncol=1)
  x = as.numeric(apply(quan, 1, qmix.exp, alpha=alpha, mu = mu))
  
  pdf = 0
  delta = y = z = c()
  rate0 = 1/mu
  for(i in 1:m0) {
    delta = cbind(delta, dexp(x,rate0[i])-dexp(x,rate0[m0]))
    y = cbind(y,(x-mu[i])/mu[i]^2*dexp(x,rate0[i]))
    z = cbind(z,(x^2-4*mu[i]*x+2*mu[i]^2)/mu[i]^4*dexp(x,rate0[i])/2)
    pdf = pdf + alpha[i]*dexp(x, rate0[i])
  }
  
  bi  = cbind(delta[,1:(m0-1)], y)/pdf
  B11 = t(bi)%*%bi
  z = z - bi%*%solve(B11)%*%t(bi)%*%z
  cor(z)
}