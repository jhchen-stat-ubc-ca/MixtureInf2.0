#' tildeB22.pois
#'
#' @description This function computes \tilde B_{22} under the Poisson mixture.
#'              \tilde B_{22} is standardized to the correlation matrix.
#'              It is used in the emtest.pois function.
tildeB22.pois <- function(alpha, theta) {
  m0 = length(alpha)
  x = qpois(0.0001, min(theta)):qpois(0.9999, max(theta))
  delta = y = z = c()
  
  dpois.m0 = dpois(x,theta[m0])
  pdf = 0
  for(i in 1:m0) {
    dpois.i = dpois(x,theta[i])
    delta = cbind(delta, dpois.i - dpois.m0)
    y = cbind(y, (x/theta[i]-1)*dpois.i)
    z = cbind(z, (x*(x-1)/theta[i]^2- 2*x/theta[i]+1)*dpois.i/2)
    pdf = pdf + alpha[i]*dpois.i
  }
  
  bi = cbind(delta[, -m0], y, z)/(pdf+1e-50)
  BB = t(bi)%*%diag(pdf)%*%bi
  B11 = BB[1:(2*m0-1),1:(2*m0-1)]
  eigen.BB = eigen(B11)[[1]]
  if(eigen.BB[1] >= eigen.BB[2*m0-1]*10^8) {
    err = T; corr = 0
  } else {
    B12 = BB[1:(2*m0-1),(2*m0):(3*m0-1)]
    B22 = BB[(2*m0):(3*m0-1),(2*m0):(3*m0-1)]
    tildeB22 = B22 -t(B12)%*%solve(B11)%*%B12
    diagB22 = diag(diag(tildeB22)^(-1/2), m0, m0)
    corr = diagB22%*%tildeB22%*%diagB22
  }
  list("corr"=corr, "Degenerate"= err)
}