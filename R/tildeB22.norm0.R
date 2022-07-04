#' tildeB22.norm0
#'
#' @description This function computes \tilde B_{22} matrix in the
#'              EM-test for one-parameter subpopulation distributions
#'              including normal with known variance.
#' @param alpha 
#' @param theta 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
#' @note See Chen and Li (2010JASA).
#'       \tilde B_{22} in this R-function is the correlation matrix.
#'       Renamed from tb2.
tildeB22.norm0 <- function(alpha, theta, N=10000) {
  #alpha:  vector of mixture probabilities.
  #theta:  vector of means of each component.
  m0 = length(alpha)
  quan=matrix((0:(N-1)+0.5)/N,ncol=1)
  qq.x = as.numeric(apply(quan,1,qmix.norm,alpha,theta,rep(1, m0)))
  ## qq.x contain quantiles of the given normal mixture.
  delta = y = z = c()
  
  dnorm.m0 = dnorm(x,theta[m0],1)
  pdf = 0
  for(i in 1:m0) {
    dnorm.i = dnorm(qq.x,theta[i],1)
    delta = cbind(delta, dnorm.i-dnorm.m0)
    y     = cbind(y,(qq.x-theta[i])*dnorm.i)
    z     = cbind(z,((qq.x-theta[i])^2-1)*dnorm.i/2 )
    pdf = pdf+alpha[i]*dnorm.i
  }
  
  bi = cbind(delta[, -m0], y, z)/pdf
  BB = t(bi)%*%bi
  B11 = BB[1:(2*m0-1),1:(2*m0-1)]
  eigen.BB = eigen(B11)[[1]]
  if(eigen.BB[1] >= eigen.BB[2*m0-1]*10^8) {
    err = T; corr = 0
    ## corr is computed only if the matrix does not degnerate.
  } else {
    err = F
    B12 = BB[1:(2*m0-1),(2*m0):(3*m0-1)]
    B22 = BB[(2*m0):(3*m0-1),(2*m0):(3*m0-1)]
    tildeB22 = B22-t(B12)%*%solve(B11)%*%B12
    diagB22  = diag(diag(tildeB22)^(-1/2), m0, m0)
    corr = diagB22%*%tildeB22%*%diagB22
  }
  list(corr, err)
}