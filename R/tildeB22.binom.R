#' tildeB22.binom
#'
#' @description This function computes the \tilde B_{22} matrix in the
#' EM-test for one-parameter subpopulation distributions,
#' including binomial. See Chen and Li (2010JASA). 
#' It is used in the emtest.binom function.
tildeB22.binom <- function(alpha, theta, size) {
  m0=length(alpha)
  x=0:size
  delta = y = z =c()
  
  temp.m0 = dbinom(x,size,theta[m0])
  for(i in 1:m0)  {
    temp.i =  dbinom(x,size,theta[i])
    delta = cbind(delta, temp.i - temp.m0)
    d1binom = (x/theta[i]-(size-x)/(1-theta[i]))*temp.i
    y = cbind(y, d1binom)
    d2binom= (x*(x-1)/theta[i]^2-2*x*(size-x)/(theta[i]*(1-theta[i]))
              +(size-x)*(size-x-1)/((1-theta[i])^2))*temp.i
    z = cbind(z,d2binom/2)
  }
  pmf.x = dmix.binom(x, size, alpha, theta)
  bi=cbind(delta[,1:(m0-1)], y, z)/pmf.x
  B = t(bi)%*%diag(pmf.x)%*% bi
  B11 = B[1:(2*m0-1),1:(2*m0-1)]
  B12 = B[1:(2*m0-1),(2*m0):(3*m0-1)]
  B22 = B[(2*m0):(3*m0-1),(2*m0):(3*m0-1)]
  tB22= B22-t(B12)%*%solve(B11)%*%B12
  ### should check the singularity of B11.
  diagB22=diag(diag(tB22)^(-1/2), m0, m0)
  corr=diagB22%*%tB22%*%diagB22
  corr
}