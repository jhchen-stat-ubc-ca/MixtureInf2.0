#' tb2.exp
#'
#' @description This function computes \tilde B_{22} under the Exponential mixture.
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the means of each component.
#'
#' @return
#' @export
#'
#' @examples
tb2.exp <-
  function(alpha,theta)
  {
    m0=length(alpha)
    theta=theta/sum(alpha*theta)
    N=10000
    quan=matrix((0:(N-1)+0.5)/N,ncol=1)
    x=as.numeric(apply(quan,1,qmixexp,alpha=alpha,theta=theta))
    pdf=0
    delta=c()
    y=c()
    z=c()
    
    for(i in 1:m0)
    {
      delta=cbind(delta, dexp(x,rate=1/theta[i])-dexp(x,rate=1/theta[m0]))
      y=cbind(y,(x-theta[i])/theta[i]^2*dexp(x,rate=1/theta[i]))
      z=cbind(z,(x^2-4*theta[i]*x+2*theta[i]^2)/theta[i]^4*dexp(x,rate=1/theta[i])/2)
      pdf=pdf+alpha[i]*dexp(x,rate=1/theta[i])
    }
    bi=cbind(delta[,1:(m0-1)],y,z)/pdf
    B=t(bi)%*%bi
    B11=B[1:(2*m0-1),1:(2*m0-1)]
    B12=B[1:(2*m0-1),(2*m0):(3*m0-1)]
    B22=B[(2*m0):(3*m0-1),(2*m0):(3*m0-1)]
    tB22=B22-t(B12)%*%solve(B11)%*%B12
    diagB22=diag(diag(tB22)^(-1/2),m0,m0)
    corr=diagB22%*%tB22%*%diagB22
    round(corr,3)
  }