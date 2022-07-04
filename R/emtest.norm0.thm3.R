#' emtest.norm0.thm3
#'
#' 
#' @description This function computes a_h values in Theorem 3 of Li and Chen (JASA2009) by simulation
#' @param tb The covariance matrix \tilde{B_{22}.
#' @param N The number of repetitions.
#' @param tol A value below which is judged as 0.
#'
#' @return
#' @export
#'
#' @examples
#' @note Need to first install the R package quadprog. 
emtest.norm0.thm3 <- function(tb, N=10000, tol = 1e-8)
{	
  m0 = dim(tb)[1]	
  eig.tb = eigen(tb)
  tb2 = eig.tb$vectors%*%diag(sqrt(eig.tb$values))%*%t(eig.tb$vectors)
  
  output = c()
  for(i in 1:N)
  {
    out=solve.QP(Dmat=tb,dvec=tb2%*%rnorm(m0,0,1),
                 Amat=diag(rep(1,m0)), bvec=rep(0,m0))
    output=c(output, sum(out$solution > tol))
  }
  table(output)/N
}