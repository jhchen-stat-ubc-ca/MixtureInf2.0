library(quadprog)
#' emtest.thm3
#'
#' @description This function computes a_h as in Theorem 3 of Li and Chen (JASA2009) by MC.
#' It is used in the emtest.exp function.
#' @export
emtest.thm3 <- function(tb, N=10000, tol = 1e-8) {	
  m0 = dim(tb)[1]	
  eig.tb = eigen(tb)
  tb.root = eig.tb[[2]]%*%diag(sqrt(eig.tb[[1]]))%*%t(eig.tb[[2]])
  
  output = c()
  for(i in 1:N) {
    out=solve.QP(Dmat=tb, dvec= tb.root%*%rnorm(m0),
                 Amat= diag(rep(1,m0)), bvec=rep(0,m0))
    output=c(output, sum(out$solution > tol))
  }
  table(output)/N
}