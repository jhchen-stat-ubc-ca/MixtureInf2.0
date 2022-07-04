#' EMtest.norm.omega
#'
#' @description Computes the overlapping probability of two sub-populations of the Gaussian mixture.
#' @param alpi 
#' @param mui 
#' @param sigi 
#' @param alpj 
#' @param muj 
#' @param sigj 
#' @param tol 
#'
#' @return
#' @export
#'
#' @examples
#' @note Chen, Li and Fu (2012 Jasa). Section 3.1, replaced over1() in the previous package.
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