#' emtest.binom
#'
#' @description This function computes the EM-test statistic and the p-value for H0: m=m0.
#' @param x The input data, either a vector or a matrix columns with two columns: count and freq. 
#' @param size The number of trials of the binomial.
#' @param m0 The order under the null hypothesis.
#' @param CC A optional tuning parameter for the EM-test procedure. 
#' @param lambda A tuning for the fitting null model.
#' @param init.val The initial values for the PMLE under the null model. 
#' @param n.init The number of initial values chosen for the EM-algorithm.
#' @param n.iter The least amount of iterations for all initial values.
#' @param max.iter The maximum amount of iterations allowed.
#' @param tol The tolerance value for the convergence of the EM-algorithm.
#' @param k The number of EM iterations in order to obtain the EM-test statistic.
#' @param rformat The format of the output. See rforma().
#'
#' @return
#' @export
#'
#' @examples alpha = c(.5, .1, .4)
#' theta = c(.1, .35, .8)
#' size = 25
#' x = rmix.binom(1000, size, alpha, theta)
#' y = as.matrix(table(x))
#' count = as.numeric(rownames(y))
#' freq = y[,1]
#' emtest.binom(x, size, m0=2, CC=NULL, lambda = 0, init.val=NULL,
#' n.init=10, n.iter=50, max.iter=200, tol=1e-6, k=3, rformat=FALSE)
emtest.binom <- function(x, size, m0, CC=NULL, lambda = 0, init.val=NULL,
                         n.init=10, n.iter=50, max.iter = 2000, tol=1e-6, k=3, rformat=FALSE)
{
  if(is.vector(x)) {
    y = as.matrix(table(x))
    count = as.numeric(rownames(y))
    freq = y[,1]
  }
  if(is.matrix(x)) {
    count = x[,1]
    freq  = x[,2]
  }
  n = sum(freq)
  
  ###MLE of parameters under the null model	
  outnull = pmle.binom.sub(count, freq, size, m0, lambda = lambda, 
                           init.val, n.init, max.iter, n.iter, tol)
  alpha=outnull[[1]]
  theta=outnull[[2]]
  t0 = rbind(alpha, theta)
  
  ###Size of penalized function
  if (m0==1) { ah=c(0.5,0.5)
  if (is.null(CC)) { CC = 0.54 }
  }
  
  if (m0==2) {
    tb2 = tildeB22.binom(alpha, theta, size)
    ah  = c(0.5-acos(tb2[1,2])/2/pi, 0.5, acos(tb2[1,2])/2/pi)
    if (is.null(CC)) {
      temp = exp(5-10.6*tb2[1,2]-123/n)
      CC = 0.5*temp/(1+temp) 
    }
  }
  
  if (m0==3) {
    tb2 = tildeB22.binom(alpha, theta, size)
    a0 =0.5-acos(tb2[1,2])/4/pi-acos(tb2[1,3])/4/pi-acos(tb2[2,3])/4/pi
    a2 = 0.5-a0
    w123=(tb2[1,2]-tb2[1,3]*tb2[2,3])/sqrt(1-tb2[1,3]^2)/sqrt(1-tb2[2,3]^2)
    w132=(tb2[1,3]-tb2[1,2]*tb2[3,2])/sqrt(1-tb2[1,2]^2)/sqrt(1-tb2[3,2]^2)
    w231=(tb2[2,3]-tb2[2,1]*tb2[3,1])/sqrt(1-tb2[2,1]^2)/sqrt(1-tb2[3,1]^2)
    a1 =0.75-acos(w123)/4/pi-acos(w132)/4/pi-acos(w231)/4/pi
    a3 = 0.5-a1
    ah=c(a0,a1,a2,a3)
    if (is.null(CC)) {
      temp = exp(3.3-5.5*tb2[1,2]-5.5*tb2[2,3]-165/n)
      CC = 0.5*temp/(1+temp) 
    }
  }
  
  if (m0 > 3) {
    if (is.null(CC))  CC=0.5  ## exist no recommended value.
    tb2 = tildeB22.binom(alpha, theta, size)
    ah  = emtest.thm3(tb2, N=10000, tol=1e-8)
    ## use Monte Carlo.
  }
  
  out = emstat.binom(count, freq, outnull, size, m0, 
                     CC, n.init, max.iter, n.iter, tol, k)
  emnk  = out[[1]]
  alpha = out[[2]]
  theta = out[[3]]
  t1 = rbind(alpha, theta)
  p=sum(ah*pchisq(emnk,0:m0,lower.tail = F))
  
  if (rformat==F)
  {
    t0=rousignif(t0)
    t1=rousignif(t1)
    emnk=rousignif(emnk)
    p =rousignif(p)
    CC=rousignif(CC)
  }
  
  list('MLE under null order = m0'= t0,
       'pMLE under the alter order = 2m0' = t1,
       'EM-test statistic' = emnk,
       'P-value' = p,
       'Level of penalty'= CC,
       'iter.n'= out[[4]])
}