#' emtest.binom
#'
#' @description This function computes the EM-test statistic and the p-value for H0: m=m0.
#' @param x The input data that can either be a vector or a matrix with two columns: count and freq. 
#' @param size The number of trials of the binomial.
#' @param m0 The order under the null hypothesis.
#' @param CC A optional tuning parameter for the EM-test procedure. 
#' @param lambda A tuning parameter for the fitted null model.
#' @param init.val The initial values for the PMLE under the null model. 
#' @param n.init The number of initial values chosen for the EM-algorithm.
#' @param n.iter The least amount of iterations for all initial values.
#' @param max.iter The maximum number of iterations allowed.
#' @param tol The tolerance value for the convergence of the EM-algorithm.
#' @param k The number of EM iterations in order to obtain the EM-test statistic.
#' @param rformat The format of the output. See rforma().
#'
#' @return Return an object of class EM-test with the following elements:
#' The MLE of the parameters under the null hypothesis (order = m0)
#' The PMLE of the parameters under the specific alternative hypothesis whose order is 2m0
#' EM-test statistic
#' P-value
#' Level of penalty
#' The number of iterations
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#' @references Chen, J. and Li, P. (2011). Tuning the EM-test for the order of finite mixture models. The Canadian Journal of Statistics. 39, 389-404.
#' Li, P. and Chen, J. (2010). Testing the order of a finite mixture model. JASA. 105, 1084-1092.
#' Li, P., Chen, J. and Marriott, P. (2009). Non-finite Fisher information and homogeneity: The EM approach. Biometrika. 96, 411-426.

#'
#' @examples theta = c(.3, .3, .6)
#' alpha = c(.5, .3, .2)
#' size = 20; n = 1000
#' xx=rmix.binom(n, size, alpha, theta)
#' y=as.matrix(table(xx))
#' count = as.numeric(rownames(y))
#' freq=y[,1]
#' init.val = c(.2, .7, .2, .8)
#' emtest.binom(cbind(count, freq), size, m0=2, CC=NULL, lambda = 0, init.val,
#'              n.init=10, n.iter=50, max.iter = 2000, tol=1e-6, k=5, rformat=FALSE)
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