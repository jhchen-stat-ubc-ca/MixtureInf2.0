#' emtest.exp
#'
#' @description This function computes the EM-test statistic and the p-value for the null hypothesis H0:m=m0.
#' @param x The input data that can be either a vector or a matrix with the 1st column being the observed values
#' and the 2nd column being the corresponding frequency. 
#' @param m0 The order under the null hypothesis.
#' @param CC A optional tuning parameter for the EM-test procedure.
#' @param init.val The initial values chosen for the EM-algorithm in order to compute the PMLE under the null model.
#' @param n.init The number of initial values chosen for the EM-algorithm.
#' @param n.iter The least amount of iterations for all initial values in the EM-algorithm.
#' @param tol The tolerance value for the convergence of the EM-algorithm.
#' @param k The number of EM iterations needed to obtain the EM-test statistic.
#' @param max.iter The maximum number of iterations.
#' @param rformat The format of the output, rformat=T means the format is determined by R.
#'		rformat=F means the format is determined by our default setting. When the output is
#'		larger than 0.001, it is determined by round(output,3); When the output is less than 0.001,
#'		it is determined by signif(output,3).
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
#' @examples n = 4000
#' mu = c(3, 9, 18)
#' alpha = c(.5, .3, .2)
#' x = rmix.exp(n, alpha, mu) 
#' n.init = 10; n.iter =50
#' max.iter = 2000
#' emtest.exp(x, m0 = 2, CC = NULL, init.val = NULL, n.init=10, n.iter=50, tol = 1e-6, k=3, rformat = FALSE)
emtest.exp <- function(x, m0 = 1, CC = NULL, init.val = NULL,
                       n.init=10, n.iter=50, tol = 1e-6, k=3, max.iter = 2000, 
                       rformat = FALSE) {
  if (is.data.frame(x)) stop("Data must be a vector or a matrix")
  
  if (is.matrix(x)) {
    xx=c()
    for (i in 1:nrow(x)) xx=c(xx,rep(x[i,1],x[i,2]))
    x=xx
  }
  
  ###MLE of parameters under the null model	
  outnull = pmle.exp.sub(x, m0, 0, init.val, n.init, max.iter, n.iter,tol)
  alpha = outnull$alpha
  theta = outnull$theta
  t0    = rbind(alpha, theta)
  
  n = length(x)
  ### Decide the size of penalized function
  if (is.null(CC)) {
    if (m0==1) { 
      CC=exp(0.74+82/n)/(1+exp(0.74+82/n))
      ah = c(0.5, 0.5)
    }
    
    if (m0==2) {
      tb2 = tildeB22.exp(alpha,theta)
      temp = exp(2.3-8.5*tb2[1,2])
      CC = temp/(1 + temp)
      ah = c(0.5-acos(tb2[1,2])/2/pi, 0.5, acos(tb2[1,2])/2/pi)
    }
    
    if (m0==3) {
      tb2 = tildeB22.exp(alpha,theta)
      temp = exp(2.2-5.9*tb2[1,2]-5.9*tb2[2,3])
      CC = temp/(1 + temp)
      
      a0=0.5-acos(tb2[1,2])/4/pi-acos(tb2[1,3])/4/pi-acos(tb2[2,3])/4/pi
      a2=0.5-a0
      w123=(tb2[1,2]-tb2[1,3]*tb2[2,3])/sqrt(1-tb2[1,3]^2)/sqrt(1-tb2[2,3]^2)
      w132=(tb2[1,3]-tb2[1,2]*tb2[3,2])/sqrt(1-tb2[1,2]^2)/sqrt(1-tb2[3,2]^2)
      w231=(tb2[2,3]-tb2[2,1]*tb2[3,1])/sqrt(1-tb2[2,1]^2)/sqrt(1-tb2[3,1]^2)
      a1=0.75-acos(w123)/4/pi-acos(w132)/4/pi-acos(w231)/4/pi
      a3=0.5-a1
      ah=c(a0,a1,a2,a3)
    }
    
    if (m0 > 3) {
      CC=0.5
      tb2 = tildeB22.exp(alpha,theta)
      ah = emtest.thm3(tb2, N=10000, tol=1e-8)
    }
  }
  
  out = emstat.exp(x, outnull, CC, n.init, n.iter, tol, k, max.iter)
  emnk=out[1]
  alpha=out[2:(2*m0+1)]
  theta=out[(2*m0+2):(4*m0+1)]
  t1=rbind(alpha,theta)
  p =sum(ah*pchisq(emnk,0:m0,lower.tail = F))
  
  if (rformat==F) 	{
    t0=rousignif(t0)
    t1=rousignif(t1)
    emnk=rousignif(emnk)
    p=rousignif(p)
    CC=rousignif(CC)
  }
  
  list('MLE of parameters under null hypothesis (order = m0)'= t0,
       'Parameter estimates under the order = 2m0' = t1,
       'EM-test Statistic' = emnk,
       'P-value'= p,
       'Level of penalty'= CC)
}