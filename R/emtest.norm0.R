#' emtest.norm0
#'
#' @description This function computes the EM-test statistic and the p-value for the hypothesis H0:m=m0. 
#' @param x The input data that can be either a vector or a matrix with the 1st column being the observed values
#'          and the 2nd column being the corresponding frequency.
#' @param var The known component variance.
#' @param m0 The order under the null hypothesis.
#' @param init.val The initial values chosen for the EM-algorithm to compute the PMLE under the null model
#' @param n.init The amount of initial values chosen for the EM-algorithm.
#' @param n.iter The least number of iterations for all initial values in the EM-algorithm. 
#' @param max.iter The maximum number of iterations.
#' @param tol The tolerance value for the convergence of the EM-algorithm.
#' @param k The amount of EM iterations in order to obtain the EM-test statistic.
#' @param rformat F means the format of output is determined by our default setting. When the output is
#'                larger than 0.001, it is determined by round(output,3); When the output is less than 0.001,
#'                it is determined by signif(output,3).
#'
#' @return
#' @export
#'
#' @examples
emtest.norm0 <-
  function(x, var.sub=1, m0=1, init.val=NULL, n.init=10, n.iter=50, 
           max.iter = 5000, tol = 1e-6, k=3, rformat=FALSE)  {
    #x: 		data, can be either a vector or a matrix 
    #         1st column being the observed values 
    #       	2nd column being the corresponding frequencies. 
    #var.sub:		known component variance.
    #m0: 		order of the null mixture.
    #lambda: 		optional tuning parameter for EM-test procedure.
    #init.val:	initial values chosen for the EM-algorithm to 
    #           compute the PMLE0 under the null model.
    #n.iter:  number of iterations for each initial value
    #n.init: 	number of initial values for the EM-algorithm.
    #tol: 	tolerance value for the convergence of the EM-algorithm.
    #k:  number of EM iterations in the EM-test statistic definition.
    #rformat	See iformat()
    
    nn = length(x)
    xx = x
    if (is.matrix(x)) {
      xx=c()
      for (i in 1:nrow(x)) xx=c(xx,rep(x[i,1],x[i,2]))
    }
    xx = xx/sqrt(var.sub)
    ### data is a plain vector and scaled.
    
    if(m0 > 1) {
      null.output = pmle.norm0(xx, 1, m0, lambda=0, init.val, 
                               n.init, n.iter, max.iter, tol)
      
      ### no penalty at null
      mle0 = rbind(null.output[[1]], null.output[[2]])
      ln0  = null.output[[3]]    ### no penalty.
      temp = tildeB22.norm0(mle0[1,], mle0[2,])
      degenerate = temp[[2]] 
      tb2 = temp[[1]]
    } else {
      mle0 = rbind(1, mean(xx))
      ln0 = sum(log(dnorm(xx, mean(xx), var.sub^.5 )))
      tb2 = 1
    }
    ## fit null model, evaluate tildeB22 for p-value
    ##    calculation later.
    if(degenerate) {
      list('MLE under null(order = m0)'= mle0,
           'MLE under order = 2*m0' = F,
           'EM-test statistic' = 0,
           'P-value'= 1,
           'Level of penalty'= NULL,
           'degenerate fitted null'= T)
    } else {
      ### skip this step in this case
      if (m0==1)  CC = 0.54
      if (m0==2) {
        C1 = 0.5*exp(5-10.6*tb2[1,2]-123/nn)
        CC = C1/(1 + 2*C1)    
      }
      if (m0==3) {
        CC = 0.5*exp(3.3-5.5*tb2[1,2]-5.5*tb2[2,3]-165/nn)
        CC = CC/(1 + 2*CC)   
      }
      if (m0 > 3) CC = 0.5
      ### recommended C value in EM-test Chen & Li (2010JASA pp1087)
      
      EM.stat = emstat.norm0(xx, mle0, ln0, m0, CC, n.init, n.iter, tol, k)
      
      emnk= EM.stat[1]
      alpha=EM.stat[2:(2*m0+1)]
      theta=EM.stat[(2*m0+2):(4*m0+1)]*sqrt(var)
      pmle1 = rbind(alpha, theta)
      ###  emstat.norm0 needs double check.
      
      if (m0==1) ah=c(0.5, 0.5)
      if (m0==2) {
        tt = acos(tb2[1,2])/2/pi
        ah = c(0.5 - tt, 0.5, tt)
      }
      if (m0==3) {
        a0 = 0.5-(acos(tb2[1,2])+acos(tb2[1,3])+acos(tb2[2,3]))/(4*pi)
        a2 = 0.5-a0
        w123=(tb2[1,2]-tb2[1,3]*tb2[2,3])/sqrt(1-tb2[1,3]^2)/sqrt(1-tb2[2,3]^2)
        w132=(tb2[1,3]-tb2[1,2]*tb2[3,2])/sqrt(1-tb2[1,2]^2)/sqrt(1-tb2[3,2]^2)
        w231=(tb2[2,3]-tb2[2,1]*tb2[3,1])/sqrt(1-tb2[2,1]^2)/sqrt(1-tb2[3,1]^2)
        a1 = 0.75 - (acos(w123) + acos(w132) + acos(w231))/(4*pi)
        a3 = 0.5 - a1
        ah = c(a0, a1, a2, a3)
      }
      
      if (m0 > 3) ah = emtest.thm3(tb2, N=10000, tol=1e-8)
      ### Numerical evaluation
      
      pp = sum(ah*pchisq(emnk, 0:m0, lower.tail = F))
      
      mle0[2,] = mle0[2,]*sqrt(var.sub)
      pmle1[2,] = pmle1[2,]*sqrt(var.sub)
      if (rformat==F) {
        mle0 = rousignif(mle0)
        pmle1 = rousignif(pmle1)
        emnk = rousignif(emnk)
        pp = rousignif(pp)
        CC = rousignif(CC)
      }
      
      list('MLE under null(order = m0)'= mle0,
           'MLE under order = 2*m0' = pmle1,
           'EM-test statistic' = emnk,
           'P-value'= pp,
           'Level of penalty'= CC,
           'degenerate null'= degenerate)
    }
  }