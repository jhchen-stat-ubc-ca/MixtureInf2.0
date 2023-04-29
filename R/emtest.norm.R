#' emtest.norm
#'
#' @param x The input data that can either be a vector or a matrix with the 1st column being the observed values
#'          and the 2nd column being the corresponding frequency.
#' @param m0 The order of the finite normal mixture model under the null hypothesis. 
#' @param pens A 2-dimensions vector being the size of the penalized functions for 
#'             the mixing proportion and the variance.
#' @param inival The initial value chosen for the EM-algorithm to compute the PMLE under the null model.
#' @param len The number of initial values chosen for the EM-algorithm.
#' @param niter The least amount of iterations for all initial values in the EM-algorithm. 
#' @param tol The tolerance value for the convergence of the EM-algorithm.
#' @param k The amount of EM iterations in order to obtain the EM-test statistic. 
#' @param rformat F means the format of the output is determined by our default setting. When the output is
#'                larger than 0.001, it is determined by round(output,3); When the output is less than 0.001,
#'                it is determined by signif(output,3).
#'
#' @return Return an object of class EM-test with the following elements:
#' The MLE of the parameters under the null hypothesis (order = m0)
#' The PMLE of the parameters under the specific alternative hypothesis whose order is 2m0
#' EM-test statistic
#' P-value
#' Level of penalty
#' The number of iterations
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#' @references Chen, J. and Li, P. (2009). Hypothesis test for normal mixture models: The EM approach. The Annals of Statistics. 37, 2523-2542.
#' Chen, J., Li, P. and Fu, Y. (2012). Inference on the order of a normal mixture. JASA. 107, 1096-1105.
#'
#' @examples data(grains)
#' out <- emtest.norm(unlist(grains),m0 = 2)
#' plotmix.norm(unlist(grains), alpha = out[[1]][1,], mu = out[[1]][2,], sigma = out[[1]][3,], m0 = 2)
emtest.norm <-
  function(x, m0=1, pen.size=NULL, init.val=NULL, n.init = 10, 
           n.iter=50, tol=1e-6, k=3, rformat=FALSE) {
    xx = x
    if (is.matrix(x)) {
      xx=c()
      for (i in 1:nrow(x)) xx=c(xx, rep(x[i,1], x[i,2]))   }
    nn = length(xx)
    
    ## PMLE and log-likelihood value under the null model	
    para0 = pmle.norm(x, m0, lambda=0, an = NULL, init.val, 
                      n.init, n.iter, tol)
    ln0 = para0[[2]]     ## likelihood under null, penalty excluded.
    sigma0 = para0[[1]][3,]
    
    ###Size of penalties: revised.
    if (is.null(pen.size)) {
      pen.size = c()
      pen.size[1] = 1
      pen.size[2] = EMtest.norm.penalty(m0, para0, nn) }
    
    ###create beta vectors for EM-test
    bbeta=c()
    for(h in 1:m0) {
      bbeta = rbind(cbind(bbeta,rep(0.1,3^{h-1})),
                    cbind(bbeta,rep(0.3,3^{h-1})),
                    cbind(bbeta,rep(0.5,3^{h-1})))
    }
    
    output=c()
    for (i in 1:(3^m0)) {
      para1 = EMtest.maxmm.norm(xx, bbeta[i,], m0, an = pen.size[2], 
                                para0, n.init, n.iter, tol)
      out = EMtest.norm.iter(xx, nn, m0, para1, sigma0, bbeta[i,], pen.size, k)
      output = rbind(output, out)  }
    
    index=which.max(output[,(6*m0+1)])
    out.para= output[index,1:(6*m0)]
    ln1 = output[index,(6*m0+1)]
    
    ###EM-test statistic
    EM.stat = 2*(ln1 - ln0)
    names(EM.stat) <- NULL
    p.value = pchisq(EM.stat, (2*m0), lower.tail=F)	
    
    ###output
    null.para = para0[[1]]
    full.para = rbind(out.para[1:(2*m0)], out.para[(2*m0+1):(4*m0)], 
                      out.para[(4*m0+1):(6*m0)])
    
    if (rformat==F) {
      null.para = rousignif(null.para)
      full.para = rousignif(full.para)
      EM.stat = rousignif(EM.stat)
      p.value = rousignif(p.value)		
      pen.size= rousignif(pen.size)	 }
    
    list(
      'MLE under null hypothesis (order = m0)'= null.para,
      'Parameter estimates under the order = 2m0'= full.para,
      'EM-test Statistics'= EM.stat,
      'P-values'= p.value,
      'Level of penalty'= pen.size)
  }