#' pmle.norm
#'
#' @description This function computes the PMLE or MLE of the parameters under a mixture of normals with unequal variance.
#' @param x The data that can either be a vector or a matrix with the 1st column being the observed values
#' and the 2nd column being the corresponding frequencies.
#' @param m0 The order of the finite mixture model, default value: m0 = 1.
#' @param lambda The size of the penalized function of the mixing distribution, default value: lambda = 1.
#' @param an
#' @param init.val NULL or a 3 X m0 matrix with rows made of mixing probability, component means, and variances.
#' @param n.init A computer generated n.init initials value.
#' @param n.iter The number of EM iterations for each initial values. The one gained the most in likelihood will be iterative further. 
#' @param tol The tolerance value for the convergence of the EM-algorithm, default value: tol = 1e-6.
#' @param rformat The format of the output. If rformat=T, 
#' it means the output format is determined by R software. 
#' If rformat=F, it means the output format is determined by our default setting. 
#' When the output is larger than 0.001, it is determined by round(output,3); 
#' When the output is less than 0.001, it is determined by signif(output,3).
#' The default value of rformat is F.
#'
#' @return  The PMLE or MLE of the parameters with order = m0 (mixing proportions, mixing means
#           and mixing variances), log-likelihood value at the PMLE or MLE and the penalized log-likelihood
#           value at the PMLE.

#' @export 
#'
#' @examples  n=100
#' alpha=c(.2,.8)
#' mu=c(8,9)
#' x=rmix.norm(n,alpha,mu)
#' pmle.norm(x, m0=2, lambda = 1, an = NULL, init.val = NULL, n.init = 10, 
#' n.iter=50, max.iter = 5000,tol=1e-8, rformat = FALSE)
pmle.norm <- function(x, m0, lambda = 1, an = NULL, init.val = NULL,
                      n.init = 10, n.iter=50, max.iter = 5000, tol=1e-8, rformat = F) {
  if(m0==1) stop("You do not need this function for MLE")
  
  if(is.data.frame(x)) stop("data format must be vector or matrix")
  
  if (is.vector(x)) xx = x
  ## plain vector case
  
  if (is.matrix(x)) {
    xx=c()
    for (i in 1:nrow(x)) xx=c(xx, rep(x[i,1],x[i,2]))
  }       ## it translates the frequency into a plain vector.
  
  if(is.null(an)) an = length(xx)^(-0.5)
  ## If not given, use the default penalty n^{-1/2}
  
  if(lambda < tol) lambda = tol
  ## avoid 0 mixing probability.
  
  out = pmle.norm.sub(xx, m0, lambda, an, init.val, n.init, 
                      n.iter, max.iter, tol)
  ## leave the computation to the other function.
  alpha.mix = out[[1]]
  mean.mix = out[[2]]
  var.mix = out[[3]]
  loglik = out[[4]]
  ploglik = out[[5]]
  
  if (!rformat) {
    alpha.mix = rousignif(alpha.mix)
    mean.mix = rousignif(mean.mix)
    var.mix = rousignif(var.mix)
    loglik = rousignif(loglik)
    ploglik = rousignif(ploglik)
  }
  
  list('PMLE of mixing distribution' = 
         rbind(alpha.mix, mean.mix, var.mix),
       'log-likelihood' = loglik,
       'Penalized log-likelihood'= ploglik,
       'iter.n' = out[[6]] )
}