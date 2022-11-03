#' pmle.pois
#'
#' @description This function computes the PMLE of parameters under a Poisson mixture
#' @param x Data that can be either a vector or a matrix of observations. 
#'          1st column: counts 
#'          2nd column: frequencies.
#' @param m0 The order of the finite mixture model.
#' @param lambda The size of the penalty function of the mixing proportions.
#' @param init.val A user provided initial values for the EM-algorithm to 
#'                 compute the PMLE under the null model.
#' @param n.init The number of initial values chosen for the EM-algorithm.
#' @param n.iter Least amount of iterations for all initial values in the EM-algorithm.
#' @param max.iter Maximum amount of iterations.
#' @param rformat A specific format, please see rousignif.R function. 
#' @return
#' @export
#'
#' @examples n = 1000
#' mu = c(9, 10)
#' alpha = c(.4, .6)
#' x = rmix.pois(200, alpha, mu)
#' pmle.pois (x, m0=2, lambda = 1, init.val = NULL, n.init=10, n.iter=50, tol=1e-6, max.iter=5000, rformat=FALSE)
pmle.pois <- function(x, m0, lambda = 1, init.val = NULL, n.init=10, 
                      n.iter=50, tol=1e-6, max.iter=5000, rformat=FALSE)
  #init,val:	a matrix with two rows: 
  #      mixing proportion, subpopulation means.
  #lambda:	size of penalty function of mixing proportions.
  #n.init: 	number of initial values chosen for the EM-algorithm.
  #n.iter:    the number of iterations for initial values.
  #tol: 	tolerance value for the convergence of the EM-algorithm.
{
  if(is.matrix(init.val)) {  dim.init = dim(init.val)
  if( dim.init[1]!=2 | dim.init[2]!=m0) 
    stop("dimension of initial value shoudl be 2, m0")
  }
  if(is.vector(x)) {
    min.x = min(x); max.x = max(x)
    count = min.x:max.x
    freq = count*0
    for(i in count) freq[i- min.x + 1]= sum(x==i)
    xx = cbind(count, freq)
  }  ## unify the format for subsequent calculation.
  
  if(is.matrix(x)) {
    if(dim(x)[2]!=2) stop("data matrix should have exactly 2 columns")
    xx = x
  }
  out = pmle.pois.sub(xx, m0, lambda, init.val, n.init, 
                      n.iter, tol, max.iter)
  alpha=out$alpha
  theta=out$theta
  loglik=out$loglik
  ploglik=out$ploglik
  iter.n = out$iter.n
  
  if (rformat==F) {
    alpha=rousignif(alpha)
    theta=rousignif(theta)
    loglik=rousignif(loglik)
    ploglik=rousignif(ploglik) }
  
  if (lambda==0 | m0==1) {
    list("MLE of mixing proportions:"=alpha,
         "MLE of component parameters:"=theta,
         "log-likelihood:" = loglik, 
         "Penalized log-likelihood:"= loglik,
         "iter.n:" = out$iter.n)
  } else {
    list(
      "PMLE of mixing proportions:"= alpha,
      "PMLE of component parameters:"= theta,
      "log-likelihood:"= loglik,
      "Penalized log-likelihood:"= ploglik,
      "iter.n:" = out$iter.n)
  }
}