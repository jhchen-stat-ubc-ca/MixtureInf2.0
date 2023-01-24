#' pmle.binom
#'
#' @description This function computes the PMLE under a binomial mixture
#' @param x The input data, either a vector or a matrix columns with two columns: count and freq.
#' @param size The number of trials of the subpopulation binomial distribution.
#' @param m0 The order of the finite mixture model to be fitted.
#' @param lambda The size of the penalty function of the mixing proportions.
#' @param init.val The initial values for the EM-algorithm. Usually a matrix with 2 rows: mixing proportions and probs of success.
#' @param n.init The number of the initial values for the EM-algorithm.
#' @param n.iter The number of the initial EM iterations for all initial values.
#' @param max.iter The maximum amount of iterations allowed. Set to 5000.
#' @param tol The tolerance value for the convergence of the EM-algorithm.
#' @param rformat The format of the output
#'
#' @return
#' @export
#'
#' @examples > size = 25
#' x = rmix.binom(1000, size, alpha, theta)
#' pmle.binom(x, size, m0=1, lambda=1, init.val = NULL, n.init=10,n.iter = 50, max.iter = 5000, tol = 1e-6, rformat=FALSE)
pmle.binom <- function(x, size, m0=1, lambda=1, init.val = NULL, n.init=10,
                       n.iter = 50, max.iter = 5000, tol = 1e-6, rformat=FALSE) {
  if(is.vector(x)) {
    y=as.matrix(table(x))
    count = as.numeric(rownames(y))
    freq=y[,1]
  }
  if(is.matrix(x)) {
    count= x[,1]
    freq =x[,2]
  }   ## handle both types of data.
  
  out = pmle.binom.sub(count, freq, size, m0, lambda, init.val, 
                       n.init, n.iter, max.iter, tol)
  alpha=out$alpha
  theta=out$theta
  loglik=out$loglik
  ploglik=out$ploglik
  
  if (rformat==F) {
    alpha=rousignif(alpha)
    theta=rousignif(theta)
    loglik=rousignif(loglik)
    ploglik=rousignif(ploglik)
  }
  
  if (lambda==0 | m0==1)
    list('MLE of mixing proportions:'=alpha,
         'MLE of component parameters:'=theta,
         'log-likelihood:'=loglik,
         'number of additional EM iterations:'= out$iter.n)
  else
    list('PMLE of mixing proportions:'=alpha,
         'PMLE of component parameters:'=theta,
         'log-likelihood:'=loglik,
         'Penalized log-likelihood:'=ploglik,
         'number of additional EM iterations:'= out$iter.n)
}