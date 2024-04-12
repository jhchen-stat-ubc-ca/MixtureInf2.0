#' pmle.norm0
#'
#' @description This function computes the PMLE or MLE of the parameters under a mixture of normals with equal and known
#' variance. When the level of penalty is 0, PMLE reduces to MLE.
#' @param x The input data that can either be a vector or a matrix with the 1st column being the observed values
#' and the 2nd column being the corresponding frequencies.
#' @param var.sub A known subpopulation variance, default at 1. 
#' @param m0 The order of the finite mixture model.
#' @param lambda The size of the penalty function of the mixing proportions.
#' @param ini.val The initial values chosen for the EM-algorithm, a 3m0-dimension vector including m0 mixing proportions, 
#' m0 component means and m0 component variances, or a matrix with 3m0 columns, default value: inival = NULL. (if not provided, random initial values are used.)
#' @param n.init The	number of initial values for the EM-algorithm.
#' @param n.iter The number of EM iterations for all initial values.
#' @param max.iter  Maximum amount of EM iterations, it stops at 5000.
#' @param tol The tolerance value for the convergence of the EM-algorithm, default value: tol = 1e-8.
#' @param iter.n This reports the number of EM-iteration employed.
#' @param rformat The format of the output. If rformat=T, 
#' it means the output format is determined by R software. 
#' If rformat=F, it means the output format is determined by our default setting. 
#' When the output is larger than 0.001, it is determined by round(output,3); 
#' When the output is less than 0.001, it is determined by signif(output,3).
#' The default value of rformat is F.
#'
#' @return It returns the PMLE or MLE of the parameters with order = m0 (mixing proportions and component parameters), 
#' log-likelihood value at the PMLE or MLE and the penalized log-likelihood value at the PMLE.
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples x <- rmix.norm(200,c(0.3,0.7),c(-1,2))
#' pmle.norm0(x,m0=2)
#' @export
pmle.norm0 <- function(x, var.sub = 1, m0, lambda = 1, 
                       init.val = NULL, n.init = 10, n.iter=50, 
                       max.iter = 5000, tol=1e-8, rformat = F) {
  if(m0==1) stop("You do not need this function for MLE")
  
  if (is.vector(x)) xx = x     ## plain vector
  if (is.matrix(x)) { xx=c()
  for (i in 1:nrow(x)) xx=c(xx, rep(x[i,1],x[i,2]))  } 
  xx = xx/sqrt(var.sub);   
  ###standardize data with known var 
  ## transform the frequency into a plain vector.
  nn = length(xx)
  min.x = min(xx); max.x = max(xx)
  
  if (is.null(init.val)) {
    init.val = list()
    for(i in 1:n.init) {
      alpha = runif(m0)
      alpha = alpha/sum(alpha)
      theta = runif(m0, min.x, max.x)
      init.val[[i]] = rbind(alpha, theta)
    }	
  } else { n.init = 1; 
  temp = init.val
  init.val = list()
  init.val[[1]]= temp
  init.val[[1]][2,] = init.val[[1]][2,]/sqrt(var.sub)
  }
  
  out = pmle.norm0.sub(xx, m0, lambda, 
                       init.val, n.init, n.iter, max.iter, tol)
  
  alpha = out$alpha
  theta = out$theta*sqrt(var.sub)
  loglik = out$loglik-nn/2*log(var.sub)
  ploglik = out$ploglik-nn/2*log(var.sub)
  iter.n = out$iter.n
  
  if (rformat==F) {
    alpha=rousignif(alpha)
    theta=rousignif(theta)
    loglik=rousignif(loglik)
    ploglik=rousignif(ploglik)
  }
  
  list('PMLE of mixing proportions:'= alpha,
       'PMLE of subpop means:'= theta,
       'log-likelihood:'= loglik,
       'Penalized log-likelihood:'= ploglik,
       'iter.n'= iter.n)
}