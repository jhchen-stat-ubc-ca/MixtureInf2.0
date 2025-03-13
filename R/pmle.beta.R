#' pmle.beta
#'
#' @description This function computes the PMLE of the parameters under a mixture of betas.
#' @param x A vector of the observed values.
#' @param m0 The order of the finite mixture model.
#' @param n.iter The number of EM iterations for each initial value. The one that gained the most in likelihood will be iterative further. 
#' @param max.iter The maximum number of iterations for the EM algorithm. 
#' @param tol The tolerance value for the convergence of the EM-algorithm, default value: tol = 1e-6.
#' @param epsilon The size of the penalized function of the mixing distribution, default value: epsilon = 1.
#' @param an A size control parameter that controls the severity of the penalty. The recommended value is n^{-3/2}.
#' @param seed Used for reproducing the same sequence of output numbers. Default is NULL.
#' @param maxit The maximum number of iterations for the Broyden or Newton methods used in the method of moments estimation.
#' 
#' @return  The PMLE of the parameters with order = m0 (mixing proportions, mixing alphas
#'          and mixing betas), log-likelihood value at the PMLE and the penalized log-likelihood
#'          value at the PMLE, classification on which subpopulation does the observed value belongs to. 
#' 
#'
#' @examples  
#' data <- c(rbeta(50,10,2),rbeta(50,3,18))
#' pmle.beta(data,2)
#' @export
#' 
#' 
pmle.beta <- function(x, m0, n.iter=10, max.iter=5000, tol=1e-6, epsilon=1, 
                      an=NULL, seed=NULL, maxit=5000) {
  if(is.null(an)) {an = length(x)^(3/2)}
  output = c()
  init_params=MoM_BMM(x,m0,seed,maxit)
  uniq_init_params=unique(do.call(rbind, lapply(1:ncol(init_params), function(i) init_params[,i]$x)))
  for (i in 1:nrow(uniq_init_params)) {
    para0 = c(uniq_init_params[i,1:(m0-1)],1-sum(uniq_init_params[i,1:(m0-1)]),uniq_init_params[i,m0:(3*m0-1)])
    for (j in 1:n.iter) {
      outpara = pmle.beta.sub(x,m0,para0,an,epsilon)
      para0 = outpara[1:(3*m0)]
    }
    output = rbind(output, outpara)
  }
  index = which.max(output[,(3*m0+2)])
  para0 = output[index,1:(3*m0)]
  ploglike0 = output[index,(3*m0+2)]
  increment = 1
  tt = 0
  
  while (increment > tol && tt < max.iter) {
    outpara = pmle.beta.sub(x,m0,para0,an,epsilon)
    para0 = outpara[1:(3*m0)]
    ploglike1 = outpara[3*m0+2]
    increment = ploglike1 - ploglike0
    ploglike0 = ploglike1
    tt = tt+1
  }
  mix_porp = para0[1:m0]
  alpha = para0[(m0+1):(2*m0)]
  beta = para0[(2*m0+1):(3*m0)]
  pdf.sub = t(mapply(function(pp, aa, bb) pp * dbeta(x, aa, bb), mix_porp, alpha, beta))
  pdf.mixture = apply(pdf.sub,2,sum) +1e-100
  ww = sweep(pdf.sub, 2, pdf.mixture, FUN = "/")
  return(list(mix_porp= rousignif(mix_porp),
              alpha = rousignif(alpha),
              beta= rousignif(beta),
              loglik=rousignif(outpara[3*m0+1]),
              ploglik=rousignif(outpara[3*m0+2]),
              iter.n=tt,
              classification=apply(t(ww), 1, which.max)))
}