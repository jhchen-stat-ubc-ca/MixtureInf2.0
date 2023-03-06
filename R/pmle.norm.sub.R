#' pmle.norm.sub
#'
#' @description It is used in the pmle.norm function, it does the actual EM-algorithm and computes the penalized MLE. 
pmle.norm.sub <- function(x, m0, lambda, an, init.val, n.init, 
                          n.iter, max.iter, tol) {
  sample.var = var(x);  n=length(x)
  if (!is.null(init.val)) n.init = 1
  
  output=c()
  for (i in 1:n.init) {
    if (is.null(init.val)) {
      ## generate random initial values if not provided.
      mu = sort(sample(x, m0))
      tmp = (mu[-1] + mu[-m0])/2
      alpha = rep(1, m0)
      for(ii in 1:(m0-1)) alpha[ii] = sum(x<tmp[ii])/n
      alpha = alpha - c(0, alpha[-m0])
      sigma = c(mu, max(x))-c(min(x), mu)
      sigma = (sigma[-1]+sigma[-m0])/4
      sigma = sigma^2
    }	else {	### if initial mixing distribution is provided
      alpha = init.val[1,]
      mu = init.val[2,]
      sigma = init.val[3,]
    }
    para0 = c(alpha, mu, sigma)  
    ## initial value created.
    for (j in 1:n.iter) {
      outpara = pmle.norm.sub.a(x, sample.var, m0, para0, lambda, an)
      para0 = outpara[1:(3*m0)]
    }  ###run n.iter EM-iterations first
    output = rbind(output, outpara[1:(3*m0+2)])	
  }
  ## pick up the winner
  index = which.max(output[,(3*m0+2)])
  para0 = output[index,1:(3*m0)]
  ploglike0 = output[index,(3*m0+2)]
  increment = 1
  tt = 0
  ### restart EM-iteration with the winner.
  while (increment > tol & tt < max.iter) {
    outpara = pmle.norm.sub.a(x, sample.var, m0, para0, lambda, an)
    para0 = outpara[1:(3*m0)]
    ploglike1 = outpara[3*m0+2]
    increment = ploglike1 - ploglike0
    ploglike0 = ploglike1
    tt = tt+1
  }
  
  list('alpha'= outpara[1:m0],
       'means' = outpara[(m0+1):(2*m0)],
       'variances'= outpara[(2*m0+1):(3*m0)],
       'loglik'=outpara[3*m0+1],
       'ploglik'=outpara[3*m0+2],
       'iter.n'=tt)
}