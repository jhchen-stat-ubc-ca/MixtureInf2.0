#' pmle.pois.sub
#'
#' @description This function is used by pmle.pois. It does the real work but not stand alone.
pmle.pois.sub <- function(x, m0, lambda, init.val, n.init, n.iter, 
                          tol, max.iter=5000) {
  count = as.numeric(x[,1])  
  freq  = as.numeric(x[,2])
  min.x = min(count);  max.x = max(count); nn = length(count)
  
  pmf.sub = ww = matrix(0, m0, nn)
  output = NULL
  
  if(m0 ==1) stop("This program is not for homogenous Poisson(m0=1)")
  
  if(is.null(init.val))  {
    ## generate initial values
    for(i in 1:(n.init)) {
      alpha.i = runif(m0)
      alpha.i = alpha.i/sum(alpha.i)
      theta.i = sort(sample(count, m0, prob = freq))
      ### randomly generate one initial value.
      
      for (ii in 1:n.iter)    {  
        ### run n.iter EM-iterations for each init.val
        for(j in 1:m0) pmf.sub[j,] = dpois(count, theta.i[j])*alpha.i[j]
        pmf     = colSums(pmf.sub)+1e-50
        for(j in 1:m0) {
          ww[j,] = pmf.sub[j,]/pmf
          alpha.i[j] = (sum(freq*ww[j,])+lambda)/(sum(freq)+m0*lambda)
          theta.i[j] = sum(freq*count*ww[j,])/sum(freq*ww[j,])
        }
      }
      ### n.iter iterations completed with the (i)th initial.
      for(j in 1:m0) pmf.sub[j,] = dpois(count, theta.i[j])*alpha.i[j]
      pmf     = colSums(pmf.sub)+1e-50
      pln = sum(freq*log(pmf)) + lambda*sum(log(alpha.i))
      output = rbind(output, c(alpha.i, theta.i, pln))
      ### output of the (i)th initial value.
    }
    index = which.max(output[,(2*m0+1)])
    alpha = output[index,1:m0]
    theta = output[index,(m0+1):(2*m0)]
    pln0 = output[index,(2*m0+1)]
    ### pick the winner.
  } else {
    if(length(init.val)!=2*m0) stop("init.val error")
    alpha = init.val[1:m0]; alpha = alpha/sum(alpha)
    theta = init.val[(m0+1):(2*m0)]  
    for(j in 1:m0) pmf.sub[j,] = dpois(count, theta[j])*alpha[j]
    pmf     = colSums(pmf.sub)+1e-50
    pln0 = sum(freq*log(pmf)) + lambda*sum(log(alpha))
  }   ### if user has supplied an initial value.
  
  err=1;  tt=0
  for(j in 1:m0) pmf.sub[j,] = dpois(count, theta[j])*alpha[j]
  pmf     = colSums(pmf.sub)+1e-50
  
  while(err > tol & tt < max.iter) {
    ###EM-iteration with the winning initial value
    for(j in 1:m0) {
      ww[j,] = pmf.sub[j,]/pmf
      alpha[j] = (sum(freq*ww[j,])+lambda)/(sum(freq)+m0*lambda)
      theta[j] = sum(freq*count*ww[j,])/sum(freq*ww[j,])
    }
    
    for(j in 1:m0) pmf.sub[j,] = dpois(count, theta[j])*alpha[j]
    pmf  = colSums(pmf.sub)+1e-50
    pln1 = sum(freq*log(pmf)) + lambda*sum(log(alpha))
    err = pln1-pln0
    pln0 = pln1
    tt = tt+1
  }    ### EM-iteration is done.
  
  ln = pln1 - lambda*sum(log(alpha))
  index  = sort(theta, index.return=TRUE)$ix
  alpha0 = alpha[index]
  theta0 = theta[index]
  ## penalized likelihood value, 
  ## report subpopulation means in the increasing order.
  
  list("alpha"=alpha0, 
       "theta"=theta0,
       "loglik"=ln,
       "ploglik"=pln1,
       "iter.n" = tt)
  ## iter.n is the number of iteration used.
  ## useful to determine if the EM has converged.
}