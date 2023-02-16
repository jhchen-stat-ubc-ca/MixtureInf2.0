#' pmle.exp.sub
#'
#' @description This function does the real work of computing the PMLE of the mixing distribution under the exponential mixture model. 
#' It is used in the pmle.exp function.
pmle.exp.sub <- function(xx, m0, lambda, init.val, n.init, 
                         n.iter, max.iter, tol) {
  n = length(xx)
  pdf.sub = w = matrix(0, m0, n)
  if (!is.null(init.val)) n.init = 1
  
  output=c()
  for (i in 1:n.init) {
    if (is.null(init.val)) {
      ## generate random initial values if not provided.
      theta = sort(sample(xx, m0))
      tmp = (theta[-1] + theta[-m0])/2
      alpha = rep(1, m0)
      for(ii in 1:(m0-1)) alpha[ii] = sum(xx<tmp[ii])/n
      alpha = alpha - c(0, alpha[-m0])
    }	else {	### if initial mixing distribution is provided
      alpha = init.val[1,]
      theta = init.val[2,]
    }
    
    for (j in 1:n.iter) {
      ###run n.iter EM-iterations first
      for(j in 1:m0) pdf.sub[j,] = dexp(xx, 1/theta[j])*alpha[j]
      pdf.mixture = apply(pdf.sub,2,sum) + 1e-100
      for(j in 1:m0) {
        w[j,] = pdf.sub[j,]/pdf.mixture 
        alpha[j] = (sum(w[j,])+lambda)/(n+m0*lambda)
        theta[j] = sum(w[j,]*xx)/sum(w[j,])
      }
    }
    
    for(j in 1:m0) pdf.sub[j,] = dexp(xx, 1/theta[j])*alpha[j]
    pdf.mixture = apply(pdf.sub,2,sum) +1e-100
    pln = sum(log(pdf.mixture))+lambda*sum(log(alpha))
    output = rbind(output, c(alpha, theta, pln))
  }
  ### pick the winner
  index = which.max(output[,(2*m0+1)])
  alpha = output[index, 1:m0]
  theta = output[index,(m0+1):(2*m0)]
  pln0  = output[index,(2*m0+1)]
  increment = 1
  tt = 0
  for(j in 1:m0) pdf.sub[j,] = dexp(xx, 1/theta[j])*alpha[j]
  pdf.mixture = apply(pdf.sub, 2, sum)
  
  ### restart EM-iteration with the winner.
  while(increment > tol & tt < max.iter) {
    for(j in 1:m0) {
      w[j,] = pdf.sub[j,]/pdf.mixture 
      alpha[j] = (sum(w[j,])+lambda)/(n+m0*lambda)
      theta[j] = sum(w[j,]*xx)/sum(w[j,])
    }
    for(j in 1:m0) pdf.sub[j,] = dexp(xx, 1/theta[j])*alpha[j]
    pdf.mixture = apply(pdf.sub, 2, sum)
    pln1 = sum(log(pdf.mixture)) + lambda*sum(log(alpha))
    increment = pln1-pln0 
    pln0 = pln1
    tt = tt + 1
  }
  ln = pln1-lambda*sum(log(alpha))
  index = sort(theta,index.return=TRUE)$ix
  alpha0 = alpha[index]
  theta0 = theta[index]
  
  list("alpha"=alpha0,
       "theta"=theta0,
       "loglik"=ln,
       "ploglik"=pln1,
       "n.iter" = tt)
} 