#' pmle.norm0.sub
#'
#' @description It is used in the pmle.norm0 function, It does the actual computing for the EM-algorithm.
#' @param x The input data that can either be a vector or a matrix with the 1st column being the observed values
#' and the 2nd column being the corresponding frequencies.
#' @param m0 The order of the finite mixture model.
#' @param lambda The size of the penalty function of the mixing proportions.
#' @param init.val The initial values chosen for the EM-algorithm, a 3m0-dimension vector including m0 mixing proportions, 
#' m0 component means and m0 component variances, or a matrix with 3m0 columns, default value: inival = NULL. (if not provided, random initial values are used.)
#' @param n.init The	number of initial values for the EM-algorithm.
#' @param n.iter The number of EM iterations for all initial values.
#' @param max.iter  Maximum amount of EM iterations, it stops at 5000.
#' @param tol The tolerance value for the convergence of the EM-algorithm, default value: tol = 1e-8.
#'  
#' @export 
pmle.norm0.sub <- function(x, m0, lambda, 
                           init.val, n.init, n.iter, max.iter, tol) {
  nn = length(x)
  output=c()
  for(i in 1:n.init)  {
    alpha = init.val[[i]][1,]
    theta = init.val[[i]][2,]
    for (j in 1:n.iter)  {   
      ### run n.iter EM-iterations to find the best init.val
      pdf.sub = apply(as.matrix(theta,ncol=1),1,dnorm,x=x,sd=1)
      pdf.sub = t(t(pdf.sub)*alpha)+1e-100
      pdf = apply(pdf.sub, 1, sum)
      ww = pdf.sub/pdf
      alpha = (apply(ww, 2, sum) + lambda)/(nn + m0*lambda)
      theta = apply(ww*x,2,sum)/apply(ww,2,sum)
    }
    pdf.sub = apply(as.matrix(theta,ncol=1),1,dnorm,x=x,sd=1)
    pdf.sub=t(t(pdf.sub)*alpha)+1e-100
    pdf = apply(pdf.sub,1,sum)
    pln = sum(log(pdf)) + lambda*sum(log(alpha))
    output = rbind(output,c(alpha,theta,pln))
  }
  index = which.max(output[,(2*m0+1)])
  alpha = output[index, 1:m0]
  theta = output[index,(m0+1):(2*m0)]
  pln0 = output[index,(2*m0+1)]
  err=1
  tt = 0
  pdf.sub = apply(as.matrix(theta,ncol=1),1,dnorm,x=x,sd=1)
  pdf.sub = t(t(pdf.sub)*alpha)+1e-100
  pdf = apply(pdf.sub,1,sum)
  while(err > tol & tt < max.iter)
    ### EM-iteration with the best init.val
  {
    ww = pdf.sub/pdf
    alpha = (apply(ww,2,sum)+lambda)/(nn + m0*lambda)
    theta = apply(ww*x,2,sum)/apply(ww,2,sum)
    pdf.sub = apply(as.matrix(theta,ncol=1),1,dnorm,x=x,sd=1)
    pdf.sub = t(t(pdf.sub)*alpha)+1e-100
    pdf = apply(pdf.sub,1,sum)
    pln1 = sum(log(pdf))+lambda*sum(log(alpha))
    err = pln1-pln0
    pln0 = pln1
    tt = tt + 1
  }
  ln = pln1 - lambda*sum(log(alpha))
  index = sort(theta,index.return=TRUE)$ix
  alpha0 = alpha[index]
  theta0 = theta[index]
  ### report outcome with increasing subpopulation means
  
  list("alpha"= alpha0,
       "theta"= theta0,
       "loglik"= ln,
       "ploglik"= pln1,
       "iter.n" = tt)
}