#' pmle.norm0.sub
#'
#' @param x Input data as a plain vector.
#' @param m0 The order of the normal mixture to be fitted.
#' @param lambda The size of penalty function on mixing proportions.
#' @param init.val The initial values for the EM-algorithm in matrix form. Each row = c(alpha, mu).
#' @param n.init The number of initial values for the EM-algorithm.
#' @param n.iter The number of EM-iterations for each initial value.
#' @param max.iter The subsequent number of EM-iterations allowed.
#' @param tol The tolerance value for the convergence of the EM-algorithm.
#'
#' @return
#' @export
#'
#' @examples
pmle.norm0.sub <- function(x, m0, lambda, 
                           init.val, n.init, n.iter, max.iter, tol)
{
  nn = length(x)
  output=c()
  for(i in 1:n.init)  {
    alpha = init.val[i, 1:m0]
    theta = init.val[i, (m0+1):(2*m0)]
    for (j in 1:n.iter)  {   
      ### run n.iter EM-iterations to find the best init.val
      pdf.component = apply(as.matrix(theta,ncol=1),1,dnorm,x=x,sd=1)
      pdf.component = t(t(pdf.component)*alpha)+1e-100/m0
      pdf = apply(pdf.component, 1, sum)
      ww = pdf.component/pdf
      alpha = (apply(ww, 2, sum) + lambda)/(nn + m0*lambda)
      theta = apply(ww*x,2,sum)/apply(ww,2,sum)
    }
    pdf.component = apply(as.matrix(theta,ncol=1),1,dnorm,x=x,sd=1)
    pdf.component=t(t(pdf.component)*alpha)+1e-100/m0
    pdf = apply(pdf.component,1,sum)
    pln = sum(log(pdf)) + lambda*sum(log(alpha))
    output = rbind(output,c(alpha,theta,pln))
  }
  index = which.max(output[,(2*m0+1)])
  alpha = output[index, 1:m0]
  theta = output[index,(m0+1):(2*m0)]
  pln0 = output[index,(2*m0+1)]
  err=1
  t=0
  pdf.component = apply(as.matrix(theta,ncol=1),1,dnorm,x=x,sd=1)
  pdf.component = t(t(pdf.component)*alpha)+1e-100/m0
  pdf = apply(pdf.component,1,sum)
  while(err > tol & t< max.iter)
    ### EM-iteration with the best init.val
  {
    ww = pdf.component/pdf
    alpha = (apply(ww,2,sum)+lambda)/(nn + m0*lambda)
    theta = apply(ww*x,2,sum)/apply(ww,2,sum)
    pdf.component = apply(as.matrix(theta,ncol=1),1,dnorm,x=x,sd=1)
    pdf.component = t(t(pdf.component)*alpha)+1e-100/m0
    pdf = apply(pdf.component,1,sum)
    pln1 = sum(log(pdf))+lambda*sum(log(alpha))
    err = pln1-pln0
    pln0 = pln1
    t = t+1
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
       "iter.n" = t)
}