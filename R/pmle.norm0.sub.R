#' pmle.norm0.sub
#'
#' @description It is used in the pmle.norm0 function, It does the actual computing for the EM-algorithm. 
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