#' pmle.binom.sub
#'
#' @description A sub function for pmle.binom, does the actual work of PMLE binomial mixture.
#' It is used in the pmle.binom function.
#' @export
pmle.binom.sub <- function(count, freq, size, m0, lambda,
                           init.val, n.init, max.iter, n.iter, tol)  {
  if(m0 > 1)	{
    if (is.matrix(init.val))  n.init = nrow(init.val)
    if (is.vector(init.val)) { 
      init.val = t(as.matrix(init.val)); n.init = 1}
    output=c()
    for(i in 1:n.init) {
      if (is.null(init.val))   {
        alpha = runif(m0,0,1)
        alpha = alpha/sum(alpha)
        theta = sort(runif(m0,0,1))
      }
      else  {
        alpha=init.val[i,1:m0]
        alpha=alpha/sum(alpha)
        theta=sort(init.val[i,(m0+1):(2*m0)])
      }
      
      for (j in 1:n.iter) ###run n.iter EM-iterations first
      {
        pdf.sub =t(t(apply(as.matrix(theta,ncol=1),1,dbinom,x=count,size=size))*alpha)+1e-100/m0
        pdf=apply(pdf.sub,1,sum)
        ww = pdf.sub/pdf
        alpha=(apply(freq*ww,2,sum)+lambda)/(sum(freq)+m0*lambda)
        theta=apply(freq*ww*count, 2, sum)/apply(freq*ww*size, 2, sum)
      }
      pdf = dmix.binom(count, size, alpha, theta)
      pln = sum(freq*log(pdf))+lambda*sum(log(alpha))
      output = rbind(output, c(alpha, theta, pln))
    }
    index = which.max(output[,(2*m0+1)])
    alpha = output[index,1:m0]
    theta = output[index,(m0+1):(2*m0)]
    pln0  = output[index,(2*m0+1)]
    ## identified the best initial value.
    
    err = 1;   tt  = 0
    pdf.sub = t(t(apply(as.matrix(theta,ncol=1),1,dbinom,x=count,size=size))*alpha)+1e-100/m0
    pdf = apply(pdf.sub,1,sum)
    while(err > tol & tt < max.iter)  {
      ww = pdf.sub/pdf
      alpha = (apply(freq*ww,2,sum)+lambda)/(sum(freq)+m0*lambda)
      theta = apply(freq*ww*count,2,sum)/apply(freq*ww*size,2,sum)
      pdf.sub =t(t(apply(as.matrix(theta,ncol=1),1,dbinom,x=count,size=size))*alpha)+1e-100/m0
      pdf   = apply(pdf.sub,1,sum)
      pln1  = sum(freq*log(pdf))+lambda*sum(log(alpha))
      err   = pln1 - pln0
      pln0  = pln1
      tt = tt+1
    }
    
    ln = pln1 - lambda*sum(log(alpha))
    index  = sort(theta,index.return=TRUE)$ix
    alpha0 = alpha[index]
    theta0 = theta[index]
    ## report output with theta sorted.
    
    list("alpha"=alpha0,"theta"=theta0,
         "loglik"=ln,"ploglik"=pln1, "iter.n"=tt)
  }
  else  {  
    theta0 = sum(freq*count)/sum(freq*size)
    list("alpha"=1, "theta"=theta0,
         "loglik" =sum(freq*log(dbinom(count,size,theta0))),
         "ploglik"=sum(freq*log(dbinom(count,size,theta0))),
         "iter.n" =0)   }
}