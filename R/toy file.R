pmle.norm <-
  function(x,m0=1,lambda=0,inival=NULL,len=10,niter=50,tol=1e-6,rformat=FALSE)
  {
    if (is.data.frame(x))
    {	
      if (ncol(x)==2)
        x=as.matrix(x)
      if (ncol(x)==1 | ncol(x)>2)
        x=x[,1]
    }
    if (is.matrix(x))
    {
      xx=c()
      for (i in 1:nrow(x))
        xx=c(xx,rep(x[i,1],x[i,2]))
      x=xx
    }
    
    out=phi.norm(x,m0,lambda,inival,len,niter,tol)
    alpha=out[[1]]
    mean=out[[2]]
    var=out[[3]]
    loglik=out[[4]]
    ploglik=out[[5]]
    
    if (rformat==F)
    {
      alpha=rousignif(alpha)
      mean=rousignif(mean)
      var=rousignif(var)
      loglik=rousignif(loglik)
      ploglik=rousignif(ploglik)
    }
    
    if (m0>1)
    {
      list('PMLE of mixing proportions'=alpha,
           'PMLE of means'=mean,
           'PMLE of variances'=var,
           'log-likelihood'=loglik,
           'Penalized log-likelihood'=ploglik)
    }
    else 
    {
      list('MLE of mixing proportions'=alpha,
           'MLE of means'=mean,
           'MLE of variances'=var,
           'log-likelihood:'=loglik)
    }
  }

pmle.norm0 <-
  function(x,var,m0=1,lambda=0,inival=NULL,len=10,niter=50,tol=1e-6,rformat=FALSE)
  {
    if (is.data.frame(x))
    {	
      if (ncol(x)==2)
        x=as.matrix(x)
      if (ncol(x)==1 | ncol(x)>2)
        x=x[,1]
    }
    if (is.matrix(x))
    {
      xx=c()
      for (i in 1:nrow(x))
        xx=c(xx,rep(x[i,1],x[i,2]))
      x=xx
    }
    
    x=x/sqrt(var)
    n=length(x)
    out=phi0.norm(x,m0,lambda,inival,len,niter,tol)
    alpha=out$alpha
    theta=out$theta*sqrt(var)
    loglik=out$loglik-n/2*log(var)
    ploglik=out$ploglik-n/2*log(var)
    
    if (rformat==F)
    {
      alpha=rousignif(alpha)
      theta=rousignif(theta)
      loglik=rousignif(loglik)
      ploglik=rousignif(ploglik)
    }
    
    if (lambda==0 | m0==1)
      list('MLE of mixing proportions:'=alpha,
           'MLE of component parameters:'=theta,
           'log-likelihood:'=loglik)
    else
      list('PMLE of mixing proportions:'=alpha,
           'PMLE of component parameters:'=theta,
           'log-likelihood:'=loglik,
           'Penalized log-likelihood:'=ploglik)
  }
