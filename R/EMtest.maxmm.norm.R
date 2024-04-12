#' EMtest.maxmm.norm
#'
#' @description This function computes the mixture's PMLE with specific structure. 
#' It is used in the emtest.norm function.
#' @export
EMtest.maxmm.norm <- 
  function(xx, beta.i, m0, an, para0, n.init, n.iter, tol) {
    alpha0 = para0[[1]][1,];  mu0 = para0[[1]][2,]; sigma0 = para0[[1]][3,]
    output=c()
    
    ### divide the real line
    eta=rep(0, (m0+1)); eta[1] = min(xx); eta[m0+1] = max(xx)
    if (m0 > 1) eta[2:m0] = (mu0[1:(m0-1)] + mu0[2:m0])/2
    
    ###initial values for specific EM-algorithm
    for (i in 0:n.init) {
      alpha = mu = sigma=c()
      if(i > 0) {
        for (j in 1:m0) {	
          alpha = c(alpha, alpha0[j]*beta.i[j], alpha0[j]*(1-beta.i[j]))
          mu = c(mu, runif(2, eta[j], eta[j+1]))
          sigma = c(sigma, runif(2, 0.25*sigma0[j], 2*sigma0[j])) }
        para = c(alpha, mu, sigma)
      } else {
        for (j in 1:m0) {	
          alpha = c(alpha, alpha0[j]*beta.i[j], alpha0[j]*(1-beta.i[j]))
          mu = c(mu, (4*mu0[j]+runif(2, eta[j], eta[j+1]))/5)
          sigma = c(sigma, runif(2, 0.9*sigma0[j], 1.1*sigma0[j])) 
        }
      }
      para = c(alpha, mu, sigma)
      ## added an extra one with the null subpop parameter values. 
      
      ###run n.iter specfic EM-iterations first
      for (j in 1:n.iter) {
        outpara = single.iter.norm(xx, m0, para, sigma0, beta.i, eta, an)
        para = outpara[1:(6*m0)] }
      output = rbind(output, outpara)	
    }  ## loop i completes here.
    
    ### choose the most performing one.
    index = which.max(output[,(6*m0+1)])
    para = output[index,1:(6*m0)]
    pln0 = output[index,(6*m0+1)]
    err = 1
    tt = 0
    
    while (err > tol & tt < 2000) {
      ### EM-iteration with the initial value with 
      ###  the largest penalized log-likelihood
      outpara = single.iter.norm(xx, m0, para, sigma0, beta.i, eta, an)
      para = outpara[1:(6*m0)]
      pln1 = outpara[6*m0+1]
      err = pln1-pln0
      pln0 = pln1
      tt = tt+1
    }
    para
  }


single.iter.norm <-
  function(xx, m0, para, sigma0, beta.i, eta, an) {
    ## replacement of iter2.norm
    alpha = para[1:(2*m0)]
    mu = para[(2*m0+1):(4*m0)]
    sigma = para[(4*m0+1):(6*m0)]; sig = sigma^0.5
    
    ###E-step of EM-algorithm
    pdf.sub = xx.sq = c()
    for (j in 1:(2*m0)) {	
      pdf.sub=cbind(pdf.sub, dnorm(xx,mu[j],sig[j])*alpha[j])
      xx.sq = cbind(xx.sq,(xx-mu[j])^2)
    }
    pdf = apply(pdf.sub,1,sum)+1e-50
    ww = pdf.sub/pdf
    
    ###M-step of EM-algorithm
    ww.tot = apply(ww,2,sum)
    alpha = ww.tot/length(xx)
    mu = apply(ww*xx,2,sum)/ww.tot
    sigma = (apply(ww*xx.sq,2,sum)+2*an*sigma0)/(ww.tot+2*an)
    
    alpha1=c(); mu1=c()
    for (j in 1:m0) {
      alpha1 = c(alpha1, (alpha[2*j-1]+alpha[2*j])*c(beta.i[j], 1-beta.i[j]))
      mu1 = c(mu1, min(max(mu[2*j-1], eta[j]), eta[j+1]),
              min(max(mu[2*j], eta[j]), eta[j+1])) }
    sig1 = sigma^.5
    ### confined within eta intervals.
    
    ###Compute the log-likelihood value
    ln = sum(log(dmix.norm(xx, alpha1, mu1, sig1) + 1e-50))
    temp = array(rbind(sigma0, sigma0))/sigma  
    ## vector with repeated entries
    pln = ln -an*sum(temp - log(temp)-1)  
    ## compare again pn before finalize.
    ## -an*( sigma0/sigma + log(sigma/sigma0) -1)
    output = c(alpha1, mu1, sigma, pln)
  }