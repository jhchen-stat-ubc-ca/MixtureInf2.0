#' EMtest.maxmm.norm0
#'
#' @param xx The input data that can be either a vector or a matrix with the 1st column being the observed values
#'          and the 2nd column being the corresponding frequency. It is used in the emstat.norm0 function.
EMtest.maxmm.norm0 <-
  function(xx, bbeta, mle0, m0, n.init, n.iter, tol) {
    eta = rep(0, m0+1);  eta[1] = min(xx);  eta[m0+1] = max(xx)
    theta0 = mle0[2,]
    if(m0>1) {
      for(i in 2:m0) eta[i]=(theta0[i-1]+theta0[i])/2
    }
    # divide space of mean into
    #  [eta_1, eta_2], [eta_2, eta_3] an so on.
    
    output=c()
    theta1 = theta2 = rep(0, m0)
    for (i in 0:n.init)  {
      ###initial values for EM-iteration
      if(i > 0) {
        alpha=runif(m0);     alpha=alpha/sum(alpha)
        alpha1=alpha*bbeta;   alpha2=alpha*(1-bbeta)
        for (l in 1:m0) {
          theta1[l]=runif(1, eta[l], eta[l+1])
          theta2[l]=runif(1, eta[l], eta[l+1]) }
      } else {
        alpha1=mle0[1,]*bbeta;   alpha2=mle0[1,]*(1-bbeta)
        for (l in 1:m0) {
          theta1[l] = (4*mle0[2, l]+runif(1, eta[l], eta[l+1]))/5
          theta2[l] = (4*mle0[2, l]+runif(1, eta[l], eta[l+1]))/5 }
        ## added extra one with the null mixing proportion 
      }
      
      for (j in 1:n.iter) {
        #pdf.sub1 = dmix.norm(x, alpha1, theta1, rep(1, m0))
        #pdf.sub2 = dmix.norm(x, alpha2, theta2, rep(1, m0))
        pdf.part1=apply(as.matrix(theta1,ncol=1),1,dnorm,x=xx)
        pdf.part2=apply(as.matrix(theta2,ncol=1),1,dnorm,x=xx)
        ## density for each subpopulation (2*m0)
        pdf.sub1=t(t(pdf.part1)*alpha1)+1e-100/m0
        pdf.sub2=t(t(pdf.part2)*alpha2)+1e-100/m0
        ## mixing weights applied
        pdf=apply(pdf.sub1+pdf.sub2,1,sum)
        w1 =pdf.sub1/pdf
        w2 =pdf.sub2/pdf
        alpha = apply(w1+w2,2,mean)
        alpha1=alpha*bbeta
        alpha2=alpha*(1-bbeta)
        theta1=apply(w1*xx,2,sum)/apply(w1,2,sum)
        theta2=apply(w2*xx,2,sum)/apply(w2,2,sum)
        for(l in 1:m0) {
          theta1[l]=max(min(theta1[l],eta[l+1]),eta[l])
          theta2[l]=max(min(theta2[l],eta[l+1]),eta[l])
        }  ### enforce the range of subpop means
      }  ## end of loop j:  fix number of iteration.
      
      pdf = dmix.norm(xx, c(alpha1, alpha2), c(theta1, theta2), 
                      rep(1, 2*m0))
      ln = sum(log(pdf))
      output=rbind(output, c(alpha1, alpha2, theta1, theta2,ln))
    }
    
    index = which.max(output[,(4*m0+1)])
    ## pick up the winner among all initial values.
    alpha1=output[index,1:m0]
    alpha2=output[index,(m0+1):(2*m0)]
    theta1=output[index,(2*m0+1):(3*m0)]
    theta2=output[index,(3*m0+1):(4*m0)]
    ln0=output[index,(4*m0+1)]
    
    err=1; tt=0
    while (err > tol & tt < 2000) {
      ### continue with the initial value with 
      ### the largest penalized log-likelihood
      pdf.part1 = apply(as.matrix(theta1,ncol=1),1,dnorm,x=xx)
      pdf.part2 = apply(as.matrix(theta2,ncol=1),1,dnorm,x=xx)
      pdf.sub1 = t(t(pdf.part1)*alpha1)+1e-100/m0
      pdf.sub2 = t(t(pdf.part2)*alpha2)+1e-100/m0
      pdf = apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
      ln1 = sum(log(pdf))
      
      w1=pdf.sub1/pdf;   w2=pdf.sub2/pdf
      alpha =apply(w1+w2,2,mean)
      alpha1=alpha*bbeta
      alpha2=alpha*(1-bbeta)
      theta1=apply(w1*xx,2,sum)/apply(w1,2,sum)
      theta2=apply(w2*xx,2,sum)/apply(w2,2,sum)
      for(l in 1:m0) {
        theta1[l]=max(min(theta1[l],eta[l+1]),eta[l])
        theta2[l]=max(min(theta2[l],eta[l+1]),eta[l])
      }
      err = ln1-ln0
      ln0 = ln1
      tt=tt+1
    }
    list("alpha"=alpha,"theta1"=theta1,"theta2"=theta2)
  }