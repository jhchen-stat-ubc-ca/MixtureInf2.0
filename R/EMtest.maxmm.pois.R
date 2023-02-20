#' EMtest.maxmm.pois
#'
#' @description This function computes the PMLE of parameters under the alternative model 
#'              for EM-test given a beta.i value.
#'              It is used in the EMstat.pois function.
EMtest.maxmm.pois <- function(count, freq, beta.i, alpha0, theta0, 
                              m0, n.init, n.iter, tol, max.iter) {	
  eta = rep(0, m0+1)
  eta[m0+1] = max(count)
  if(m0>1) { for(i in 2:m0) eta[i]=(theta0[i-1]+theta0[i])/2 }
  ### divide the parameter space of theta
  
  output=c()
  for (i in 0:n.init) {
    ###create initial values for EM-algorithm
    if(i < 1) {
      alpha = runif(m0)
      alpha = alpha/sum(alpha)
      alpha1= alpha*beta.i
      alpha2= alpha*(1-beta.i)
      theta1 = theta2 = rep(0, m0)
      for (l in 1:m0) {
        theta1[l]=runif(1, eta[l], eta[l+1])
        theta2[l]=runif(1, eta[l], eta[l+1]) 
      }
    } else {
      alpha1 = alpha0*beta.i;   alpha2 = alpha0*(1-beta.i)
      for (l in 1:m0) {
        theta1[l] = (4*theta0[l]+runif(1, eta[l], eta[l+1]))/5 
        theta2[l] = (4*theta0[l]+runif(1, eta[l], eta[l+1]))/5 }
      ## added extra initial with the null mixing proportion
    }
    
    ###run n.iter EM-iterations to find a winner.
    for (j in 1:n.iter) {
      pdf.part1 =apply(as.matrix(theta1, ncol=1), 1, dpois, x=count)
      pdf.part2 =apply(as.matrix(theta2, ncol=1), 1, dpois, x=count)
      pdf.sub1 =t(t(pdf.part1)*alpha1)+1e-50/m0
      pdf.sub2 =t(t(pdf.part2)*alpha2)+1e-50/m0
      pdf = apply(pdf.sub1, 1, sum) + apply(pdf.sub2, 1, sum)
      w1 = pdf.sub1/pdf
      w2 = pdf.sub2/pdf
      alpha = apply(freq*(w1+w2), 2, sum)/sum(freq)
      alpha1 = alpha*beta.i
      alpha2 = alpha*(1-beta.i)
      theta1=apply(freq*w1*count,2,sum)/apply(freq*w1,2,sum)
      theta2=apply(freq*w2*count,2,sum)/apply(freq*w2,2,sum)
      for(l in 1:m0) {
        theta1[l]=max(min(theta1[l], eta[l+1]), eta[l])
        theta2[l]=max(min(theta2[l], eta[l+1]), eta[l]) }
    }  ### enforce the restriction on the subpop means.
    ### completed the fixed number of iterations.
    
    pdf = dmix.pois(count, c(alpha1, alpha2), c(theta1, theta2))
    ln = sum(freq*log(pdf))   ## log likelihood value
    output=rbind(output, c(alpha1, alpha2, theta1, theta2, ln))
  }
  
  ## pick the winner among initial values.
  index=which.max(output[,(4*m0+1)])
  alpha1 = output[index,1:m0]
  alpha2 = output[index,(m0+1):(2*m0)]
  theta1 = output[index,(2*m0+1):(3*m0)]
  theta2 = output[index,(3*m0+1):(4*m0)]
  ln0=output[index,(4*m0+1)]
  
  ## real iteration with winning initial value starts.
  err=1; tt=0
  pdf.part1 =apply(as.matrix(theta1,ncol=1),1,dpois,x=count)
  pdf.part2 =apply(as.matrix(theta2,ncol=1),1,dpois,x=count)
  pdf.sub1 =t(t(pdf.part1)*alpha1)+1e-50/m0
  pdf.sub2 =t(t(pdf.part2)*alpha2)+1e-50/m0
  pdf=apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
  while (err>tol & tt< max.iter)  {
    w1 =pdf.sub1/pdf
    w2 =pdf.sub2/pdf
    alpha=apply(freq*w1+freq*w2,2,sum)/sum(freq)
    alpha1=alpha*beta.i
    alpha2=alpha*(1-beta.i)
    theta1=apply(freq*w1*count,2,sum)/apply(freq*w1,2,sum)
    theta2=apply(freq*w2*count,2,sum)/apply(freq*w2,2,sum)
    for(l in 1:m0) {
      theta1[l]=max(min(theta1[l],eta[l+1]),eta[l])
      theta2[l]=max(min(theta2[l],eta[l+1]),eta[l]) 
      ### pull back to intervals defined by eta.
    }
    pdf = dmix.pois(count, c(alpha1, alpha2), c(theta1, theta2))
    ln1 = sum(freq*log(pdf))
    err = ln1-ln0
    ln0 = ln1
    tt = tt+1
  }
  list("alpha1"=alpha1, 
       "alpha2"=alpha2, 
       "theta1"=theta1,
       "theta2"=theta2)
}