#' maxmm.binom
#'
#' @description This function is meant for em-test. It gets the PMLE under a restricted alternative given a beta.j value.
#' It is used in the emtest.binom function.
maxmm.binom <- function(count, freq, size, m0, beta.j, theta0, 
                        n.init, max.iter, n.iter, tol)   {
  ### the dividing points of the parameter space of theta)
  eta=rep(0, m0+1)
  eta[1]=0
  eta[m0+1]=1
  if(m0 > 1) { for(i in 2:m0) eta[i]=(theta0[i-1]+theta0[i])/2 }
  
  output=c()
  theta1 = theta2 = rep(0,m0)
  for (i in 1:n.init) {
    ###initial values for EM-algorithm
    alpha = runif(m0);    alpha = alpha/sum(alpha)
    alpha1 = alpha*beta.j;   alpha2 = alpha*(1-beta.j)
    for (l in 1:m0) {
      theta1[l] = runif(1,eta[l],eta[l+1])
      theta2[l] = runif(1,eta[l],eta[l+1])
    }
    
    for (j in 1:n.iter) ###run n.iter EM-iterations first
    {
      pdf.part1=apply(as.matrix(theta1,ncol=1),1,dbinom,x=count,size=size)
      pdf.part2=apply(as.matrix(theta2,ncol=1),1,dbinom,x=count,size=size)
      pdf.sub1=t(t(pdf.part1)*alpha1)+1e-100/m0
      pdf.sub2=t(t(pdf.part2)*alpha2)+1e-100/m0
      pdf = apply(pdf.sub1,1,sum) + apply(pdf.sub2,1,sum)
      w1=pdf.sub1/pdf;   	 w2=pdf.sub2/pdf
      alpha=apply(freq*w1+freq*w2,2,sum)/sum(freq)
      alpha1= alpha*beta.j; alpha2=alpha*(1-beta.j)
      theta1 = apply(freq*w1*count,2,sum)/apply(freq*w1*size,2,sum)
      theta2 = apply(freq*w2*count,2,sum)/apply(freq*w2*size,2,sum)
      for(l in 1:m0) {
        theta1[l]=max(min(theta1[l],eta[l+1]), eta[l])
        theta2[l]=max(min(theta2[l],eta[l+1]), eta[l])
      }   ## apply the restriction for em-test null fit.
    }
    pdf = dmix.binom(count, size, c(alpha1, alpha2), c(theta1, theta2))
    ln  = sum(freq*log(pdf+1e-100))
    output = rbind(output,c(alpha1,alpha2,theta1,theta2,ln))
  }
  
  index=which.max(output[,(4*m0+1)])
  alpha1 = output[index,1:(m0)]
  alpha2 = output[index,(m0+1):(2*m0)]
  theta1 = output[index,(2*m0+1):(3*m0)]
  theta2 = output[index,(3*m0+1):(4*m0)]
  ln0 = output[index,(4*m0+1)]
  ## identify the best initial 
  
  pdf.part1=apply(as.matrix(theta1,ncol=1),1,dbinom,x=count,size=size)
  pdf.part2=apply(as.matrix(theta2,ncol=1),1,dbinom,x=count,size=size)
  pdf.sub1=t(t(pdf.part1)*alpha1)+1e-100/m0
  pdf.sub2=t(t(pdf.part2)*alpha2)+1e-100/m0
  pdf = apply(pdf.sub1,1,sum) + apply(pdf.sub2,1,sum)
  
  err = 1; tt = 0
  while (err > tol & tt < max.iter) {
    ###EM-iteration the best initial value.
    w1 = pdf.sub1/pdf;            w2 = pdf.sub2/pdf
    alpha = apply(freq*w1+freq*w2,2,sum)/sum(freq)
    alpha1 = alpha*beta.j;        alpha2=alpha*(1-beta.j)
    theta1 = apply(freq*w1*count,2,sum)/apply(freq*w1*size,2,sum)
    theta2 = apply(freq*w2*count,2,sum)/apply(freq*w2*size,2,sum)
    for(l in 1:m0) {
      theta1[l]=max(min(theta1[l],eta[l+1]),eta[l])
      theta2[l]=max(min(theta2[l],eta[l+1]),eta[l])
    }   ### the restriction stays.
    pdf.part1=apply(as.matrix(theta1,ncol=1),1,dbinom,x=count,size=size)
    pdf.part2=apply(as.matrix(theta2,ncol=1),1,dbinom,x=count,size=size)
    pdf.sub1=t(t(pdf.part1)*alpha1)+1e-100/m0
    pdf.sub2=t(t(pdf.part2)*alpha2)+1e-100/m0
    pdf = apply(pdf.sub1,1,sum) + apply(pdf.sub2,1,sum)
    ln1 = sum(freq*log(pdf))
    err = ln1 -ln0
    ln0 = ln1
    tt = tt + 1
  }
  list("alpha"=alpha,"theta1"=theta1,"theta2"=theta2, "iter.n"=tt)
}