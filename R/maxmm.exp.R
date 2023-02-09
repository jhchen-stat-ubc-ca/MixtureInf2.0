#' maxmm.exp
#'
#' @description This function computes the PMLE of parameters under the alternative model given a beta value. 
#' @param x The input data that can be either a vector or a matrix with the 1st column being the observed values
#' and the 2nd column being the corresponding frequency.
#' @param bbeta The fixed mixing proportions.
#' @param theta0 The subpopulations mean fitted under the null hypothesis.
#' @param n.init The number of initial values chosen for the EM-algorithm.
#' @param n.iter The initial number of EM iterations.
#' @param tol The tolerance value for the convergence of the EM-algorithm.
#' @param max.iter The maximum number of iterations.
#'
#' @return
#' @export
#'
#' @examples
maxmm.exp <- function(x, bbeta, theta0, n.init, n.iter, tol, max.iter) {
  ### Divide the parameter space of subpop mean)
  m0 = length(bbeta)
  eta = rep(0,m0+1);  eta[1] = min(x); eta[m0+1] = max(x)
  if(m0>1) {
    for(i in 2:m0) eta[i]=(theta0[i-1]+theta0[i])/2
  }
  
  #initial values for EM-algorithm
  output=c()
  for (i in 1:n.init) {
    alpha = runif(m0);   alpha = alpha/sum(alpha)
    alpha1= alpha*bbeta; alpha2 = alpha*(1-bbeta)
    rate1 = rate2 = rep(0, m0)
    for (l in 1:m0) {
      rate1[l] = 1/runif(1,eta[l],eta[l+1])
      rate2[l] = 1/runif(1,eta[l],eta[l+1])
    }
    
    ### first run n.iter EM-iterations for each initial		
    for (j in 1:n.iter)  {
      pdf.part1 = apply(as.matrix(rate1,ncol=1),1, dexp, x=x)
      pdf.part2 = apply(as.matrix(rate2,ncol=1),1, dexp, x=x)
      pdf.sub1 = t(t(pdf.part1)*alpha1)+1e-100/m0
      pdf.sub2 = t(t(pdf.part2)*alpha2)+1e-100/m0
      pdf = apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
      w1 = pdf.sub1/pdf
      w2 = pdf.sub2/pdf
      alpha = apply(w1+w2,2,mean)
      alpha1= alpha*bbeta
      alpha2= alpha*(1-bbeta)
      theta1 = apply(w1*x,2,sum)/apply(w1,2,sum)
      theta2 = apply(w2*x,2,sum)/apply(w2,2,sum)
      for(l in 1:m0) {
        theta1[l]=1/max(min(theta1[l],eta[l+1]),eta[l])
        theta2[l]=max(min(theta2[l],eta[l+1]),eta[l])
      }
      rate1 = 1/theta1; rate2 = 1/theta2
    }
    
    pdf.part1= apply(as.matrix(rate1,ncol=1),1,dexp,x=x)
    pdf.part2= apply(as.matrix(rate2,ncol=1),1,dexp,x=x)
    pdf.sub1 = t(t(pdf.part1)*alpha1)+1e-100/m0
    pdf.sub2 = t(t(pdf.part2)*alpha2)+1e-100/m0
    pdf = apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
    ln  = sum(log(pdf))
    output= rbind(output,c(alpha, rate1, rate2, ln))
  }
  
  index = which.max(output[,(3*m0+1)])
  alpha = output[index,1:m0]
  rate1 = output[index,(m0+1):(2*m0)]
  rate2 = output[index,(2*m0+1):(3*m0)]
  ln0   = output[index,(3*m0+1)]
  err = 1
  tt  = 0
  pdf.part1 = apply(as.matrix(rate1,ncol=1),1,dexp,x=x)
  pdf.part2 = apply(as.matrix(rate2,ncol=1),1,dexp,x=x)
  pdf.sub1  = t(t(pdf.part1)*alpha1)+1e-100/m0
  pdf.sub2  = t(t(pdf.part2)*alpha2)+1e-100/m0
  pdf=apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
  
  while (err > tol*.05 & tt < max.iter) {
    ###EM-iteration with the initial value with the 
    # largest penalized log-likelihood
    w1 = pdf.sub1/pdf
    w2 = pdf.sub2/pdf
    alpha = apply(w1+w2,2,mean)
    alpha1= alpha*bbeta;  alpha2= alpha*(1-bbeta)
    theta1= apply(w1*x,2,sum)/apply(w1,2,sum)
    theta2= apply(w2*x,2,sum)/apply(w2,2,sum)
    for(l in 1:m0)  {
      theta1[l]=max(min(theta1[l],eta[l+1]),eta[l])
      theta2[l]=max(min(theta2[l],eta[l+1]),eta[l])
    }
    
    rate1 = 1/theta1; rate2 = 1/theta2
    pdf.part1= apply(as.matrix(rate1,ncol=1),1,dexp,x=x)
    pdf.part2= apply(as.matrix(rate2,ncol=1),1,dexp,x=x)
    pdf.sub1 = t(t(pdf.part1)*alpha1)+1e-100/m0
    pdf.sub2 = t(t(pdf.part2)*alpha2)+1e-100/m0
    pdf = apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
    ln1 = sum(log(pdf))
    err = ln1 - ln0
    ln0 = ln1
    tt  = tt + 1
  }
  list("alpha"=alpha, "theta1"=theta1, "theta2"=theta2)
}