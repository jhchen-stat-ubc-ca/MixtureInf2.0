#' emstat.binom
#'
#' @description This function computes the EM-test statistics for a binomial mixture.
#' It is used in the emtest.binom function.
emstat.binom <- function(count, freq, outnull, size, m0, CC, 
                         n.init, max.iter, n.iter, tol, k) {
  bbeta=c()
  for(h in 1:m0)  {
    bbeta=rbind(cbind(bbeta, rep(0.1,3^{h-1})),
                cbind(bbeta, rep(0.3,3^{h-1})),
                cbind(bbeta, rep(0.5,3^{h-1})))
  }
  
  pln1=c();   para=c(); iter.n = 0
  ##For each beta, calculate the statistic m_n^{(k)}
  for(j in 1:(3^m0)) {
    beta.j = bbeta[j,]
    para0 = maxmm.binom(count, freq, size, m0, beta.j, outnull[[2]], 
                        n.init, max.iter, n.iter, tol)
    alpha = para0$alpha
    theta1= para0$theta1
    theta2= para0$theta2
    iter.n = iter.n + para0$iter.n
    
    ###Iteration starts from here###
    for(i in 1:k) {
      ###E-step of EM-algorithm
      alpha1 = alpha*beta.j	
      alpha2 = alpha*(1-beta.j)
      pdf.part1=apply(as.matrix(theta1,ncol=1),1,dbinom,x=count,size=size)
      pdf.part2=apply(as.matrix(theta2,ncol=1),1,dbinom,x=count,size=size)
      pdf.sub1 =t(t(pdf.part1)*alpha1)+1e-100/m0
      pdf.sub2 =t(t(pdf.part2)*alpha2)+1e-100/m0
      pdf = apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
      w1 =pdf.sub1/pdf
      w2 =pdf.sub2/pdf
      alpha = apply(freq*w1+freq*w2,2,sum)/sum(freq)
      
      ###M-step of EM-algorithm
      for(h in 1:m0) {
        temp1 = sum(freq*w1[,h]);  temp2 = sum(freq*w2[,h])
        if (temp1 < temp2)
          beta.j[h]=min((temp1+CC)/(temp1+temp2+CC),0.5)
        else
          beta.j[h]=max((temp1+CC)/(temp1+temp2+CC),0.5)
      }    ## penalty for mixing proportion
      alpha1= alpha*beta.j
      alpha2= alpha*(1-beta.j)
      ## split the mixing weights.
      theta1= apply(freq*w1*count,2,sum)/apply(freq*w1*size,2,sum)
      theta2= apply(freq*w2*count,2,sum)/apply(freq*w2*size,2,sum)
      ## compute the mean.
    }		
    ###Compute the penalized log-likelihood value and EM-test statistic
    pdf = dmix.binom(count, size, c(alpha1,alpha2),c(theta1,theta2))
    para=rbind(para,c(alpha1,alpha2,theta1,theta2))
    pln1[j]=sum(freq*log(pdf))+CC*sum(log(1-abs(1-2*beta.j)))
  }
  index = which.max(pln1)
  para1 = para[index,]
  emnk  = 2*(pln1[index]-outnull[[3]])
  list(emnk, para1[1:(2*m0)], 
       para1[(2*m0+1):(4*m0)], iter.n)
}